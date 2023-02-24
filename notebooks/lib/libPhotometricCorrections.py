# libPhotometricCorrections.py
#
# Author          : Sylvie Dagoret-Campagne
# Affiliaton      : IJCLab/IN2P3/CNRS
# Creation Date   : 2023/02/23
# Last update     : 2023/02/22
#
# A python tool to calculate Photometric Correction
# 
#
#
import os
import sys
from pathlib import Path

import atmosphtransmemullsst
from atmosphtransmemullsst.simpleatmospherictransparencyemulator import SimpleAtmEmulator

from scipy import interpolate
import numpy as np

from rubin_sim.phot_utils import Bandpass, Sed

#README.md        darksky.dat      filter_r.dat     hardware_g.dat   hardware_y.dat   lens3.dat        total_g.dat      total_y.dat
#README_SOURCE.md detector.dat     filter_u.dat     hardware_i.dat   hardware_z.dat   m1.dat           total_i.dat      total_z.dat
#atmos_10.dat     filter_g.dat     filter_y.dat     hardware_r.dat   lens1.dat        m2.dat           total_r.dat      version_info
#atmos_std.dat    filter_i.dat     filter_z.dat     hardware_u.dat   lens2.dat        m3.dat           total_u.dat
hardware_filenames = ["hardware_u.dat","hardware_g.dat","hardware_r.dat","hardware_i.dat","hardware_z.dat","hardware_y.dat"] 
filter_filenames = ["filter_u.dat","filter_g.dat","filter_r.dat","filter_i.dat","filter_z.dat","filter_y.dat" ]
total_filenames = ["total_u.dat","total_g.dat","total_r.dat","total_i.dat","total_z.dat","total_y.dat" ]
filter_tagnames = ["u","g","r","i","z","y"]
Filter_tagnames = ["U","G","R","I","Z","Y"]
filtercolor_tagnames = ["u-g","g-r","r-i","i-z","z-y"]
Filtercolor_tagnames = ["U-G","G-R","R-I","I-Y","Z-Y"]
filter_color = ["b","g","r","orange","grey","k"]
NFILT=len(filter_filenames)

WLMIN=300.
WLMAX=1100.
WLBIN=1.
NWLBIN=int((WLMAX-WLMIN)/WLBIN)
WL=np.linspace(WLMIN,WLMAX,NWLBIN)

#FILTERWL: precalculated array containing center, boundaries and width of each filter.
#index 0 : minimum wavelength of filter border
#index 1 : minimum wavelength of filter border
#index 2 : center wavelength of filter
#index 3 : filter width


FILTERWL = np.array([[ 324.03003755,  402.12765957,  363.59690349,   78.09762203],
       [ 392.11514393,  561.32665832,  473.54069923,  169.21151439],
       [ 542.3028786 ,  700.50062578,  619.49926767,  158.19774718],
       [ 681.47684606,  827.65957447,  752.01084117,  146.18272841],
       [ 808.63579474,  932.79098874,  868.488419  ,  124.15519399],
       [ 914.76846058, 1044.93116395,  969.10570859,  130.16270338]])


def fII0(self,wl,s):
        return np.trapz(s/wl,wl)
      
def fII1(wl,phi,wlb):
    return np.trapz(phi*(wl-wlb),wl)
  

print("libPhotometricCorrections.py :: Use atmosphtransmemullsst.__path__[0],'../data/simplegrid as the path to data")
data_path = os.path.join(atmosphtransmemullsst.__path__[0],'../data/simplegrid')
print(f"libPhotometricCorrections :: data_path = {data_path}")


       
# The emulator as a global variable
emul_atm = SimpleAtmEmulator(data_path)


class PhotometricCorrections:
  def __init__(self,am0=1.2,pwv0=5.0,oz0=300,ncomp0=1,tau0=0.04,beta0=-1):
        """
        Constructor
        """
        global emul_atm,sed
        self.WL = emul_atm.GetWL()
        
        # standard atmosphere parameters
        self.am0 = am0
        self.pwv0 = pwv0
        self.oz0 = oz0
        self.ncomp0 = ncomp0
        self.tau0 = tau0
        self.beta0 = beta0
        
        # standard atmosphere
        if ncomp0 == 1:
          taus = [tau0]
          betas = [beta0]    
          self.atm_std = emul_atm.GetAllTransparencies(self.WL,am0,pwv0,oz0,ncomp=ncomp0, taus=taus, betas=betas, flagAerosols=True)
        else:
          self.atm_std = emul_atm.GetAllTransparencies(self.WL,am0,pwv0,oz0)
          
        # instrumental filter
        self.bandpass_inst = {} 
        path_rubin_sim_throughput=os.path.join(os.getenv("HOME"),"rubin_sim_data/throughputs/baseline")
        for index,filename in enumerate(hardware_filenames):
          fullfilename=os.path.join(path_rubin_sim_throughput,filename)
          arr= np.loadtxt(fullfilename)
          # interpolate  filter transmission
          ff = interpolate.interp1d(x=arr[:,0], y=arr[:,1],fill_value="extrapolate")
          fname = filter_tagnames[index]
          self.bandpass_inst[fname] = Bandpass(wavelen=self.WL,sb=ff(self.WL))
          
        # total filter (instrumental x atmosphere)  
        self.bandpass_total_std = {} 
        for index,f in enumerate(filter_tagnames):
          self.bandpass_total_std[f] = Bandpass(wavelen=self.WL,sb=self.bandpass_inst[f].sb * self.atm_std)
         
        # Normalized response 
        self.phiArray_std, _ = Sed().setup_phi_array([self.bandpass_total_std[f] for f in filter_tagnames])
        
        # Integrals IIb0(std) and IIb1(std)
        self.all_II0_std = {}
        self.all_II1_std = {}
        for index,f in enumerate(filter_tagnames):
    
          the_II0 = self.fII0(self.bandpass_total_std[f].wavelen,self.bandpass_total_std[f].sb)
          self.all_II0_std[f] = the_II0
          the_II1 = self.fII1(self.WL,self.phiArray_std[index,:],FILTERWL[index,2])
          self.all_II1_std[f] = the_II1
          
          
        # Non standard calculations will be calculated later after initialisation
        self.am = 1.2
        self.pwv = 0
        self.oz = 0
        self.ncomp = 1
        self.tau = 0.04
        self.beta = -1
        
        self.atm_nonstd = None
        self.bandpass_total_nonstd = None
        self.phiArray_nonstd = None
        self.all_II0_nonstd = None
        self.all_II1_nonstd = None
        self.all_II0ratio_nonstd = None
        self.all_II1sub_nonstd = None
        
        self.coll_atm_nonstd = None
        self.coll_bandpass_total_nonstd = None
        self.coll_phiArray_nonstd = None
        self.coll_all_II0_nonstd = None
        self.coll_all_II1_nonstd = None
        self.coll_all_II0ratio_nonstd = None
        self.coll_all_II1sub_nonstd = None
        
        
         
          
  def fII0(self,wl,s):
        return np.trapz(s/wl,wl)
      
  def fII1(self,wl,phi,wlb):
        return np.trapz(phi*(wl-wlb),wl)
      
  def CalculateObs(self,am=1.2,pwv=5.0,oz=300,ncomp=1,tau=0.04,beta=-1):
        """
        """
        self.am = am
        self.pwv = pwv 
        self.oz = oz
        self.ncomp = ncomp
        self.tau = tau
        self.beta = beta
        
        # non standard atmosphere
        if ncomp == 1:
          taus = [tau]
          betas = [beta]    
          self.atm_nonstd = emul_atm.GetAllTransparencies(self.WL,am,pwv,oz,ncomp=ncomp, taus=taus, betas=betas, flagAerosols=True)
        else:
          self.atm_std = emul_atm.GetAllTransparencies(self.WL,am,pwv,oz)
          
        # non standard total filter (instrumental x atmosphere)  
        self.bandpass_total_nonstd = {} 
        for index,f in enumerate(filter_tagnames):
          self.bandpass_total_nonstd[f] = Bandpass(wavelen=self.WL,sb=self.bandpass_inst[f].sb * self.atm_nonstd)
        
        # Non standard Normalized response 
        self.phiArray_nonstd, _ = Sed().setup_phi_array([self.bandpass_total_nonstd[f] for f in filter_tagnames])
        
        # Integrals IIb0(non std) and IIb1(non std)
        self.all_II0_nonstd = {}
        self.all_II1_nonstd = {}
        self.all_II0ratio_nonstd = {}
        self.all_II1sub_nonstd = {}
        
        for index,f in enumerate(filter_tagnames):
    
          the_II0 = self.fII0(self.bandpass_total_nonstd[f].wavelen,self.bandpass_total_nonstd[f].sb)
          self.all_II0_nonstd[f] = the_II0
          self.all_II0ratio_nonstd[f] = the_II0/self.all_II0_std[f]
          the_II1 = self.fII1(self.WL,self.phiArray_nonstd[index,:],FILTERWL[index,2])
          self.all_II1_nonstd[f] = the_II1
          self.all_II1sub_nonstd[f] = self.all_II1_std[f] - the_II1
          
  def CalculateMultiObs(self,am,pwv,oz,ncomp,tau,beta):
        """
        """
        
        self.coll_atm_nonstd = []
        self.coll_bandpass_total_nonstd = []
        self.coll_phiArray_nonstd = []
        self.coll_all_II0_nonstd = []
        self.coll_all_II1_nonstd = []
        self.coll_all_II0ratio_nonstd = []
        self.coll_all_II1sub_nonstd = []
        
        if isinstance(am, list) or isinstance(am, np.ndarray):
          if isinstance(am, list):
            all_am = np.array(am)
          else:
            all_am = am
          for am in all_am:
            self.CalculateObs(am,pwv,oz,ncomp,tau,beta)
            self.coll_atm_nonstd.append(self.atm_nonstd)
            self.coll_bandpass_total_nonstd.append(self.coll_bandpass_total_nonstd)
            self.coll_phiArray_nonstd.append(self.phiArray_nonstd)
            self.coll_all_II0_nonstd.append(self.all_II0_nonstd) 
            self.coll_all_II1_nonstd.append(self.all_II1_nonstd)
            self.coll_all_II0ratio_nonstd.append(self.all_II0ratio_nonstd)
            self.coll_all_II1sub_nonstd.append(self.all_II1sub_nonstd)
                
                    
        elif isinstance(pwv, list) or isinstance(pwv, np.ndarray): 
          if isinstance(pwv, list):
            all_pwv = np.array(pwv)
          else:
            all_pwv = pwv
            
          for pwv in all_pwv:
            self.CalculateObs(am,pwv,oz,ncomp,tau,beta)
            self.coll_atm_nonstd.append(self.atm_nonstd)
            self.coll_bandpass_total_nonstd.append(self.coll_bandpass_total_nonstd)
            self.coll_phiArray_nonstd.append(self.phiArray_nonstd)
            self.coll_all_II0_nonstd.append(self.all_II0_nonstd) 
            self.coll_all_II1_nonstd.append(self.all_II1_nonstd)
            self.coll_all_II0ratio_nonstd.append(self.all_II0ratio_nonstd)
            self.coll_all_II1sub_nonstd.append(self.all_II1sub_nonstd)
                
         
        elif isinstance(oz, list) or isinstance(oz, np.ndarray): 
          if isinstance(oz, list):
            all_oz = np.array(oz)
          else:
            all_oz = oz
            
          for oz in all_oz:
            self.CalculateObs(am,pwv,oz,ncomp,tau,beta)
            self.coll_atm_nonstd.append(self.atm_nonstd)
            self.coll_bandpass_total_nonstd.append(self.coll_bandpass_total_nonstd)
            self.coll_phiArray_nonstd.append(self.phiArray_nonstd)
            self.coll_all_II0_nonstd.append(self.all_II0_nonstd) 
            self.coll_all_II1_nonstd.append(self.all_II1_nonstd)
            self.coll_all_II0ratio_nonstd.append(self.all_II0ratio_nonstd)
            self.coll_all_II1sub_nonstd.append(self.all_II1sub_nonstd)
        
        else:
          print("Not implemented yet")             
          
        
          
      
           
          
                  
            


        
################################################################################################        


def main():
    print("============================================================")
    print("Photometric Corrections                                     ")
    print("============================================================")
    
    # retrieve the path of data
    path_data =  os.path.join(atmosphtransmemullsst.__path__[0],'../data/simplegrid')
    # create emulator  
    emul = SimpleAtmEmulator(path = path_data)
    wl = [400.,800.,900.]
    am=1.2
    pwv =4.0
    oz=300.
    transm = emul.GetAllTransparencies(wl,am,pwv,oz)
    print("wavelengths (nm) \t = ",wl)
    print("transmissions    \t = ",transm)
    
    

if __name__ == "__main__":
    main()
