# libPhotometricCorrections_auxtel.py
#
# Author          : Sylvie Dagoret-Campagne
# Affiliaton      : IJCLab/IN2P3/CNRS
# Creation Date   : 2024/01/03
# Last update     : 2024/04/06
#
# A python tool to calculate Photometric Correction
# 
# We use Auxtel transmission with SDSS filters
#
#
import os
import sys
from pathlib import Path

from getObsAtmo import ObsAtmo

from scipy import interpolate
import numpy as np

from rubinsimphot.phot_utils import Bandpass, Sed
from rubinsimphot.data.data_sets import  get_data_dir
from rubinsimphot.phot_utils import PhotometricParameters 

#README.md        darksky.dat      filter_r.dat     hardware_g.dat   hardware_y.dat   lens3.dat        total_g.dat      total_y.dat
#README_SOURCE.md detector.dat     filter_u.dat     hardware_i.dat   hardware_z.dat   m1.dat           total_i.dat      total_z.dat
#atmos_10.dat     filter_g.dat     filter_y.dat     hardware_r.dat   lens1.dat        m2.dat           total_r.dat      version_info
#atmos_std.dat    filter_i.dat     filter_z.dat     hardware_u.dat   lens2.dat        m3.dat           total_u.dat

#/Users/dagoret/MacOSX/GitHub/LSST/AtmosphericSimulation/rubinsimphot/src/rubin_sim_data/throughputs/auxtel>ls
#auxtel_sdss_g.dat                                                                auxtel_sdss_u.dat
#auxtel_sdss_i.dat                                                                auxtel_sdss_z.dat
#auxtel_sdss_r.dat                                                                multispectra_holo4_003_HD142331_20230802_AuxTel_doGainsPTC_v3.0.3_throughput.txt


filter_tagnames = ["u","g","r","i","z","y"]
Filter_tagnames = ["U","G","R","I","Z","Y"]
filtercolor_tagnames = ["u-g","g-r","r-i","i-z","z-y"]
Filtercolor_tagnames = ["U-G","G-R","R-I","I-Z","Z-Y"]
filter_color = ["b","g","r","orange","grey","k"]

# New version with prime sdss filters (March-April 2024)
hardware_filenames = ["auxtel_sdss_up_total.dat","auxtel_sdss_gp_total.dat","auxtel_sdss_rp_total.dat","auxtel_sdss_ip_total.dat","auxtel_sdss_zp_total.dat","auxtel_sdss_yp_total.dat"] 
filter_filenames = ["auxtel_sdss_up.dat","auxtel_sdss_gp.dat","auxtel_sdss_rp.dat","auxtel_sdss_ip.dat","auxtel_sdss_zp.dat" ,"auxtel_sdss_yp.dat"]
total_filenames = ["auxtel_sdss_up_total.dat","auxtel_sdss_gp_total.dat","auxtel_sdss_rp_total.dat","auxtel_sdss_ip_total.dat","auxtel_sdss_zp_total.dat","auxtel_sdss_yp_total.dat"]
filter_tagnames = ["u","g","r","i","z","y"]
Filter_tagnames = ["U","G","R","I","Z","Y"]
filtercolor_tagnames = ["u-g","g-r","r-i","i-z","z-y"]
Filtercolor_tagnames = ["U-G","G-R","R-I","I-Z","Z-Y"]
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



FILTERWL = np.array([[ 353.        ,  385.        ,  369.        ,   32.        ],
                  [ 393.        ,  560.        ,  476.5       ,  167.        ],
                  [ 557.        ,  703.        ,  630.        ,  146.        ],
                  [ 688.        ,  859.        ,  773.5       ,  171.        ],
                  [ 812.        ,  938.        ,  875.76271186,  126.        ],
                  [ 934.        , 1060.        ,  997.        ,  126.        ]])


F0 = 3631.0 # Jy 1, Jy = 10^{-23} erg.cm^{-2}.s^{-1}.Hz^{-1}
Jy_to_ergcmm2sm1hzm1 = 1e-23
DT = 30.0 # seconds
gel = 1.08269375
#hP = 6.62607015E-34 # J⋅Hz−1
hP = 6.626196E-27
A  = 9636.0 # cm2
pixel_scale = 0.1 #arcsec/pixel
readnoise = 8.96875

#ZPT_cont =  2.5 \log_{10} \left(\frac{F_0 A \Delta T}{g_{el} h} \right)
ZPTconst = 2.5*np.log10(F0*Jy_to_ergcmm2sm1hzm1*A*DT/gel/hP)

def set_photometric_parameters(exptime, nexp, readnoise=readnoise):
    # readnoise = None will use the default (8.8 e/pixel). Readnoise should be in electrons/pixel.
    photParams = PhotometricParameters(exptime=exptime, nexp=nexp, readnoise=readnoise)
    return photParams

photoparams = set_photometric_parameters(DT, 1 , readnoise=readnoise )
photoparams._gain = gel
photoparams._exptime = DT
photoparams._effarea = A
photoparams._platescale = pixel_scale


def fII0(wl,s):
  return np.trapz(s/wl,wl)
      
def fII1(wl,phi,wlb):
  return np.trapz(phi*(wl-wlb),wl)
  
def ZPT(wl,s):
  return 2.5*np.log10(fII0(wl,s)) + ZPTconst

#print("libPhotometricCorrections.py :: Use atmosphtransmemullsst.__path__[0],'../data/simplegrid as the path to data")
#data_path = os.path.join(atmosphtransmemullsst.__path__[0],'../data/simplegrid')
#print(f"libPhotometricCorrections :: data_path = {data_path}")


       
# The emulator as a global variable
#emul_atm = SimpleAtmEmulator(data_path)
emul_atm = ObsAtmo()


class PhotometricCorrections:
  def __init__(self,am0=1.2,pwv0=5.0,oz0=300.,tau0=0.0,beta0=1.2):
        """
        Constructor
        """
        global emul_atm,sed
        self.WL = emul_atm.GetWL()
        
        # standard atmosphere parameters
        self.am0 = am0
        self.pwv0 = pwv0
        self.oz0 = oz0
        self.tau0 = tau0
        self.beta0 = beta0
        
        # standard atmosphere
        self.atm_std = emul_atm.GetAllTransparencies(self.WL,am0,pwv0,oz0, tau0, beta0)
        
          
        # instrumental filter
        self.bandpass_inst = {} 
        #path_rubin_sim_throughput=os.path.join(os.getenv("HOME"),"rubin_sim_data/throughputs/auxtel")
        fdir = get_data_dir()
        path_rubin_sim_throughput = os.path.join(fdir, 'throughputs', 'auxtel')

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
        self.all_ZP = {}

        # loop on filters
        for index,f in enumerate(filter_tagnames):
          the_II0 = self.fII0(self.bandpass_total_std[f].wavelen,self.bandpass_total_std[f].sb)
          self.all_II0_std[f] = the_II0
          the_II1 = self.fII1(self.WL,self.phiArray_std[index,:],FILTERWL[index,2])
          self.all_II1_std[f] = the_II1
          self.all_ZP[f] = 2.5*np.log10(the_II0) + ZPTconst
          
          
        # Non standard calculations will be calculated later after initialisation
        self.am = 1.2
        self.pwv = 0
        self.oz = 0
        self.tau = 0.04
        self.beta = 1
        
        self.atm_nonstd = None
        self.bandpass_total_nonstd = None
        self.phiArray_nonstd = None
        self.all_II0_nonstd = None
        self.all_II1_nonstd = None
        self.all_II0ratio_nonstd = None
        self.all_II1sub_nonstd = None
        self.all_ZP_nonstd  = None
        
        self.coll_atm_nonstd = None
        self.coll_bandpass_total_nonstd = None
        self.coll_phiArray_nonstd = None
        self.coll_all_II0_nonstd = None
        self.coll_all_II1_nonstd = None
        self.coll_all_II0ratio_nonstd = None
        self.coll_all_II1sub_nonstd = None
        self.coll_allZP_nonstd = None
        
        self.allparameters = None
        self.allcollperfilter = None
        
        
         
          
  def fII0(self,wl,s):
        return np.trapz(s/wl,wl)
      
  def fII1(self,wl,phi,wlb):
        return np.trapz(phi*(wl-wlb),wl)
      
  def CalculatePerfilter(self):
        """
        Given a set of non standard parameter values, the integrals are calculated for each
        Filter.
        """
        
        self.allcollperfilter = {} # init the main dictionary

        # loop on filters
        for f in filter_tagnames:
          list_II0_nonstd = []
          list_II1_nonstd = []
          list_II0ratio_nonstd = []
          list_II1sub_nonstd = []
          list_ZPT_nonstd = []
          
          # loop on the different parameter set conditions to build a list    
          for idx,param in enumerate(self.allparameters):
            list_II0_nonstd.append(self.coll_all_II0_nonstd[idx][f])
            list_II1_nonstd.append(self.coll_all_II1_nonstd[idx][f])
            list_II0ratio_nonstd.append(self.coll_all_II0ratio_nonstd[idx][f])
            list_II1sub_nonstd.append(self.coll_all_II1sub_nonstd[idx][f])
            list_ZPT_nonstd.append(self.coll_allZP_nonstd[idx][f])

           

          # create a dictionnary of lists , the keys are the filter values  
          # fill the dictionnary of integrals for that filter     
          filter_dict = {}
          filter_dict["II0_nonstd"] = np.array(list_II0_nonstd)
          filter_dict["II1_nonstd"] = np.array(list_II1_nonstd)
          filter_dict["II0ratio_nonstd"] = np.array(list_II0ratio_nonstd)
          filter_dict["II1sub_nonstd"] = np.array(list_II1sub_nonstd)
          filter_dict["ZPT_nonstd"] = np.array(list_ZPT_nonstd)
          
          self.allcollperfilter[f] = filter_dict     
                    
              
              
      
  def CalculateObs(self,am=1.2,pwv=5.0,oz=300,tau=0.0,beta=1.2):
        """
        Calculate the integrals for that atmospheric parameter condition
        """
        self.am = am
        self.pwv = pwv 
        self.oz = oz
        self.tau = tau
        self.beta = beta
        
        # non standard atmosphere
        #if tau >0.0 :
        self.atm_nonstd = emul_atm.GetAllTransparencies(self.WL,am,pwv,oz, tau,beta)
        #else:
        #  self.atm_std = emul_atm.GetAllTransparencies(self.WL,am,pwv,oz)
          
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
        self.all_ZP_nonstd = {}
        
        # loop on filters
        for index,f in enumerate(filter_tagnames):    
          the_II0 = self.fII0(self.bandpass_total_nonstd[f].wavelen,self.bandpass_total_nonstd[f].sb)
          self.all_II0_nonstd[f] = the_II0
          self.all_II0ratio_nonstd[f] = the_II0/self.all_II0_std[f]
          the_II1 = self.fII1(self.WL,self.phiArray_nonstd[index,:],FILTERWL[index,2])
          self.all_II1_nonstd[f] = the_II1
          self.all_II1sub_nonstd[f] = self.all_II1_std[f] - the_II1
          self.all_ZP_nonstd[f] =  2.5*np.log10(the_II0) + ZPTconst
          
  def CalculateMultiObs(self,am,pwv,oz,tau,beta):
        """
        At least on parameters is an array for varying conditions 
        - loop on CalculateObs
        """
        
        # reset the lists of integrals per varying condtion parameter
        self.coll_atm_nonstd = []
        self.coll_bandpass_total_nonstd = []
        self.coll_phiArray_nonstd = []
        self.coll_all_II0_nonstd = []
        self.coll_all_II1_nonstd = []
        self.coll_all_II0ratio_nonstd = []
        self.coll_all_II1sub_nonstd = []
        self.coll_allZP_nonstd = []
        
        # test if airmass is the varying parameter
        if isinstance(am, list) or isinstance(am, np.ndarray):
          if isinstance(am, list):
            all_am = np.array(am)
          else:
            all_am = am
            
          self.allparameters = all_am
          
          # loop on airmass
          for am in all_am:

            #calculate all the integrals 
            self.CalculateObs(am,pwv,oz,tau,beta)
            # copy the dictionaries (filters key) of calculated quantities
            self.coll_atm_nonstd.append(self.atm_nonstd)
            self.coll_bandpass_total_nonstd.append(self.bandpass_total_nonstd)
            self.coll_phiArray_nonstd.append(self.phiArray_nonstd)
            self.coll_all_II0_nonstd.append(self.all_II0_nonstd) 
            self.coll_all_II1_nonstd.append(self.all_II1_nonstd)
            self.coll_all_II0ratio_nonstd.append(self.all_II0ratio_nonstd)
            self.coll_all_II1sub_nonstd.append(self.all_II1sub_nonstd)
            self.coll_allZP_nonstd.append(self.all_ZP_nonstd) 
                
        # test if pwv is the varying parameter            
        elif isinstance(pwv, list) or isinstance(pwv, np.ndarray): 
          if isinstance(pwv, list):
            all_pwv = np.array(pwv)
          else:
            all_pwv = pwv
            
          self.allparameters = all_pwv
            
          # loop on pwv
          for pwv in all_pwv:
            self.CalculateObs(am,pwv,oz,tau,beta)
            self.coll_atm_nonstd.append(self.atm_nonstd)
            self.coll_bandpass_total_nonstd.append(self.bandpass_total_nonstd)
            self.coll_phiArray_nonstd.append(self.phiArray_nonstd)
            self.coll_all_II0_nonstd.append(self.all_II0_nonstd) 
            self.coll_all_II1_nonstd.append(self.all_II1_nonstd)
            self.coll_all_II0ratio_nonstd.append(self.all_II0ratio_nonstd)
            self.coll_all_II1sub_nonstd.append(self.all_II1sub_nonstd)
            self.coll_allZP_nonstd.append(self.all_ZP_nonstd) 
                
        # test if ozone is the varying parameter   
        elif isinstance(oz, list) or isinstance(oz, np.ndarray): 
          if isinstance(oz, list):
            all_oz = np.array(oz)
          else:
            all_oz = oz
            
          self.allparameters = all_oz
            
          # loop on ozone
          for oz in all_oz:
            self.CalculateObs(am,pwv,oz,tau,beta)
            self.coll_atm_nonstd.append(self.atm_nonstd)
            self.coll_bandpass_total_nonstd.append(self.bandpass_total_nonstd)
            self.coll_phiArray_nonstd.append(self.phiArray_nonstd)
            self.coll_all_II0_nonstd.append(self.all_II0_nonstd) 
            self.coll_all_II1_nonstd.append(self.all_II1_nonstd)
            self.coll_all_II0ratio_nonstd.append(self.all_II0ratio_nonstd)
            self.coll_all_II1sub_nonstd.append(self.all_II1sub_nonstd)
            self.coll_allZP_nonstd.append(self.all_ZP_nonstd) 

        # test if aerosols VAOD is the varying parameter    
        elif isinstance(tau, list) or isinstance(tau, np.ndarray): 
          if isinstance(tau, list):
            all_tau = np.array(tau)
          else:
            all_tau = tau
          
          
          self.allparameters = all_tau
          
          # loop on tau
          for tau in all_tau:
            self.CalculateObs(am,pwv,oz,tau,beta)
            self.coll_atm_nonstd.append(self.atm_nonstd)
            self.coll_bandpass_total_nonstd.append(self.bandpass_total_nonstd)
            self.coll_phiArray_nonstd.append(self.phiArray_nonstd)
            self.coll_all_II0_nonstd.append(self.all_II0_nonstd) 
            self.coll_all_II1_nonstd.append(self.all_II1_nonstd)
            self.coll_all_II0ratio_nonstd.append(self.all_II0ratio_nonstd)
            self.coll_all_II1sub_nonstd.append(self.all_II1sub_nonstd)
            self.coll_allZP_nonstd.append(self.all_ZP_nonstd)  
        
        # test if aerosols beta is the varying parameter    
        elif isinstance(beta, list) or isinstance(beta, np.ndarray): 
          if isinstance(beta, list):
            all_beta = np.array(beta)
          else:
            all_beta = beta
          
          
          self.allparameters = all_beta
          
          # loop on beta
          for beta in all_beta:
            self.CalculateObs(am,pwv,oz,tau,beta)
            self.coll_atm_nonstd.append(self.atm_nonstd)
            self.coll_bandpass_total_nonstd.append(self.bandpass_total_nonstd)
            self.coll_phiArray_nonstd.append(self.phiArray_nonstd)
            self.coll_all_II0_nonstd.append(self.all_II0_nonstd) 
            self.coll_all_II1_nonstd.append(self.all_II1_nonstd)
            self.coll_all_II0ratio_nonstd.append(self.all_II0ratio_nonstd)
            self.coll_all_II1sub_nonstd.append(self.all_II1sub_nonstd)
            self.coll_allZP_nonstd.append(self.all_ZP_nonstd) 

        else:
          print("Not implemented yet")  

        # copy quantities into list per atm condition  
        self.CalculatePerfilter() 
          
                  
          
        
          
      
           
          
                  
            


        
################################################################################################        


def main():
    print("============================================================")
    print("Photometric Corrections for Auxtel                          ")
    print("============================================================")
    
  
    # create emulator  
    # from getObsAtmo.getObsAtmo import ObsAtmo

    emul =  ObsAtmo("AUXTEL")
   
    wl = [400.,800.,900.]
    am=1.2
    pwv =4.0
    oz=300.
    transm = emul.GetAllTransparencies(wl,am,pwv,oz)
    print("wavelengths (nm) \t = ",wl)
    print("transmissions    \t = ",transm)
    
    

if __name__ == "__main__":
    main()
