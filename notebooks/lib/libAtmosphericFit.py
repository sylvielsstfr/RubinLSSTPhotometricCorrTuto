# libAtmosphericFit
#
# Author          : Sylvie Dagoret-Campagne
# Affiliaton      : IJCLab/IN2P3/CNRS
# Creation Date   : 2023/02/18
#
# A python tool to fit quickly atmospheric transmission of Auxtel Spectra with scipy.optimize and
# using the atmospheric emulator atmosphtransmemullsst
# 
#
#
import os
import sys
from pathlib import Path

#import atmosphtransmemullsst
#from atmosphtransmemullsst.simpleatmospherictransparencyemulator import SimpleAtmEmulator

from scipy.optimize import curve_fit,least_squares
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import pickle

dir_path_data = "/../data/simplegrid" 


file_data_dict = {
    "info_training" :"atmospherictransparencygrid_params_training.pickle",
    "info_test" : "atmospherictransparencygrid_params_test.pickle",
    "data_rayleigh_training" : "atmospherictransparencygrid_rayleigh_training.npy",
    "data_rayleigh_test" : "atmospherictransparencygrid_rayleigh_test.npy",
    "data_o2abs_training" : "atmospherictransparencygrid_O2abs_training.npy",
    "data_o2abs_test" : "atmospherictransparencygrid_O2abs_test.npy",
    "data_pwvabs_training": "atmospherictransparencygrid_PWVabs_training.npy",
    "data_pwvabs_test":"atmospherictransparencygrid_PWVabs_test.npy",
    "data_ozabs_training" : "atmospherictransparencygrid_OZabs_training.npy",
    "data_ozabs_test" : "atmospherictransparencygrid_OZabs_test.npy"
}

def find_data_path():
    """
    Search the path for the atmospheric emulator
    """
    
    print("relativ data path ",dir_path_data)
    
    dir_file_abspath = os.path.dirname(os.path.abspath(__file__))
    print("abspath = ",dir_file_abspath)

    dir_file_realpath = os.path.dirname(os.path.realpath(__file__))
    print("realpath = ",dir_file_realpath)

    dir_file_sys = Path(sys.path[0])
    print("syspath = ",dir_file_sys)

    dir_file_dirname = os.path.dirname(__file__) 
    print("dirname = ",dir_file_dirname)
    
    path_data = dir_file_realpath + dir_path_data
    
    for key, filename in file_data_dict.items():
        file_path = os.path.join(path_data ,  filename)
        flag_found = os.path.isfile(file_path)  
        if flag_found :
            print(f"found data file {file_path}")
        else:
            print(f">>>>>>>>>> NOT found data file {file_path}")
            
    return path_data 


class SimpleAtmEmulator:
    """
    Emulate Atmospheric Transparency above LSST from a data grids
    extracted from libradtran and analytical functions for aerosols.
    There are 3 grids:
    - 2D grid Rayleigh transmission vs (wavelength,airmass)
    - 2D grid O2 absorption vs  (wavelength,airmass)
    - 3D grid for PWV absorption vs (wavelength,airmass,PWV)
    - 3D grid for Ozone absorption vs (wavelength,airmass,Ozone)
    - Aerosol transmission for any number of components
    
    """
    def __init__(self,path='../data/simplegrid'):
        """
        Initialize the class for data point files from which the 2D and 3D grids are created.
        Interpolation are calculated from the scipy RegularGridInterpolator() function
        """
        
        self.path = path
        self.fn_info_training = file_data_dict["info_training"]
        self.fn_info_test = file_data_dict["info_test"]
        self.fn_rayleigh_training = file_data_dict["data_rayleigh_training"]
        self.fn_rayleigh_test = file_data_dict["data_rayleigh_test"]
        self.fn_O2abs_training = file_data_dict["data_o2abs_training"]
        self.fn_O2abs_test = file_data_dict["data_o2abs_test"]
        self.fn_PWVabs_training = file_data_dict["data_pwvabs_training"]
        self.fn_PWVabs_test = file_data_dict[ "data_pwvabs_test"]
        self.fn_OZabs_training = file_data_dict["data_ozabs_training"]
        self.fn_OZabs_test = file_data_dict["data_ozabs_test"]
        

        self.info_params_training = None
        self.info_params_test = None
        self.data_rayleigh_training = None
        self.data_rayleigh_test = None
        self.data_O2abs_training = None
        self.data_O2abs_test = None
        self.data_PWVabs_training = None
        self.data_PWVabs_test = None
        self.data_OZabs_training = None
        self.data_OZabs_test = None
        
        self.loadtables()
        
        self.WLMIN = self.info_params_training["WLMIN"]
        self.WLMAX = self.info_params_training["WLMAX"]
        self.WLBIN = self.info_params_training["WLBIN"]
        self.NWLBIN = self.info_params_training['NWLBIN']
        self.WL = self.info_params_training['WL']
        
        self.AIRMASSMIN = self.info_params_training['AIRMASSMIN']
        self.AIRMASSMAX = self.info_params_training['AIRMASSMAX']
        self.NAIRMASS = self.info_params_training['NAIRMASS']
        self.DAIRMASS = self.info_params_training['DAIRMASS']
        self.AIRMASS = self.info_params_training['AIRMASS']
        
        self.PWVMIN = self.info_params_training['PWVMIN']
        self.PWVMAX = self.info_params_training['PWVMAX'] 
        self.NPWV = self.info_params_training['NPWV']
        self.DPWV = self.info_params_training['DPWV'] 
        self.PWV = self.info_params_training['PWV']
        
        
        self.OZMIN =  self.info_params_training['OZMIN']
        self.OZMAX = self.info_params_training['OZMAX']
        self.NOZ = self.info_params_training['NOZ']
        self.DOZ =  self.info_params_training['DOZ'] 
        self.OZ = self.info_params_training['OZ']
        
        
        self.lambda0 = 550.
        self.tau0 = 1.


        self.func_rayleigh_train = RegularGridInterpolator((self.WL,self.AIRMASS),self.data_rayleigh_training)
        self.func_O2abs_train = RegularGridInterpolator((self.WL,self.AIRMASS),self.data_O2abs_training)
        self.func_PWVabs_train = RegularGridInterpolator((self.WL,self.AIRMASS,self.PWV),self.data_PWVabs_training)
        self.func_OZabs_train = RegularGridInterpolator((self.WL,self.AIRMASS,self.OZ),self.data_OZabs_training)

        
        
    def loadtables(self):
        """
        Load files into grid arrays
        """
        
        filename=os.path.join(self.path,self.fn_info_training)     
        with open(filename, 'rb') as f:
            self.info_params_training = pickle.load(f)
            
        filename=os.path.join(self.path,self.fn_info_test)     
        with open(filename, 'rb') as f:
            self.info_params_test = pickle.load(f)        
        
        filename=os.path.join(self.path,self.fn_rayleigh_training)
        with open(filename, 'rb') as f:
            self.data_rayleigh_training=np.load(f)
            print("data_rayleigh_training",self.data_rayleigh_training.shape)
            
        filename=os.path.join(self.path,self.fn_rayleigh_test)
        with open(filename, 'rb') as f:
            self.data_rayleigh_test=np.load(f)
            print("data_rayleigh_test",self.data_rayleigh_test.shape)
            
        filename=os.path.join(self.path,self.fn_O2abs_training)
        with open(filename, 'rb') as f:
            self.data_O2abs_training=np.load(f)
            
        filename=os.path.join(self.path,self.fn_O2abs_test)
        with open(filename, 'rb') as f:
            self.data_O2abs_test=np.load(f)
                  
        filename=os.path.join(self.path,self.fn_PWVabs_training)
        with open(filename, 'rb') as f:
            self.data_PWVabs_training=np.load(f)
            
        filename=os.path.join(self.path,self.fn_PWVabs_test)
        with open(filename, 'rb') as f:
            self.data_PWVabs_test=np.load(f)
            
            
        filename=os.path.join(self.path,self.fn_OZabs_training)
        with open(filename, 'rb') as f:
            self.data_OZabs_training=np.load(f)
            
        filename=os.path.join(self.path,self.fn_OZabs_test)
        with open(filename, 'rb') as f:
            self.data_OZabs_test=np.load(f)
            
            
    def GetWL(self):
        return self.WL
            
    def GetRayleighTransparencyArray(self,wl,am):
        #pts = np.array([ np.array([the_wl,am],dtype=object)  for the_wl in wl ],dtype=object) 
        #return np.array([ self.func_rayleigh_train(pt) for pt in pts])
        pts = [ (the_wl,am) for the_wl in wl ]
        pts = np.array(pts)
        return self.func_rayleigh_train(pts)
    
    
    def GetO2absTransparencyArray(self,wl,am):
        #pts = np.array([ np.array([the_wl,am],dtype=object)   for the_wl in wl],dtype=object)
        #return np.array([ self.func_O2abs_train(pt) for pt in pts])
        pts = [ (the_wl,am) for the_wl in wl ]
        pts = np.array(pts)
        return self.func_O2abs_train(pts)
    
    
    def GetPWVabsTransparencyArray(self,wl,am,pwv):
        #pts = np.array([ np.array([the_wl,am,pwv],dtype=object)   for the_wl in wl],dtype=object)
        #return  np.array([self.func_PWVabs_train(pt) for pt in pts])
        pts = [ (the_wl,am,pwv) for the_wl in wl ]
        pts = np.array(pts)
        return self.func_PWVabs_train(pts)
    
    
    def GetOZabsTransparencyArray(self,wl,am,oz):
        #pts = np.array([ np.array([the_wl,am,oz],dtype=object)   for the_wl in wl],dtype=object)
        #return np.array([self.func_OZabs_train(pt) for pt in pts])
        pts = [ (the_wl,am,oz) for the_wl in wl ]
        pts = np.array(pts)
        return self.func_OZabs_train(pts)
    
    
    def GetGriddedTransparencies(self,wl,am,pwv,oz,flagRayleigh=True,flagO2abs=True,flagPWVabs=True,flagOZabs=True):
        """
        Emulation of libradtran simulated transparencies. Decomposition of the
        total transmission in different processes:
        - Rayleigh scattering
        - O2 absorption
        - PWV absorption
        - Ozone absorption
        
        inputs:
        - wl : wavelength array or list
        - am :the airmass,
        - pwv : the precipitable water vapor (mm)
        - oz : the ozone column depth in Dobson unit
        - flags to activate or not the individual interaction processes
        
        outputs:
        - 1D array of atmospheric transmission (save size as wl)
        
        """
        

        if flagRayleigh:
            transm = self.GetRayleighTransparencyArray(wl,am)
        else:
            transm = np.ones(len(wl))
            
        if flagO2abs:
            transm *= self.GetO2absTransparencyArray(wl,am)
            
        if flagPWVabs:
            transm *= self.GetPWVabsTransparencyArray(wl,am,pwv)
            
        if flagOZabs:
            transm *= self.GetOZabsTransparencyArray(wl,am,oz)
            
        # if pts are array of array, RegularGridInterpolator returns an array of (n,1)
        # In that case, one must flatten the array
        if transm.ndim >= 2:
            transm =transm.squeeze()
        return transm
            
    def GetAerosolsTransparencies(self,wl,am,ncomp,taus=None,betas=None):
        """
        Compute transmission due to aerosols:
        
        inputs:
        - wl : wavelength array
        - am : the airmass
        - ncomp : the number of aerosol components
        - taus : the vertical aerosol depth of each component at lambda0 vavelength
        - betas : the angstrom exponent. Must be negativ.
        
        
        outputs:
        - 1D array of atmospheric transmission (save size as wl)
        
        """
        if not isinstance(wl, np.ndarray):   
            wl = np.array(wl,dtype=object)
        NWL=wl.shape[0]
        
        transm = np.ones(NWL)
        
        if ncomp <=0:
            return transm
        else:
            taus=np.array(taus)
            betas=np.array(betas)
            
            NTAUS=taus.shape[0]
            NBETAS=betas.shape[0]
        
            assert ncomp<=NTAUS
            assert ncomp<=NBETAS     
        
            for icomp in range(ncomp):            
                exponent = (taus[icomp]/self.tau0)*np.exp(betas[icomp]*np.log(wl/self.lambda0))*am
                transm *= np.exp(-exponent)
            
            return transm
        
        
    def GetAllTransparencies(self,wl,am,pwv,oz,ncomp=0, taus=None, betas=None, flagRayleigh=True,flagO2abs=True,flagPWVabs=True,flagOZabs=True,flagAerosols=False):
        """
        Combine interpolated libradtran transmission with analytical expression for the
        aerosols
        
        inputs:
        - wl : wavelength array or list
        - am :the airmass,
        - pwv : the precipitable water vapor (mm)
        - oz : the ozone column depth in Dobson unit
        - ncomp : number of aerosols components,
        - taus & betas : arrays of parameters for aerosols
        - flags to activate or not the individual interaction processes
        
        outputs:
        - 1D array of atmospheric transmission (save size as wl)
        
        """
        
        
        transm = self.GetGriddedTransparencies(wl,am,pwv,oz,flagRayleigh=flagRayleigh,flagO2abs=flagO2abs,flagPWVabs=flagPWVabs,flagOZabs=flagOZabs)
        
        if flagAerosols:
            transmaer = self.GetAerosolsTransparencies(wl,am,ncomp,taus,betas)
            transm *=transmaer
           
            
        return transm
            
print("Path at terminal when executing this file")
print(os.getcwd() + "\n")
current_path = os.getcwd()

print("This file path, relative to os.getcwd()")
print(__file__ + "\n")

print("This file full path (following symlinks)")
full_path = os.path.realpath(__file__)
print(full_path + "\n")

print("This file directory and name")
path, filename = os.path.split(full_path)
print(path + ' --> ' + filename + "\n")

print("This file directory only")
print(os.path.dirname(full_path))

print("The data path is")
data_path = current_path + "/../data/simplegrid"   
print(data_path)
        
# global variable
#emul = SimpleAtmEmulator(os.path.join(atmosphtransmemullsst.__path__[0],'../data/simplegrid'))
#emul = SimpleAtmEmulator(path='/Users/sylvie/MacOSX/GitHub/LSST/AuxTelComm/notebooks_usdf/FitAtmosphericParameters/data/simplegrid')
emul = SimpleAtmEmulator(path=data_path)

#emul = SimpleAtmEmulator(path='/Users/dagoret/MacOSX/GitHub/LSST/AuxTelComm/notebooks_usdf/FitAtmosphericParameters/data/simplegrid')
 
 
    
# =========    Functions for the fit grey factor, pwv, ozone ========================
def fluxpred_greypwvo3(params,*arg):
    """
    Model prediction y = f(x) where x is the array of wavelength and y the measured flux
    
    inputs:
      - params : the unknown atmospheric parameters to fit (pwv,oz)
      - *args  : additionnal data provided to compute the model : 
                   1) the wavelength array, 
                   3) sedxthroughput,
                   2) the airmass
      
    output:
      - the array of predicted flux values at each wavelength data point
      
    It computes the production of the atmospheric transparency by the SED - throughout product  
    It makes use of the atmospheric model transparency emulated in the external global object emul 
    which depend on atmospheric parameters to fit
    
    """
    
    global emul
    
    alpha,pwv,oz = params    
    
    # have to build the wavelength array as follow (I don't know why)
    #wl = []
    #for el in arg[0]:
    #    wl.append(el)
    
   
    wl = arg[0]
    the_sedxthroughput = arg[2]
    airmass = arg[1]
        

    fl = alpha*emul.GetAllTransparencies(wl ,airmass,pwv,oz,ncomp=0,flagAerosols=False)
    
    fl *= the_sedxthroughput
    return fl


def func_residuals_greypwvo3(params,*arg):
    """
    This function is called by the scipy.optimize.least_squares function to compute the normalized residuals
    at each data point at each wavelength (yi-f(xi,thetaj)/sigmai
    
    input:
      - params : the unknown atmospheric parameters to fit (pwv,oz)
      - *args  : additionnal data provided to compute the mode : 
                 to the model
                   1) the wavelength array, 
                   2) flux data : the measured fluxes
                   3) dataerr : the error on fluxes
                   4) the airmass
                   5) sedxthroughput,
                  
    
    return:
      - the array of residuals for the given set of parameters required by the least_squares function 
    """
    
    
    
    the_wl = arg[0]         # wavelength array required by the model function pred2
    the_sedthr = arg[4]     # the product of the SED by the throughput required by the model function pred2
    the_airmass = arg[3]    # the airmass required by the model function pred2 
    the_data = arg[1]       # the observed fluxes array
    the_dataerr = arg[2]    # the error on the observed fluxes array
    
     
    alpha,pwv , oz = params   # decode the parameters
    
    
    flux_model = fluxpred_greypwvo3(params,the_wl,the_airmass,the_sedthr)
    
    
    residuals = (the_data - flux_model)/the_dataerr
    
 
    
    return residuals

# =========    Functions for the fit grey factor, pwv, No ozone ========================
def fluxpred_greypwv(params,*arg):
    """
    Model prediction y = f(x) where x is the array of wavelength and y the measured flux
    
    inputs:
      - params : the unknown atmospheric parameters to fit (pwv,oz)
      - *args  : additionnal data provided to compute the model : 
                   1) the wavelength array, 
                   3) sedxthroughput,
                   2) the airmass
      
    output:
      - the array of predicted flux values at each wavelength data point
      
    It computes the production of the atmospheric transparency by the SED - throughout product  
    It makes use of the atmospheric model transparency emulated in the external global object emul 
    which depend on atmospheric parameters to fit
    
    """
    
    global emul
    
    alpha,pwv = params   
    oz=300. 
    
    # have to build the wavelength array as follow (I don't know why)
    #wl = []
    #for el in arg[0]:
    #    wl.append(el)
    
   
    wl = arg[0]
    the_sedxthroughput = arg[2]
    airmass = arg[1]
        

    fl = alpha*emul.GetAllTransparencies(wl ,airmass,pwv,oz,ncomp=0,flagAerosols=False)
    
    fl *= the_sedxthroughput
    return fl


def func_residuals_greypwv(params,*arg):
    """
    This function is called by the scipy.optimize.least_squares function to compute the normalized residuals
    at each data point at each wavelength (yi-f(xi,thetaj)/sigmai
    
    input:
      - params : the unknown atmospheric parameters to fit (pwv,oz)
      - *args  : additionnal data provided to compute the mode : 
                 to the model
                   1) the wavelength array, 
                   2) flux data : the measured fluxes
                   3) dataerr : the error on fluxes
                   4) the airmass
                   5) sedxthroughput,
                  
    
    return:
      - the array of residuals for the given set of parameters required by the least_squares function 
    """
    
    
    
    the_wl = arg[0]         # wavelength array required by the model function pred2
    the_sedthr = arg[4]     # the product of the SED by the throughput required by the model function pred2
    the_airmass = arg[3]    # the airmass required by the model function pred2 
    the_data = arg[1]       # the observed fluxes array
    the_dataerr = arg[2]    # the error on the observed fluxes array
    
     
    alpha,pwv = params      # decode the parameters
    
    
    flux_model = fluxpred_greypwv(params,the_wl,the_airmass,the_sedthr)
    
    
    residuals = (the_data - flux_model)/the_dataerr
    
    return residuals



# 
class FitAtmosphericParams:
    """
    Class to handle a buch of different methods for fitting
    """
    def __init__(self):
        print("Init FitAtmosphericParams")


    def fit_greypwvo3(self,params0,xdata,ydata,yerrdata,airmass,sedxthroughput):
        """
        Fit a grey term, precipitable water varpor and ozone

        """
            
        res_fit = least_squares(func_residuals_greypwvo3, params0,bounds=([0.1,0.001,50.],[2,9.5,550.]),args=[xdata,ydata,yerrdata,airmass,sedxthroughput])
        
        alpha_fit,pwv_fit,oz_fit = res_fit.x
    
    
        # results fo the fit
        popt= res_fit.x
        J = res_fit.jac
        pcov = np.linalg.inv(J.T.dot(J))
        sigmas = np.sqrt(np.diagonal(pcov))
        cost=res_fit.cost
        residuals = res_fit.fun
        chi2=np.sum(residuals**2)
        ndeg = J.shape[0]- popt.shape[0]
        chi2_per_deg = chi2/ndeg
        pwve = sigmas[1]*np.sqrt(chi2_per_deg)
        oze =sigmas[2]*np.sqrt(chi2_per_deg)
        greye = sigmas[0]*np.sqrt(chi2_per_deg)
        
        fit_dict = {"chi2":chi2,"ndeg":ndeg,"chi2_per_deg":chi2_per_deg,"popt":popt,"sigmas":sigmas,"pwv_fit":pwv_fit,"oz_fit":oz_fit,"grey_fit":alpha_fit,"pwve":pwve,"oze":oze,"greye":greye}

        return res_fit,fit_dict
    
    
    def pred_greypwvo3(self,params,xdata,airmass,sedxthroughput):
        """
        """
        flux_model = fluxpred_greypwvo3(params,xdata,airmass,sedxthroughput)
        return flux_model
    
#--------
    
    def fit_greypwv(self,params0,xdata,ydata,yerrdata,airmass,sedxthroughput):
        """
        Fit a grey term, precipitable water varpor and ozone

        """
            
        res_fit = least_squares(func_residuals_greypwv, params0,bounds=([0.1,0.001],[2,9.5]),args=[xdata,ydata,yerrdata,airmass,sedxthroughput])
        
        alpha_fit,pwv_fit = res_fit.x
    
    
        # results fo the fit
        popt= res_fit.x
        J = res_fit.jac
        pcov = np.linalg.inv(J.T.dot(J))
        sigmas = np.sqrt(np.diagonal(pcov))
        cost=res_fit.cost
        residuals = res_fit.fun
        chi2=np.sum(residuals**2)
        ndeg = J.shape[0]- popt.shape[0]
        chi2_per_deg = chi2/ndeg
        pwve = sigmas[1]*np.sqrt(chi2_per_deg)
        greye = sigmas[0]*np.sqrt(chi2_per_deg)
        
        fit_dict = {"chi2":chi2,"ndeg":ndeg,"chi2_per_deg":chi2_per_deg,"popt":popt,"sigmas":sigmas,"pwv_fit":pwv_fit,"grey_fit":alpha_fit,"pwve":pwve,"greye":greye}

        return res_fit,fit_dict
    
    
    def pred_greypwv(self,params,xdata,airmass,sedxthroughput):
        """
        """
        flux_model = fluxpred_greypwv(params,xdata,airmass,sedxthroughput)
        return flux_model
        
        


def main():
    print("============================================================")
    print("Simple Atmospheric emulator for Rubin-LSST observatory")
    print("============================================================")
    
    # retrieve the path of data
    path_data =  find_data_path()  
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