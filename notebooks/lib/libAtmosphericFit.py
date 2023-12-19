# libAtmosphericFit
#
# Author          : Sylvie Dagoret-Campagne
# Affiliaton      : IJCLab/IN2P3/CNRS
# Creation Date   : 2023/02/18
# Last update     : 2023/12/19
#
# A python tool to fit quickly atmospheric transmission of Auxtel Spectra with scipy.optimize and
# using the atmospheric emulator atmosphtransmemullsst
# 
#
#
import os
import sys
from pathlib import Path

#from ObsAtmo import ObsAtmo
from getObsAtmo import ObsAtmo




#import atmosphtransmemullsst
#from atmosphtransmemullsst.simpleatmospherictransparencyemulator import SimpleAtmEmulator

from scipy.optimize import curve_fit,least_squares
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import pickle




#print("libAtmosphericFit.py :: Use atmosphtransmemullsst.__path__[0],'../data/simplegrid as the path to data")
#data_path = os.path.join(atmosphtransmemullsst.__path__[0],'../data/simplegrid')
#print(f"libAtmosphericFit.py :: data_path = {data_path}")

       
# The emulator as a global variable
#emul = SimpleAtmEmulator(data_path)
emul = ObsAtmo()

#######################     Fuctions for fitting   ###########################################
    
# =================    Functions for the fit grey factor, pwv, ozone ========================
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

# ================    Functions for the fit grey factor, pwv, No ozone ========================
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

# ================    Functions for the atmospheric profile airmass, pwv, ozone, aerosols ========
def transmpred_ampwvozaer(params,*arg):
    """
    Model prediction y = f(x) where x is the array of wavelength and y the measured flux
    
    inputs:
      - params : the unknown atmospheric parameters to fit (am,pwv,oz,vaod,beta)
      - *args  : additionnal data provided to compute the model : 
                   1) the wavelength array, 
                  
      
    output:
      - the array of predicted air transmission  at each wavelength data point
      
    """
    
    global emul
    
    # decode parameters
    am,pwv,oz,vaod,beta = params   
    
    # decode the argument
    wl = arg[0]
    
    transm = emul.GetAllTransparencies(wl ,am,pwv,oz,ncomp=1,taus=[vaod],betas=[beta],flagAerosols=True)
     
    return transm


def func_residuals_tramampwvozaer(params,*arg):
    """
    This function is called by the scipy.optimize.least_squares function to compute the normalized residuals
    at each data point at each wavelength (yi-f(xi,thetaj)/sigmai
    
    input:
      - params : the unknown atmospheric parameters to fit (pwv,oz)
      - *args  : additionnal data provided to compute the mode : 
                 to the model
                   1) the wavelength array, 

                  
    
    return:
      - the array of residuals for the given set of parameters required by the least_squares function 
    """
    
    
    
    the_wl = arg[0]         # wavelength array required by the model function pred2
    the_data = arg[1]       # the observed fluxes array
    
    
     
    am,pwv,oz,vaod,beta = params      # decode the parameters
    
    
    transm_model = transmpred_ampwvozaer(params,the_wl)
    
    
    residuals = (the_data - transm_model)
    
    return residuals

########################## Fitting class ############################################ 
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
    
    
    #-------------------

    def fit_transmampwvozaer(self,params0,xdata,ydata):
        """
        Fit transmission

        """            
        res_fit = least_squares(func_residuals_tramampwvozaer, params0,bounds=([1.001,0.001,0.001,0,-3],[2.499,9.5,599,0.5,0]),args=[xdata,ydata])
        
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
        
        
        fit_dict = {"chi2":chi2,"ndeg":ndeg,"chi2_per_deg":chi2_per_deg,"popt":popt,"sigmas":sigmas}

        return res_fit,fit_dict
    
    
    def pred_transmampwvozaer(self,params,xdata):
        """
        """
        transm_model = transmpred_ampwvozaer(params,xdata)
        return transm_model
        
################################################################################################        


def main():
    print("============================================================")
    print("Simple Atmospheric emulator for Rubin-LSST observatory")
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
