# RubinLSSTPhotometricCorrTuto

**Photometric Correction tutorail for virtual DESC collaboration meeting 27-Feb - 3 March 2023**

- authors : Sylvie Dagoret-Campagne & Martin Rodriguez Monroy
- affiliation : IJCLab/IN2P3/CNRS
- Creaiton date : February 2023
- Last verification : March 3rd 2023

## Goal

Understand the impact of atmospheric transparency variations on photometry, the necessity to estimate those variations and the necessity to correct them up to some level.


## Program

### Part 01 : Estimation of atmospheric parameters

After a general introduction, an example of fitting atmospheric parameters from spectra are show.

First steps in this tutorial, may be skipped.

### Part 02 : Photometric Corrections
The core of this tutorial are how to calculate Photometric correction.

It is possible to exploit SED-color dependence on observed fluxes in order to guess the SED shape (average slope) in each band.


### Part 03 : Applications
Study the impact on science after applied photometric corrections.

Still under developpement.



## Code Installation

### a) Installation of this tutorial `RubinLSSTPhotometricCorrTuto`

    git clone https://github.com/sylvielsstfr/RubinLSSTPhotometricCorrTuto
      

### b) Installation of the package `atmosphtransmemullsst`


Package to emulate atmospheric transparency


    git clone https://github.com/sylvielsstfr/git@github.com:LSSTDESC/atmosphtransmemullsst.git
    cd atmosphtransmemullsst
    python -m pip install -e .
    
    
- note this package contains data in the path `atmosphtransmemullsst/data/simplegrid`
This path should be provided during runtime at initialisation.

(see in the notebooks)


### c) Installation of rubin_sim package

Package to simulate Rubin-LSST, including Filters, SED

    
    git clone https://github.com/lsst/rubin_sim.git
    
Follow the installation instructions including the corresponding data: `rubin_sim_data`


    
Moreover there is a `rubin_sim` tutorial which may be read in order to
use the  `rubin_sim` library:

    https://github.com/lsst/rubin_sim_notebooks
      
In particular, magnitude calculation from the SED and filter function and the photoelectron statistic and sky background errors. 


For the SED, I found them here

    https://github.com/rhiannonlynne/photometry_sample
    


       
       
## Notebooks
 
        
        