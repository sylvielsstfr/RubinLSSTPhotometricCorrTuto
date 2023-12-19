# RubinLSSTPhotometricCorrTuto

**Photometric Correction tutorail for virtual DESC collaboration meeting 27-Feb - 3 March 2023**

- authors : Sylvie Dagoret-Campagne & Martin Rodriguez Monroy
- affiliation : IJCLab/IN2P3/CNRS
- Creaiton date : February 2023
- Last verification : December 19 2023

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


    git clone https://github.com/LSSTDESC/getObsAtmo.git
    cd getObsAtmo
    python setup.py install
    
    

### c) Installation of rubin_sim package

Package to simulate Rubin-LSST, including Filters, SED, a lite version of the rubin_sim package

    
    git clone https://github.com:sylvielsstfr/rubinsimphot.git
    cd rubinsimphot
    pip install -e .'[dev]'

Spectro-photometric Data are provided inside this lite version of rubin_sim
    

 

       
       
## Notebooks

   see (README)[README.ipynb]
   
 
        
        