# RubinLSSTPhotometricCorrTuto

**Photometric Correction tutorail for virtual DESC collaboration meeting 27-Feb - 3 March 2023**

- authors : Sylvie Dagoret-Campagne & Martin Rodriguez Monroy
- affiliation : IJCLab/IN2P3/CNRS
- Last verification : February 2021

## Goal

Understand the impact of atmospheric transparency variations on photometry, the necessity to estimate those variations and the necessity to correct them up to some level.

###1. Estimation of atmospheric variation by measuring spectra fro calibration stars with the Rubin-LSST auxiliary telescope

###2. LSST Flux corrections

###3. Impact on science
- example Photoz 

## Program


![Workflow][def]


[def]: workflow_Feb152023.jpeg


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


    




       
       
## Notebooks
 
        
        