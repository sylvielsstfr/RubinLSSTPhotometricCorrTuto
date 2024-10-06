# README.md

- author : Sylvie Dagoret-Campagne
- affiliation : IJCLab.in2p3.fr
- creation date : 2024-10-03
- last update : 2024-10-05


## Base tool to read gaia spectra  

- https://github.com/corentinravoux/gaiaspec/blob/main/gaiaspec/getGaia.py
- gaiaspec : https://github.com/corentinravoux/gaiaspec/blob/main/gaiaspec/getGaia.py
- git clone git@github.com:corentinravoux/gaiaspec.git


## notebooks

### Test Gaia access

- **TestGaia.ipynb** : ????

### Browse Gaia and CALSPEC catalog

- **BrowseGaiaCatalog.ipynb** : can work from usdf, not elasewhere probably because of acess-rights
- **BrowseGaiaCalspecCatalogs.ipynb**	

### Save spectra in hdf5 file 
- to read it back later from anywhere : **SaveGaiaSpectra_tohdf5.ipynb**
- **SaveGaiaAndCalspecSpectra_tohdf5.ipynb**


### Read spectra from hdf5 file
- **ReadGaiaSpectra_fromhdf5.ipynb**
- **ReadGaiaCalspecSpectra_fromhdf5.ipynb**	




### Smooth Gaia spectra TBD
- **SmoothGaiaSpectra.ipynb** : find the good way to interpolate and extend the gaia spectrum (2024-10-05)

### Check the calibration of spectra by computing magnitudes in LSST bands (2024-10-06)


- **CompareMagsGaiaCalspecPickleSpectra_fromhdf5.ipynb** : learn how to compute magnitudes in LSST bands, and the need to extend the gaia spectra inside LSST band range definition. 

### Comparson of Gaia spectra with Pickles SED
- **CompareMagsGaiaCalspecSpectra_fromhdf5.ipynb** : Compare flambda and fnu

### Match Gaia spectra with  Pickles
