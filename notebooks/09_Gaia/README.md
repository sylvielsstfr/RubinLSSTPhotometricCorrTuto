# README.md

- author : Sylvie Dagoret-Campagne
- affiliation : IJCLab.in2p3.fr
- creation date : 2024-10-03
- last update : 2024-10-07


## Base tool to read gaia spectra  

- https://github.com/corentinravoux/gaiaspec/blob/main/gaiaspec/getGaia.py
- gaiaspec : https://github.com/corentinravoux/gaiaspec/blob/main/gaiaspec/getGaia.py
- git clone git@github.com:corentinravoux/gaiaspec.git


## notebooks

### Test Gaia access

- **TestGaia.ipynb** : ????

### Browse Gaia and CALSPEC catalog to get some understanding of how to access the info

- **BrowseGaiaCatalog.ipynb** : can work from usdf, not elsewhere probably because of acess-rights
- **BrowseGaiaCalspecCatalogs.ipynb**	 : can work from usdf, not elsewhere probably because of acess-rights.

### Save spectra in hdf5 file.
Those notebooks produce a single hdf5 file with all spectra. Thoe spectra are extracted from gaiaspec package (Corentin Ravoux). Then these spectra in hdf5 can be read from anywhere. Note Spectra are converted in unit erg/cm2/nm/s. 
- **SaveGaiaSpectra_tohdf5.ipynb** : Save only Gaia spectra in an hdf5 file.
- **SaveGaiaAndCalspecSpectra_tohdf5.ipynb** : Save both Gaia and Calspec spectra in HDF5 file.


### Read spectra from hdf5 file.
The HDF5 file written in  **SaveGaiaSpectra_tohdf5.ipynb** and **SaveGaiaAndCalspecSpectra_tohdf5.ipynb** are read-back. The advantage is that we don't need Gaia access read.
The Spectra are plotted in unit of erg/cm2/nm/s. 
- **ReadGaiaSpectra_fromhdf5.ipynb** : View Gaia spectra Only. Spectra are red from hdf5 file.
- **ReadGaiaCalspecSpectra_fromhdf5.ipynb**	 : View all Gaia Spectra and compare it with its calspec correspondant. File are red from hdf5 file.




### Smooth Gaia spectra (2024-10-05)
Exercice to check how the interpolation and extrapolation of Gaia Spectra can be performed.
Moreover it shows how to smooth a Gaia Spectra given a filtering window.
It is for later use namely when calculating the magnitudes.
- **SmoothGaiaSpectra.ipynb** : find the good way to interpolate and extend the gaia spectrum, including smoothing.

### Check the calibration of spectra by computing magnitudes in LSST bands (2024-10-06)
 Learn how to compute magnitudes in LSST bands from Spectra:
-**CompareMagsGaiaCalspecSpectra_fromhdf5.ipynb** : Show the histogram of magnitudes of Gaia spectra in the LSST bands. And it computes the magnitude difference between Gaia and Calspec Spectra.
Note it is needed to extend the gaia spectra spectrum at its borders inside the whole wavelength LSST band range definition. 

### Comparison of Gaia spectra with Pickles SED (2024-10-06)
-**CompareMagsGaiaCalspecPickleSpectra_fromhdf5.ipynb** : Compare flambda and fnu of one Gaia to Extreme SED pickles colors. It needs to be developed or corrected. Better work in next notebook.

### Match Gaia spectra with  Pickles by using nearest neighbourg (2024-10-06)
The goal is to find the neared Pickle SED to a Gaia Spectrum. The nearest Neighbourg is done after renormalizing Pickle Spectra to match Z band magnitude to that of Gaia. Then the Nearest Neighbourg is done by Matching the magnitudes in G,R,I (ad Z) by definition.

-**FindNearestNeighbourgsGaiaCalspecPickleSpectra_fromhdf5.ipynb**: Match on one Gaia Spectrum.
-**LoopFindNearestNeighbourgsGaiaCalspecPickleSpectra_fromhdf5.ipynb**: Match on all Gaia Spectra.
