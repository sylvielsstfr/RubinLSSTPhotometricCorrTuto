# README.md : 04_FitSED

- author : Sylvie Dagoret-Campagne
- affiliation : IJCLab/IN2P3/CNRS
- creation date : 2024-01-18
- update : 2024-02-25

## purpose

Find a method to estimate the SED from magnitudes or integrated flux in a known band.
We propose to use gaussian process to interpolate the SED function between point-fluxes in bands.
The point-fluxes has to be fit from observed or catalog magnitudes or in-band fluxes.


- Note that a method to add unknown extra-points outside LSST filter range is recommented to stabilize
the SED fitted by the Gaussian Process.


- Note that the fit of the SED through the SED points allow to calculate the derivative if the SED usefull to compute the color correction.

## Notebook


- **04a_FitSEDshape.ipynb** : Just example showing how to access to SED, no fit
                      
- **04b_FitSEDshape-GaussianProcess.ipynb** : Fit using Gaussian Process Regression (scikit-learn)

- **04c_FitSEDshape-GaussianKernelRegression.ipynb** : Fit using Gaussian Basis Function Regression (astroml)
- 
- **04d_FitSEDshape-GaussianBasisRegression.ipynb** : Fit using Gaussian Kernel Regression (astroml)
