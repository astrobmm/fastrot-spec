-----------------------
   SYNTHETIC SPECTRA
-----------------------

This folder contains the models required to compute the composite 
spectrum of a rapidly totating star. 

- General properties of the synthetic spectra:

  Range: 360 - 550 nm
  v sini: 0 km/s
  R=lambda/Delta lambda=100000 (at lambda=450 nm)

  Col. 1: Wavelength (A)
  Col. 2: F_phot (erg cm-2 s-1 A-1)

- The files are packed in zip files with sizes < 20 Mb, and
  separated by abundances.

- This is an inventory of all available models:

   T: 6500 - 26000 K, step 100 K
   log g: 3,0, 3.5, 4.0, 4.5
   [M/H]=-2.5, .2.0, -1.5, -1.0, -0.5, 0.0, +0.2, +0.5
   
   N=((26000-6500)/100+1)*4*8    _______  6272 models

   T: 26100 - 31000 K, step 100 K
   log g: 3.5, 4.0, 4.5
   [M/H]=-2.5, .2.0, -1.5, -1.0, -0.5, 0.0, +0.2, +0.5

   N=((31000-26100)/100+1)*3*8   -------- 1200 models

   T: 31100 - 37000 K, step 100 K
   log g: 4.0, 4.5
   [M/H]=-2.5, .2.0, -1.5, -1.0, -0.5, 0.0, +0.2, +0.5

   N=((37000-31100)/100+1)*2*8   --------  960 models

                                ---------------------
                                TOTAL --- 8432 models
                                ---------------------
