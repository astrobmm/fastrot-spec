                -------------------------
                GRID OF SYNTHETIC SPECTRA
                -------------------------
______________________
CONTENTS OF THE FOLDER

This folder contains the models required to compute the composite 
spectrum of a rapidly rotating star. 

GENERAL PROPERTIES OF THE MODELS

  Range: 360 - 550 nm
  v sini: 0 km/s
  R=lambda/Delta lambda=100000 (at lambda=450 nm)

  Col. 1: Wavelength (A)
  Col. 2: F_phot (erg cm-2 s-1 A-1)

_________
ZIP FILES

  The models are packed in zip files with sizes < 20 Mb, and
  separated by abundances. Due to the github limits to the file
  size, for each abundance the models were stored in several 
  zip files named:

  models_[abundance code]_part[NN].zip

  where the abundance codes m25, m20, m15, m10, m05, p00, p02, p05
  correspond to [M/H]=-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, +0.2, +0.5

  NN goes from 01-17 (m25, m20, m15, m10, m05, p00) and 
  01-18 (p02, p05)

  When unpacking the models, the names of the individual files are

  t[TTTTT]g[GGG][AAA]_360_550_v0.m

  TTTTT is the temperature (e.g. 06500, 23200)
  GGG is log g multiplied by 100 (e.g. 300, 450)
  AAA is the abundance code as above

  An example of filename would be: t19400g450m25_360_550_v0.m
_________
INVENTORY

  This is the inventory of all available models:

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
