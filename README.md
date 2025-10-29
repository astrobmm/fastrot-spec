# fastrot-spec
This suite of codes computes the surface parameters -radius, temperature, and effective gravity as a function of the latitude- of rapidly rotating stars with radiative envelopes, and builds their synthetic spectra. 

The parameterisation of the stellar surface is done using the semianalytical approach of the ESTER model presented by Espinosa-Lara and Rieutord (2011); that approach avoids any discussion about the gravity darkening exponent. The computation of the synthetic spectrum is done using the ```atlas9``` and ```synthe``` suite of codes by Kurucz (2014). Three main codes are provided in folder ```00_codes/```, namely ```star3d```, ```highrespec3d```, and ```lowrespec3d```. ```star3d``` computes the parameterization of the stellar surface that will be used by the other two programs to compute the high-resolution (3600-5500 A, R=100,000) and low-resolution (90-1.0e+6 A, variable with wavelength, according to the Castelli-Kurucz low resolution models). 

This procedure is described in the paper "Surface parameterisation and spectral synthesis of rapidly rotating stars. Vega as a test bed" (Montesinos, B., 2024, Astronomy and Astrophysics, 688, A97). 

October 2025: Updates and improvements have been added with respect to the first version of the codes, made publicly available in 2024. Use the new versions of the codes and scripts in case the original ones are used.

  - Program ```spec3d``` has been splitted into ```highrespec3d and ```lowrespec3d```.
  - Program ```star3d``` includes a small modification required for the new code ```lowrespec3d```.
  - The computation of the limb darkening coefficients has been improved, attempting to gain speed in the calculations. 

## Installation & Requirements
 
Clone the following structure of folders on your computer. Your must have at the same level these five directories:

```
00_codes/
00_ld/
00_models/
00_pckfiles/
11_results/
```
Download and unzip the 16 zip files in folder ```00_models/``` which contain the high-resolution synthetic spectra computed for eight values of the metallicity, namely [M/H]= -2.5, -2.0, -1.5, -1.0, -0.5, +0.0, +0.2, +0.5 (Kurucz's way of labelling this files has been used: m25, m20, m15, m05, p00, p02, p05 ('m', 'p' stand for 'minus', 'plus'). These are used to compute the emergent high-resolution spectrum of the target object using the programme ```highrespec3d```.WARNING: Mac computers unzip zip files in two different ways: is you use ```unzip``` from the terminal, the uncompressed files will be stored in the current directory, however, if you decompress the zip file clicking on the icon, a folder with the name of the zip will be created in the current directory and the files will be stored there. Check that the models are in ```00_models/``` before running the code.

Download and unzip the two zip files in folder ```00_pckfiles/``` which contain the collection of low-resolution Castelli-Kurucz models for the eigh metallicities specified in the above paragraph. These models are used to compute the low-resolution spectrum using the program ```lowrespec3d```.

## How to run the codes

The codes are written in fortran77, so you must have a fortran compiler in your system. Go to folder ```00_codes/``` and compile the programs; in linux platforms this works:

```
gfortran -std=legacy star3d.f -o star3d.exe
gfortran -std=legacy highrespec3d.f -o highrespec3d.exe
gfortran -std=legacy lowrespec3d.f -o lowrespec3d.exe
gfortran -std=legacy normspec.f -o normspec.exe
```

Some tests were made on Macs laptops, a fortran compiler such is the one in 

[https://github.com/fxcoudert/gfortran-for-macOS/releases]

seems to work, using the same commmand as above.

Make sure the scripts ```fastrot-highres-spec.sh``` and ```fastrot-lowres-spec``` have  the option to be executable, to check that type: 

```
ls -l *.sh
```

and see if you get something like ```-rwx------```,  if that is not the case, type:

```
chmod u+x *.sh
```

To run the codes, use the input files given in the directory ```00_codes/```, just type:

```
./fastrot-highres-spec.sh
```
or
```
./fastrot-lowres-spec.sh
```

And that's it! Your results will be in folder ```11_results/``` 

Another option is to use the two python scripts stored in ```00_codes/```. ```run_from_python.py``` runs the whole code from the terminal, and ```functions_python.py``` allows the user to call the code from an independent python code, e.g. if a grid of models with different input parameters needs to be built, that function would link the python code written by the user, which generates a set of input parameters, with the fortran codes.

A description of the input and output files can be found in ```00_codes/README_inputs_outputs.txt```
