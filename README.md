# fastrot-spec
Computes the surface parameters -radius, temperature, and effective gravity as a function of the latitude- of rapidly rotating stars with radiative envelopes, and builds their synthetic spectra. The parameterisation of the stellar surface is done using the semianalytical approach of the ESTER model presented by Espinosa-Lara and Rieutord (2011); that approach avoid any discussion about the gravity darkening exponent. The computation of the synthetic spectra is done using the atlas9 and synthe suite of codes by Kurucz (2014).

## Installation & Requirements
 
Clone the structure of folders on your computer. Your must have at the same level these four directories:

```
00_codes/
00_ld/
00_models/
11_results/
```
Unzip the eight zip files in folder ```00_models``` which contain the synthetic spectra computed for eight values of the metallicity, namely [M/H]= -2.5, -2.0, -1.5, -1.0, -0.5, +0.0, +0.2, +0.5 (Kurucz's way of labelling this files has been used: m25, m20, m15, m05, p00, p02, p05 ('m', 'p' stand for 'minus', 'plus'). The contents of folders ```00_ld/``` and ```00_models/``` should not be modified or touched.

## How to run the code

The codes are written in fortran77, so you must have a fortran compiler in your system. Go to folder 00_codes and compile the programs by typing:

```
gfortran -std=legacy star3d.f -o star3d.exe
gfortran -std=legacy spec3d.f -o spec3d.exe
gfortran -std=legacy normspec.f -o normspec.exe
```

Make sure the script ```fastrot-spec.sh``` has the option to be an executable, if not, type:

```
chmod u+x fastrot-spec.sh
```

To run the codes, using the input files given in the directory just type:

```
./fastrot-spec.sh
```
And that's it! Your results will be in folder ```11_results/``` 

A description of the input and output files can be found in ```00_codes/README_inputs_outputs.txt```
