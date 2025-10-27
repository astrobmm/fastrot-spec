# fastrot-spec
Computes the surface parameters -radius, temperature, and effective gravity as a function of the latitude- of rapidly rotating stars with radiative envelopes, and builds their synthetic spectra. The parameterisation of the stellar surface is done using the semianalytical approach of the ESTER model presented by Espinosa-Lara and Rieutord (2011); that approach avoids any discussion about the gravity darkening exponent. The computation of the synthetic spectra is done using the atlas9 and synthe suite of codes by Kurucz (2014).

## Installation & Requirements
 
Clone the above structure of folders on your computer. Your must have at the same level these five directories:

```
00_codes/
00_ld/
00_models/
00_python/
11_results/
```
Unzip the 16 zip files in folder ```00_models/``` which contain the synthetic spectra computed for eight values of the metallicity, namely [M/H]= -2.5, -2.0, -1.5, -1.0, -0.5, +0.0, +0.2, +0.5 (Kurucz's way of labelling this files has been used: m25, m20, m15, m05, p00, p02, p05 ('m', 'p' stand for 'minus', 'plus'). WARNING: Mac computers unzip zip files in two different ways: is you use ```unzip``` from the terminal, the uncompressed files will be stored in the current directory, however, if you decompress the zip file clicking on the icon, a folder with the name of the zip will be created in the current directory and the files will be stored there. Check that the models are in ```00_models/``` before running the code.

## How to run the code

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

Make sure the scripts ```fastrot-highres-spec.sh``` and ```fastrot-lowres-spec``` have  the option to be executable, type: 

```
ls -l *.sh
```

and see if you get something like ```-rwx------```,  if not, type:

```
chmod u+x *.sh
```

To run the codes, using the input files given in the directory, just type:

```
./fastrot-highres-spec.sh
```
or
```
./fastrot-lowres-spec.sh
```

And that's it! Your results will be in folder ```11_results/``` 

Another option is to use the two python scripts stored in ```00_python/```. ```run_from_python.py``` runs the whole code from the terminal, and ```functions_python.py``` allows the user
to call the code from an independent python code, e.g. if a grid of models with different input parameters needs to be built, that function would link the python code written by the user, 
which generates a set of input parameters, with the fortran codes.





A description of the input and output files can be found in ```00_codes/README_inputs_outputs.txt```
