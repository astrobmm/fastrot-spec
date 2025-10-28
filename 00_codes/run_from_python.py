import os
import argparse

#-----------------------------------------#
# Inputs to create the file "input.star3d"
#-----------------------------------------#
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inclination", default = 5.0, type = float, help = "Inclination in deg, i = 90 is equator on")
parser.add_argument("-ms", "--mstar", default = 2.15, type = float, help = "Stellar mass in solar masses")
parser.add_argument("-t", "--polar_temperature", default = 10000, type = float, help = "Polar temperature of the star in K")
parser.add_argument("-rs", "--rstar", default = 2.726, type = float, help = "Stellar radius in the equatior in solar radii")
parser.add_argument("-o", "--omega", default = 0.632, type = float, help = "Omega/Omega_Kepler")
parser.add_argument("-met", "--metallicity", default = '-0.5', type = str, help = "Metallicity. Options: -2.5,-2.0,-1.5,-1.0,-0.5,+0.0,+0.2,+0.5. Important: use '+' sign for +0.0")
args = parser.parse_args()

path_input_file = 'input.star3d'
with open(path_input_file, 'w') as file:
    file.write('* inclination [deg, i=90 equator on]' + '\n')
    file.write('* Mstar       [solar masses]' + '\n')
    file.write('* T(pole)     [polar temperature (K)]' + '\n')
    file.write('* R(equator)  [solar radii]' + '\n')
    file.write('* omega       [Omega/Omega_Kepler]' + '\n')
    file.write('* Metallicity [-2.5,-2.0,-1.5,-1.0,-0.5,+0.0,+0.2,+0.5]' + '\n')
    file.write('*             [Important: use '+' sign for +0.0]' + '\n')
    file.write('*' + '\n')
    file.write(f'{args.inclination}' + '\n')
    file.write(f'{args.mstar}' + '\n')
    file.write(f'{args.polar_temperature}' + '\n')
    file.write(f'{args.rstar}' + '\n')
    file.write(f'{args.omega}' + '\n')
    file.write(f'{args.metallicity}' + '\n')

#-----------------------------------------#
# Inputs to create the file "input.grid3d"
#-----------------------------------------#
parser = argparse.ArgumentParser()
parser.add_argument("-dphi", "--deltaphi",   default = 2.0, type = float, help = "Tile separation in longitude, in deg")
parser.add_argument("-dthe", "--deltatheta", default = 2.0, type = float, help = "Tile separation in latitude, in deg")
args = parser.parse_args()    

path_input_file = 'input.grid3d'
with open(path_input_file, 'w') as file:
    file.write('* deltaphi    [deg, tile separation]' + '\n')
    file.write('* deltatheta  [deg, tile separation]' + '\n')
    file.write(f'{args.deltaphi}' + '\n'}
    file.write/f'{args.deltatheta}' + '\n'}
    
#----------------------#
# Run the fortran codes
#----------------------#
os.system('./fastrot-highres-spec.sh')
os.system('./fastrot-lowres-spec.sh')
