import os

def create_input_star_file(i, ms, t, rs, o, met):
    
    """
    i:   Inclination in deg, i = 90 is equator on
    ms:  Stellar mass in solar masses
    t:   Polar temperature of the star in K
    rs:  Stellar radius in the equatior in solar radii
    o:   Omega/Omega_Kepler
    met: Metallicity. Options: -2.5,-2.0,-1.5,-1.0,-0.5,+0.0,+0.2,+0.5. Important: use '+' sign for +0.0
    directory: path to locate the input file
    """
    
    path_input_file = '../00_codes/input.star3d'
    with open(path_input_file, 'w') as file:
        file.write('* inclination [deg, i=90 equator on]' + '\n')
        file.write('* Mstar       [solar masses]' + '\n')
        file.write('* T(pole)     [polar temperature (K)]' + '\n')
        file.write('* R(equator)  [solar radii]' + '\n')
        file.write('* omega       [Omega/Omega_Kepler]' + '\n')
        file.write('* Metallicity [-2.5,-2.0,-1.5,-1.0,-0.5,+0.0,+0.2,+0.5]' + '\n')
        file.write('*             [Important: use '+' sign for +0.0]' + '\n')
        file.write('*' + '\n')
        file.write(f'{i}' + '\n')
        file.write(f'{ms}' + '\n')
        file.write(f'{t}' + '\n')
        file.write(f'{rs}' + '\n')
        file.write(f'{o}' + '\n')
        file.write(f'{met}' + '\n')


def run_fastrot_spec(i = 5.0, ms = 2.15, t = 10000, rs = 2.726, o = 0.632, met = '-0.5'):
    create_input_star_file(i, ms, t, rs, o, met)
    os.system('../00_codes/./fastrot-spec.sh')