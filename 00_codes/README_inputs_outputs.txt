------------------------------------------------------------
  DESCRIPTION OF THE PROGRAMS, INPUTS AND RELEVANT OUTPUTS
------------------------------------------------------------

N.B: Some auxiliary files are created by the programs, but
are not relevant regarding the outputs of interest to the 
user, therefore are not listed here.

**************
PROGRAM star3d
**************

+ Input files:

------------
input.star3d
------------
Input file containing the basic inputs to parameterise
the stellar surface. Contents are self-explanatory.
Input must be written after the lines started by
an asterisk.

This is an example:

* inclination [deg, i=90 equator on]    
* Mstar       [solar masses]            
* T(pole)     [polar temperature (K)]    
* R(equator)  [solar radii]              
* omega       [Omega/Omega_Kepler]      
* Metallicity [-2.5,-2.0,-1.5,-1.0,-0.5,+0.0,+0.2,+0.5]        
*             [Important: use '+' sign for +0.0] 
*
5.0     
2.15    
10000   
2.726   
0.632   
-0.5    

------------
input.grid3d
------------
Input file containing the information on the tile
separation in latitude and longitude. Contents are
self explanatory. Input must be written just after
the lines started by an asterisk. Values of 2.0
degrees in both directions are recommended. Larger
values could produce numerical noise. The file looks
like this:

* deltaphi    [deg, tile separation]
* deltatheta  [deg, tile separation]
2.0
2.0

+ Output files:

------------
input.spec3d [to be read by both highrespec3d and lowrespec3d]
------------

Input file to be used by programs

highrespec3d
lowrespec3d

It contains four lines specifying the metallicity of the models
([M/H]), the name of the output file where the synthetic spectrum
will be stored, the number of latitudinal strips from pole to pole
in which the star has been divided, and the number of cells that
are visible to the observer.

The name of the output spectrum from each code is composed as
follows: Taking as example the one shown above

in5.0_m2.15_tpole10000_req2.726_mh-0.5_om0.632_v021.a

in:    inclination in degrees (0.0: pole on, 90.0: equator on)
m:     stellar mass (solar units)
tpole: polar temperature (K)
req:   equatorial radius (solar units)
mh:    metallicity ([M/H])
om:    ratio Omega/Omega_Kepler
v:     v sin i (km/s)

This will be the name of the output file for highrespec3d
(see below).

lowrespec3d will give as output:

low-res-in5.0_m2.15_tpole10000_req2.726_mh-0.5_om0.632_v021.a

The program itself takes the name from input.spec3d and adds
the prefix 'low-res-' (see below)

-------
tiles.a [to be read by both highrespec3d and lowrespec3d]
-------
For each visible cell it contains:

Col. 1: x coordinate (in units of R_equator)
Col. 2: y coordinate (in units of R_equator)
Col. 3: latitude (degrees)
Col. 4: longitude (degrees)
Col. 5: projected area as seen by the observer
Col. 6: radial velocity as seen by the observer (km/s)
Col. 7: cosine of the angle between the normal to each
        cell and the line of sight
	
-----------
colat-r-t.a [output for plotting purposes]
-----------
Col. 1: colatitude (rad)
Col. 2: stellar radius (normalized to R_equatorial)
Col. 3: photospheric temperature (K)

----------------------------------------------------------
in5.0_m2.15_tpole10000_req2.726_mh-0.5_om0.632_v021.params
----------------------------------------------------------

Contains a summary of the input parameters specified in
file input.star3d  and some of the relevant stellas parameters
provided by the omega-model. Contents are self-explanatory

WARNING: star3d gives an estimate of L_star computed by adding
the values of the flux for each cell computed under the
approximation sigma*T**4. Use program lowrespec3d for a
more exact calculation (see below).

********************
PROGRAM highrespec3d
********************

This program computes the high resolution (R=100,000)
spectrum of the star in the interval 3600-5500 A, according
to the input parameters specified in input.star3d and
input.grid3d

+ Input file:

------------
input.spec3d
------------
File created by program star3d, see above

+ Output file:

-----------------------------------------------------
in5.0_m2.15_tpole10000_req2.726_mh-0.5_om0.632_v021.a
-----------------------------------------------------

[Name given as example] Composite final spectrum
of the oblate- rapidly-rotating star:

Col. 1: Wavelength (angstroms)
Col. 2: Flux (erg cm-1 s-1 A-1) at the stellar surface

********************
PROGRAM lowrespec3d
********************

This program computes the low resolution spectrum of the
star in the interval covered by the Castelli-Kurucz low
resolution synthetic models, 90.0-1.6e+6 A. The resolution
is variable accross the spectrum. The input parameters
specified in input.star3d and input.grid3d are used.

+ Input file:

------------
input.spec3d
------------
File created by program star3d, see above

+ Output file:

-------------------------------------------------------------
low-res-in5.0_m2.15_tpole10000_req2.726_mh-0.5_om0.632_v021.a
-------------------------------------------------------------

[Name given as example] Composite final spectrum
of the oblate- rapidly-rotating star:

Col. 1: Wavelength (angstroms)
Col. 2: Flux (erg cm-1 s-1 A-1) at the stellar surface

This program must be used to compute in a more exact manner
the stellar luminosity. As mentioned above, star3d gives
an estimate of L_star computed by adding the values of the
flux for each cell under the approximation that the integrated
flux is sigma*T**4. Instead, lowrespec3d uses the corresponding
Kurucz synthetic spectrum for finding the total output
flux from each cell.

****************
PROGRAM normspec
****************

+ Input file: 

Composite final spectrum in physical units computed
by program highrespec3d (see above)

+ Output file:

Composite final spectrum, continuum normalized to 1.0

------------------------------------------------------
in5.0_m2.15_tpole10000_req2.726_mh-0.5_om0.632_v021.n1
------------------------------------------------------
[Name given as example] 

Col. 1: Wavelength (angstroms)
Col. 2: Normalized intensity (continuum I=1.0)

**************
PROGRAM grid10
**************

Creates a grid of parallels and meridians separared by
10 degrees. It reads the relevant parameters from input.star3d
and gives as outputs the files mer10.a and par10.a that are 
used for plotting putposes (see the information and python 
program plot-4panels.py in folder 22_for_plotting/)

******************
PROGRAM mergefiles
******************

Merges the contents of the file xyparams.a, created by star3d
and limbdark.a created by either lowrespec3d and highrespec3d
and creates a file called 4plots.a used for plotting purposes
(see the information and python program plot-4panels.py in 
folder 22_for_plotting/). The file contains six columns:

Col. 1: x coordinate (in units of R_equator)
Col. 2: y coordinate (in units of R_equator)
Col. 3: temperature (K)
Col. 4: log g (g: effective gravity in cgs units)
Col. 5: projected velocity of each cell as seen by 
        the observer (km/s)
Col. 6: Limb darkening coefficient at 5000 A

