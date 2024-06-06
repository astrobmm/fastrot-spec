DESCRIPTION OF THE OUTPUT FILES FROM EACH PROGRAM

**************
PROGRAM star3d
**************

+ Input files:

------------
input.star3d
------------
Input file containing the basic inputs to parameterise
the stellar surface. Contents are self-explanatory.
Input must be written just after the vertical bar | .

This is an example, with the parameters for Vega:

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
the vertical bar | . Values of 2.0 degrees in both
directions are recommended. Larger values can produce
numerical noise. The file looks like this:

* deltaphi    [deg, tile separation]
* deltatheta  [deg, tile separation]
2.0
2.0

+ Output files:

------------
input.spec3d [to be read out by spec3d]
------------
Small input file to be used by program spec3d. It contains
four lines specifying the metallicity of the models ([M/H]),
the name of the output file where the synthetic spectrum will
be stored, the number of latitudinal strips from pole to pole
in which the star has been divided, and the number of cells
that are visible to the observer.

Concerning the name of the output spectrum, it is composed
as follows:

Taking as example the one for Vega:

in5.0_m2.15_tpole10000_req2.726_mh-0.5_om0.632_v021.a

in:    inclination in degrees (0.0: pole on, 90.0: equator on)
m:     stellar mass (solar units)
tpole: polar temperature (K)
req:   equatorial radius (solar units)
mh:    metallicity ([M/H])
om:    ratio Omega/Omega_Kepler
v:     v sin i (km/s)

-------
tiles.a [to be read by spec3d]
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

----------
xyparams.a [output for plotting purposes]
----------
File useful if one wishes to plot several 3D quantities
of the stellar surface in 2D graphs. For each cell, it
contains:

Col. 1: x coordinate (in units of R_equator)
Col. 2: y coordinate (in units of R_equator)
Col. 3: temperature (K)
Col. 4: log g (g: effective gravity in cgs units)
Col. 5: angle between the normal to each tile and the line
        of sight
Col. 6: projected area as seen by the observer
Col. 7: angle between the normal to the differential surface
        element and the z direction (COMPUTED before the star
	is tilted by an angle i)

----------------------------------------------------------
in5.0_m2.15_tpole10000_req2.726_mh-0.5_om0.632_v021.params
----------------------------------------------------------
Contains a summary of the input parameters specified in
file input.star3d  and some of the relevant stellas parameters
provided by the omega-model. Contents are self-explanatory


**************
PROGRAM spec3d
**************

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


****************
PROGRAM normspec
****************

+ Input file: 

Composite final spectrum in physical units computed
by program spec3d (see above)


+ Output file:

------------------------------------------------------
in5.0_m2.15_tpole10000_req2.726_mh-0.5_om0.632_v021.n1
------------------------------------------------------
[Name given as example] Composite final spectrum,
continuum normalized to 1.0

Col. 1: Wavelength (angstroms)
Col. 2: Normalized intensity (continuum I=1.0)
