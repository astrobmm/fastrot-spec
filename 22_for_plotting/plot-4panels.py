# Import packages

import numpy as np
import math
import matplotlib
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from   matplotlib import cm, colors, colorbar
from   matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                                 AutoMinorLocator)
from   matplotlib.ticker import ScalarFormatter
#

####
# Define the font size, and the type of font
#
font = {'family' : 'sans-serif', 'size' : 12}
# font = {'family' : 'sans-serif', 'weight' : 'bold', 'size' : 13}
#
plt.rc('font', **font)
plt.rc('axes', linewidth=1.25)
####

####
# Open file with the measurements
#
print(' ')
print('Reading file with these columns:')
print(' ')
print('C1: xr [Normalized to R_eq]   ')
print('C2: yr [Normalized to R_eq]   ')
print('C3: Temperature [K]           ')
print('C4: log (g [cgs]              ')
print('C5: vrot_z [km/s]             ')
print('C6: Limb darkening coefficient')

sample_file=open('4plots.a')

#
#    cmap=plt.get_cmap('rainbow')
#    label = 'Temperature [K]'
#
#    cmap=plt.get_cmap('rainbow')
#    label = 'log (g [cgs])'
#
#    cmap=plt.get_cmap('rainbow')
#    label = 'vrot component in the line of sight [km/s]'
#
#    cmap=plt.get_cmap('Reds')
#    label = 'Limb darkening coefficient'
#

####

####
#
# Skip lines with comments
datalines =sample_file.readlines()
sample_file.close()
####

####
# Initialize the variables
#
posx = []
posy = []
var1 = []
var2 = []
var3 = []
var4 = []
####

####
# Read file
#
for row in datalines:
   this_data = row.split()
   posx.append(float(this_data[0]))
   posy.append(float(this_data[1]))
   var1.append(float(this_data[2]))
   var2.append(float(this_data[3]))
   var3.append(float(this_data[4]))
   var4.append(float(this_data[5]))
####
#
# Read files with parallels and meridians
#
sample_file=open('mer10.a')
datalines =sample_file.readlines()
sample_file.close()
merx = []
mery = []
for row in datalines:
   this_data = row.split()
   merx.append(float(this_data[0]))
   mery.append(float(this_data[1]))
#
sample_file=open('par10.a')
datalines =sample_file.readlines()
sample_file.close()
parx = []
pary = []
for row in datalines:
   this_data = row.split()
   parx.append(float(this_data[0]))
   pary.append(float(this_data[1]))
####
   
####    
# Create a figure object, if both numbers are equal, the plot is squared
#
fig = plt.figure(figsize=(9., 9.))
plt.rc('axes', linewidth=1.25)
####

####
# PLOT 1: TEMPERATURE
####
ax = fig.add_subplot(2,2,1) # two rows, two columns, first plot
#
# Define the normalization for color code in terms of the temperature
norm = colors.Normalize(vmin=np.min(var1) , vmax=np.max(var1))
#
# Data points
#
cmap=plt.get_cmap('rainbow')
label = "Temperature [K]"
ax.scatter(posx, posy, marker='s', s=8., c=var1, cmap=cmap)
ax.set_aspect('equal')
#
# Meridians and parallels
### plt.plot(merx,mery, 'ro', color='black', markersize=0.1)
### plt.plot(parx,pary, 'ro', color='black', markersize=0.1)

plt.plot(merx,mery, 'o', color='black', markersize=0.1)
plt.plot(parx,pary, 'o', color='black', markersize=0.1)
#
# Colorbar, the numbers in the square refer to
# [left, bottom, width, height], where the coordinates
# are fractions that go from 0 to 1 of the plotting area. 
#
cbaxes = fig.add_axes([0.14, 0.94, 0.32, 0.02]) 
cb1 = colorbar.ColorbarBase(cbaxes, cmap=cmap,
                                    norm=norm,
                                    orientation='horizontal')
#
# Define labels for axes, number of ticks, limits, etc
#
cb1.set_label(label,fontsize=10)
cb1.ax.tick_params(labelsize=8)
ax.set_ylabel(r'Y')
ax.set_xlim([-1.3, 1.3])
ax.set_ylim([-1.3, 1.3])
ax.tick_params(which='both', width=1, direction="in")
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.locator_params(axis = 'x', nbins = 6)
ax.locator_params(axis = 'y', nbins = 6)
####

####
# PLOT 2: log g
####
ax = fig.add_subplot(2,2,2) # two rows, two columns, plot 2
#
# Define the normalization for color code 
norm = colors.Normalize(vmin=np.min(var2) , vmax=np.max(var2))
#
# Data points
#
cmap=plt.get_cmap('rainbow')
label = "log (g [cgs])"
ax.scatter(posx, posy, marker='s', s=8., c=var2, cmap=cmap)
ax.set_aspect('equal')
#
# Meridians and parallels
plt.plot(merx,mery, 'o', color='black', markersize=0.1)
plt.plot(parx,pary, 'o', color='black', markersize=0.1)
#
# Colorbar, the numbers in the square refer to
# [left, bottom, width, height], where the coordinates
# are fractions that go from 0 to 1 of the plotting area. 
#
cbaxes = fig.add_axes([0.56, 0.94, 0.32, 0.02]) 
cb1 = colorbar.ColorbarBase(cbaxes, cmap=cmap,
                                    norm=norm,
                                    orientation='horizontal')
#
# Define labels for axes, number of ticks, limits, etc
#
cb1.set_label(label,fontsize=10)
cb1.ax.tick_params(labelsize=8)
ax.set_xlim([-1.3, 1.3])
ax.set_ylim([-1.3, 1.3])
ax.tick_params(which='both', width=1, direction="in")
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.locator_params(axis = 'x', nbins = 6)
ax.locator_params(axis = 'y', nbins = 6)
####

####
# PLOT 3: V SIN I
####
ax = fig.add_subplot(2,2,3) # two rows, two columns, plot 3
#
# Define the normalization for color code
norm = colors.Normalize(vmin=np.min(var3) , vmax=np.max(var3))
#
# Data points
#
cmap=plt.get_cmap('rainbow')
label = r"v$_{\mathrm{rot}}$ component in the line of sight [km/s]"
ax.scatter(posx, posy, marker='s', s=8., c=var3, cmap=cmap)
ax.set_aspect('equal')
#
# Meridians and parallels
plt.plot(merx,mery, 'o', color='black', markersize=0.1)
plt.plot(parx,pary, 'o', color='black', markersize=0.1)
#
# Colorbar, the numbers in the square refer to
# [left, bottom, width, height], where the coordinates
# are fractions that go from 0 to 1 of the plotting area. 
#
cbaxes = fig.add_axes([0.14, 0.49, 0.32, 0.02]) 
cb1 = colorbar.ColorbarBase(cbaxes, cmap=cmap,
                                    norm=norm,
                                    orientation='horizontal')
#
# Define labels for axes, number of ticks, limits, etc
#
cb1.set_label(label,fontsize=10)
cb1.ax.tick_params(labelsize=8)
ax.set_xlabel(r'X')
ax.set_ylabel(r'Y')
ax.set_xlim([-1.3, 1.3])
ax.set_ylim([-1.3, 1.3])
ax.tick_params(which='both', width=1, direction="in")
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.locator_params(axis = 'x', nbins = 6)
ax.locator_params(axis = 'y', nbins = 6)
####

####
# PLOT 4: LIMB DARKENING COEFFICIENT 
####
ax = fig.add_subplot(2,2,4) # two rows, two columns, plot 4
#
# Define the normalization for color code
norm = colors.Normalize(vmin=np.min(var4) , vmax=np.max(var4))
#
# Data points
#
cmap=plt.get_cmap('Reds')
label = r"Limb darkening coefficient at 5000 $\mathrm{\AA}$"
ax.scatter(posx, posy, marker='s', s=8., c=var4, cmap=cmap)
ax.set_aspect('equal')
#
# Meridians and parallels
plt.plot(merx,mery, 'o', color='white', markersize=0.1)
plt.plot(parx,pary, 'o', color='white', markersize=0.1)
# Colorbar, the numbers in the square refer to
# [left, bottom, width, height], where the coordinates
# are fractions that go from 0 to 1 of the plotting area. 
#
cbaxes = fig.add_axes([0.56, 0.49, 0.32, 0.02]) 
cb1 = colorbar.ColorbarBase(cbaxes, cmap=cmap,
                                    norm=norm,
                                    orientation='horizontal')
#
# Define labels for axes, number of ticks, limits, etc
#
cb1.set_label(label,fontsize=10)
cb1.ax.tick_params(labelsize=8)
ax.set_xlabel(r'X')
ax.set_ylabel(r'Y')
ax.set_xlim([-1.3, 1.3])
ax.set_ylim([-1.3, 1.3])
ax.tick_params(which='both', width=1, direction="in")
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.locator_params(axis = 'x', nbins = 6)
ax.locator_params(axis = 'y', nbins = 6)
####

fig.subplots_adjust(hspace=0.45)

####
# When the plot appears in the display one can stretch it, and when closed,
# the pdf file built by the program will have the aspect ratio of the
# modified plot
#
plt.show()
fig.savefig('4plots.pdf', dpi=300, bbox_inches='tight')
####
