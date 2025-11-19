import numpy as np
import matplotlib.pyplot as plt

# -----------------------
# Load 2-column data file
# -----------------------

# File *.a : Col. 1: lambda (A)
#            Col. 2: Flux (erg cm-2 s-1 A-1)
#
# File *.n1: Col. 1: lambda (A)
#            Col. 2: Normalized intensity

data1 = np.loadtxt("in40.0_m2.15_tpole10000_req2.726_mh-0.5_om0.632_v158.a")
data2 = np.loadtxt("in40.0_m2.15_tpole10000_req2.726_mh-0.5_om0.632_v158.n1")

wave1 = data1[:, 0]
flux1 = data1[:, 1]
wave2 = data2[:, 0]
flux2 = data2[:, 1]

fmin1 = flux1.min()
fmax1 = flux1.max()
fmargin1 = 0.01*(fmax1 - fmin1)

fmin2 = flux2.min()
fmax2 = flux2.max()
fmargin2 = 0.01*(fmax2 - fmin2)

####    
# Create a figure object, if both numbers are equal, the plot is squared
#
fig = plt.figure(figsize=(7., 7.))
plt.rc('axes', linewidth=1.25)
####

####
# TOP: Spectrum in physical units 
####

ax = plt.subplot(2,1,1)
ax.set_xlim(3600.0,5500.0)
ax.set_ylim(fmin1 - fmargin1, fmax1 + fmargin1)
ax.set_ylabel(r"Flux (erg cm$^{-2}$ s$^{-1}$ $\mathrm{\AA}^{-1}$)", fontsize=12)
ax.plot(wave1, flux1, linewidth=1.0)
ax.tick_params(axis="y")

####
# BOTTOM: Normalized spectrum 
####

ax = plt.subplot(2,1,2)
ax.set_xlim(3600.0,5500.0)
ax.set_ylim(fmin2 - fmargin2, 1.05)
ax.set_xlabel(r"Wavelength ($\mathrm{\AA}$)", fontsize=12)
ax.set_ylabel(r"Normalized flux", fontsize=12)
ax.plot(wave2, flux2, linewidth=1.0)
ax.tick_params(axis="y")

plt.tight_layout()
plt.show()
fig.savefig('highres-spectrum.pdf', dpi=300, bbox_inches='tight')
