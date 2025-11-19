import numpy as np
import matplotlib.pyplot as plt

# -----------------------
# Load 3-column data file
# -----------------------
# Example filename: "data.txt"
# Format per line: Col. 1: colatitude (rad)
#                  Col. 2: radius (in units of equatorial radius)
#                  Col. 3: temperature (K)

# Adjust delimiter if needed (space, comma, tab)

data = np.loadtxt("colat-r-t.a")

colat  = data[:, 0]
radius = data[:, 1]
temp   = data[:, 2]

rmin = radius.min()
rmax = radius.max()
rmargin = 0.05*(rmax - rmin)

tmin = temp.min()
tmax = temp.max()
tmargin = 0.05*(tmax - tmin)

# --------------------------
# Create figure
# --------------------------
fig, ax1 = plt.subplots(figsize=(8, 5))

ax1.set_xlim(0.05,1.570796)
ax1.set_xlabel(r"Colatitude $\theta$ (rad)", fontsize=13)

# --- LEFT AXIS (Temperature) ---

color_temp = "red"
ax1.set_ylabel(r"Temperature (K)", color=color_temp, fontsize=13)
ax1.plot(colat, temp, color=color_temp, linewidth=1.25)
ax1.tick_params(axis="y", colors=color_temp)
ax1.set_ylim(tmin - tmargin, tmax + tmargin)      # <<< set scale here
fig.text(0.12, 0.08, r"$\theta=0$", ha='center', va='bottom', fontsize=12)
fig.text(0.12, 0.027, "Pole", ha='center', va='bottom', fontsize=12)
fig.text(0.88, 0.08, r"$\theta=\pi/2$", ha='center', va='bottom', fontsize=12)
fig.text(0.88, 0.027, "Equator", ha='center', va='bottom', fontsize=12)

# --- RIGHT AXIS (Radius) ---

ax2 = ax1.twinx()  # second y-axis
color_rad = "blue"
ax2.set_ylabel(r"R / R$_\mathrm{equatorial}$", color=color_rad, fontsize=13)
ax2.plot(colat, radius, color=color_rad, linewidth=1.25)
ax2.tick_params(axis="y", colors=color_rad)
ax2.set_ylim(rmin - rmargin, rmax + rmargin)       # <<< set scale here

plt.title("Temperature and radius vs colatitude", fontsize=14)
plt.tight_layout()
plt.show()
fig.savefig('temp-radius-vs.colatitude.pdf', dpi=300, bbox_inches='tight')
