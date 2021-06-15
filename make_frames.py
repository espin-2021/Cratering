"""Make frames for gif of time-evolution."""

import os
import numpy as np
from crater_functions import do_cratering, make_noisy_surface
import matplotlib.pyplot as plt

# fix random seed
np.random.seed(1023)

# Set the colourmap
cmap = "Greys_r"  # old craters in gray (reversed color map)
cmap1 = "Greys" # (Greys, as is - i.e. not reversed)
cmap2 = "Spectral_r" # (Spectral, reversed)

# Define some Variables!
size = 1000
spacing = 20

# Set up a LandLab model grid object, random topographic noise ~O(zfactor)
xy = int(size/spacing)
mg = make_noisy_surface(xy,spacing,rf=1) # rf is randomness factor

# Set the "distribution" as weights, and choose a value from those weights
Kx = 1.0  # Scaling coefficient
delta = 2.0  # km, scaling exponent
Ncraters = 2 # Number of craters to add
#Ncraters = int(input("How many craters do you want to add per timestep? ")) # see console to add your own number!
minD = int(spacing * 3) # want it to be greater than spacing, three times the smallest cell
print("Minimum possible crater size:", minD, "km")
maxD = int((size)/4) #max diameter (m), a quarter of the total domain width
print("Maximum possible crater size:", maxD, "km")
NDs = []
for D in range(minD, maxD):
    ND = Kx * D ** - delta
    NDs.append(ND)

# make output folder if needed
if os.path.isdir('figs') is False:
    os.mkdir('figs')

# Define where section is taken along y (in km)
stk_km = 500
stk = int(stk_km/spacing) #node location of stk_km
# Number of timesteps
nsteps = 10

# set vertical limits of plots (should be based on last timestep elevations)
zmin = -100
zmax = 20


old_arr = np.zeros((nsteps, xy))

for i in range(1, nsteps):
    # crater landscape
    mg = do_cratering(Ncraters, NDs, minD, maxD, xy, mg, spacing)

    # update plot
    fig, ax = plt.subplots(1, 2, dpi=250, facecolor='w', figsize=(8, 3))
    topo = mg.field_values('node', 'topographic__elevation').reshape((xy, xy)) #create an array of elevation values across the mg domain
    hs = mg.calc_hillshade_at_node(elevs='topographic__elevation') #create a hillshade array
    hill = np.reshape(hs, (xy, xy)) #reshape the hillshade array to be the same shape as the topo array
    img1 = ax[0].imshow(hill, cmap=cmap1, alpha=1, extent = [0,size, 0, size], vmin=0, vmax=1)
    img2 = ax[0].imshow(topo, cmap=cmap2, alpha=0.6, extent = [0,size, 0, size], vmin=zmin, vmax=zmax)
    ax[0].set_title("Cratered surface evolution (t = %i)" %i)
    cbar = plt.colorbar(img2, fraction=0.045, ax=ax[0], label = "Elevation [m??????]")
    ax[0].set_xlabel('X [km]')
    ax[0].set_ylabel('Y [km] ')
    #stk = int(stk_km/spacing)
    ax[0].plot(np.linspace(0, size), np.ones_like(np.linspace(0, size))*stk_km,
               c='k', linestyle='--')
    ax[0].set_xlim([0, size])
    ax[0].set_ylim([0, size])

    ax[1].plot(topo[stk, :], c=[0, 0, 0], zorder=10)
    for j in range(old_arr.shape[0]):
        if np.sum(old_arr[j, :]) != 0:
            ax[1].plot(old_arr[j, :], c=[1-j/nsteps, 1-j/nsteps, 1-j/nsteps])
    ax[1].set_title('Topographic section at Y = ' + str(stk_km) + " km")
    ax[1].set_ylabel('Elevation [m???]')
    ax[1].set_xlabel('Distance along X')
    ax[1].set_ylim() # I think making a dynamic scale is kind of fun, shows how much it changes? Instead of [zmin, zmax]

    plt.tight_layout()
    plt.savefig('figs/' + '0'.zfill(4) + '.png', bbox_inches='tight')    

    # retain old array section
    old_arr[i, :] = topo[stk, :]
