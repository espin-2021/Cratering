"""Make frames for gif of time-evolution."""

import os
import numpy as np
from landlab import RasterModelGrid, NodeStatus, values
from crater_functions import weighted_choice_sub, crater_depth, do_cratering
import matplotlib.pyplot as plt

# fix random seed
np.random.seed(1023)

# Set the colourmap
cmap = "Greys_r"  # old craters in gray

# Define some Variables!
xy = 200  # Set the number of nodes in both x and y space
spacing = 1
zfactor = 10  # set factor to multiply random noise by

# Set up a LandLab model grid object, random topographic noise ~O(zfactor)
mg = RasterModelGrid((xy, xy), xy_spacing=spacing)  # see above for variables
# create an array of zeros for each node of the model grid
z = mg.add_zeros('topographic__elevation', at='node')
# add noise to the surface
noise = values.random(mg, "topographic__elevation", at='node',
                      where=NodeStatus.CORE, distribution='uniform',
                      high=z.max()*0.75, low=z.max()*0.25)
mg.at_node["topographic__elevation"] *= zfactor

# Set the "distribution" as weights, and choose a value from those weights
Kx = 1.0  # Scaling coefficient
delta = 2.0  # km, scaling exponent
Ncraters = 100  # Number of craters to add
# int(spacing*3) # min diameter (m), xy_spacing times 3
# i.e. more tha 3 cells in diameter
minD = 1
# int((mg.number_of_nodes)/4) #max diameter (m)
# i,e, a quarter of the total domain width
maxD = 100
NDs = []
for D in range(minD, maxD):
    ND = Kx * D ** - delta
    NDs.append(ND)

# make output folder if needed
if os.path.isdir('figs') is False:
    os.mkdir('figs')

# define where section is taken along y
stk = 100

# number of timesteps
nsteps = 10

# set vertical limits of plots (should be based on last timestep elevations)
zmin = -100
zmax = 20

# loop to make time-series frames
fig, ax = plt.subplots(1, 2, dpi=250, facecolor='w', figsize=(8, 3))
_arr = mg.field_values('node', 'topographic__elevation').reshape((xy, xy))
img = ax[0].imshow(_arr, cmap=cmap, vmin=zmin, vmax=zmax)
ax[0].set_title('Topography')
cbar = plt.colorbar(img, fraction=0.045, ax=ax[0])
cbar.set_label('Elevation [m]')
ax[0].set_xlabel('X')
ax[0].set_ylabel('Y')
ax[0].plot(np.linspace(0, xy), np.ones_like(np.linspace(0, xy))*stk,
           c='r', linestyle='--')
ax[0].set_xlim([0, xy])
ax[0].set_ylim([0, xy])

ax[1].plot(_arr[stk, :], c=[0, 0, 0])
ax[1].set_title('Topographic Section at Y = ' + str(stk))
ax[1].set_ylabel('Topography [m]')
ax[1].set_xlabel('Distance along X')
ax[1].set_ylim([zmin, zmax])

plt.tight_layout()
plt.savefig('figs/' + '0'.zfill(4) + '.png', bbox_inches='tight')

# store old topo sections
old_arr = np.zeros((nsteps, xy))
old_arr[0, :] = _arr[:, stk]

for i in range(1, nsteps):
    # crater landscape
    mg = do_cratering(Ncraters, NDs, minD, maxD, xy, mg)

    # update plot
    fig, ax = plt.subplots(1, 2, dpi=250, facecolor='w', figsize=(8, 3))
    _arr = mg.field_values('node', 'topographic__elevation').reshape((xy, xy))
    img = ax[0].imshow(_arr, cmap=cmap, vmin=zmin, vmax=zmax)
    ax[0].set_title('Topography after ' + str(Ncraters*i) + ' impacts')
    cbar = plt.colorbar(img, fraction=0.045, ax=ax[0])
    cbar.set_label('Elevation [m]')
    ax[0].set_xlabel('X')
    ax[0].set_ylabel('Y')
    ax[0].plot(np.linspace(0, xy), np.ones_like(np.linspace(0, xy))*stk,
               c='r', linestyle='--')
    ax[0].set_xlim([0, xy])
    ax[0].set_ylim([0, xy])

    ax[1].plot(_arr[stk, :], c=[0, 0, 0], zorder=10)
    for j in range(old_arr.shape[0]):
        if np.sum(old_arr[j, :]) != 0:
            ax[1].plot(old_arr[j, :], c=[1-j/nsteps, 1-j/nsteps, 1-j/nsteps])
    ax[1].set_title('Topographic Section at Y = ' + str(stk))
    ax[1].set_ylabel('Topography [m]')
    ax[1].set_xlabel('Distance along X')
    ax[1].legend(['New Topo', 'Old Topo'])
    ax[1].set_ylim([zmin, zmax])

    plt.tight_layout()
    plt.savefig('figs/' + str(i).zfill(4) + '.png', bbox_inches='tight')

    # retain old array section
    old_arr[i, :] = _arr[stk, :]
