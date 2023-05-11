# -*- coding: utf-8 -*-
"""
Created on Thu May 11 07:28:28 2023

@author: erb2734
"""

from crater_functions import *

### Set up parameters of the grid
grid_size = 500 #grid length distance in km
spacing = 10 #size of a cell, (km)
xy = int(grid_size/spacing) #grid length distance, in number of cells
N_craters = 100 # number of craters to add
cc_d = 400 # diameter of the central crater ('cc'), (km)
minD = spacing
maxD = int(grid_size/2)
Kx = 1.0    #Scaling coefficient (Howard, 2007)
delta = 2.0 #km, scaling exponent (Howard, 2007)

#### Plot the rough background with a slope:
mg = make_noisy_surface(xy, spacing, rf=1);
# plot_topo_profile(mg, xy, Title = 'a');
plot_grid(mg, xy, grid_size, Title='1');

NDs = []
for D in range(minD, maxD):
    ND = Kx * D **-delta
    NDs.append(ND)

## Add craters to the background surface: 
mg = do_cratering(N_craters, NDs, minD, maxD, xy, mg, spacing);
# plot_topo_profile(mg, xy, Title = 'b');
plot_grid(mg, xy, grid_size, Title='2');

## Add a central crater
mg = central_crater(mg, 250, xy, spacing, rim = True);
plot_grid(mg, xy, grid_size, Title='3');