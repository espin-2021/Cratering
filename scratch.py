# -*- coding: utf-8 -*-
"""
Created on Thu May 11 07:28:28 2023

@author: erb2734
"""

from crater_functions import *

### Set up parameters of the grid
grid_size = 1000 #grid length distance in km
cell_size = 5 #size of a cell, (km)
N_craters = 500 # number of craters to add
cc_d = 400 # diameter of the central crater ('cc'), (km)
minD = 15
maxD = 200
Kx = 1.0    #Scaling coefficient (Howard, 2007)
delta = 2.0 #km, scaling exponent (Howard, 2007)

#### Plot the rough background:
mg = make_noisy_surface(grid_size, cell_size, rf=1);
# plot_topo_profile(mg, xy, Title = 'a');
# plot_grid(mg, grid_size, cell_size, Title='1.1');

## Add craters to the background surface: 
NDs = []
for D in range(minD, maxD):
    ND = Kx * D **-delta;
    NDs.append(ND);

mg = do_cratering(mg, grid_size, cell_size, N_craters, NDs, minD, maxD, rim = True);
# plot_topo_profile(mg, xy, Title = 'b');
plot_grid(mg, grid_size, cell_size, Title='1.2');

# ## Add a central crater:
# mg = central_crater(mg, 250, grid_size, cell_size, rim = True);
# plot_grid(mg, grid_size, cell_size, Title='1.3');


### Testing the rimless crater code:
mg2 = make_noisy_surface(grid_size, cell_size, rf=1);
# plot_topo_profile(mg, xy, Title = 'a');
# plot_grid(mg2, grid_size, cell_size, Title='2.1');

## Add craters to the background surface: 
NDs = []
for D in range(minD, maxD):
    ND = Kx * D **-delta
    NDs.append(ND)
    
mg2 = do_cratering(mg2, grid_size, cell_size, N_craters, NDs, minD, maxD, rim = True);
# plot_topo_profile(mg, xy, Title = 'b');
plot_grid(mg2, grid_size, cell_size, Title='2.2');

# # ## Add a central crater:
# # mg2 = central_crater(mg2, 250, grid_size, cell_size, rim = False);
# # plot_grid(mg2, grid_size, cell_size, Title='2.3');
