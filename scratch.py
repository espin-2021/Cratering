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
minD = 5
maxD = 300
Kx = 1.0    #Scaling coefficient (Howard, 2007)
delta = 2.0 #km, scaling exponent (Howard, 2007)

#### Plot the rough background:
###############################
mg = make_noisy_surface(grid_size, cell_size, rf=1);
# plot_topo_profile(mg, xy, Title = 'a');
# plot_grid(mg, grid_size, cell_size, Title='1.1');


## Add Craters:
###############
# ## Add craters to the background surface (new method): 
# time_interval = [4.0, 1.9]
# size_interval = [minD, maxD]
# poisson_intervals = True;     

# mg = add_craters2(mg, time_interval, size_interval, grid_size, cell_size, 
#                   poisson_intervals=True, rim = True);
# plot_grid(mg, grid_size, cell_size, Title='1.2');


# ## Add craters to the background surface (old method): 
NDs = []
for D in range(minD, maxD):
    ND = Kx * D **-delta;
    NDs.append(ND);

mg2 = add_craters1(mg, grid_size, cell_size, N_craters, NDs, minD, maxD, rim = True)
plot_grid(mg2, grid_size, cell_size, Title='1.3');
    
# ## Add a central crater:
###########################
# mg = central_crater(mg, 250, grid_size, cell_size, rim = True);
# plot_grid(mg, grid_size, cell_size, Title='1.3');

