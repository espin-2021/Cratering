# -*- coding: utf-8 -*-
"""
Created on Thu May 11 07:28:28 2023

@author: erb2734
"""

from crater_functions import *

### Set up parameters of the grid
grid_size = 1000 #grid length distance in km
cell_size = 5 #size of a cell, (km)
slope = 0.001
N_craters = 500 # number of craters to add
cc_d = 400 # diameter of the central crater ('cc'), (km)
minD = 5
maxD = 300
time_interval = [4.0, 1.9]
size_interval = [minD, maxD]
poisson_intervals = True;     

#### Plot the rough background:
###############################
mg = make_noisy_surface(grid_size, cell_size, slope = slope, rf=1);
# plot_topo_profile(mg, xy, Title = 'a');
plot_grid(mg, grid_size, cell_size, Title='1.1');

#### Add Craters
mg = add_craters2(mg, time_interval, size_interval, grid_size, cell_size, 
                  poisson_intervals=True, rim = True);
plot_grid(mg, grid_size, cell_size, Title='background');

#### Add a Central Crater with NO RIM:
###########################
mg = central_crater(mg, 500, grid_size, cell_size, rim = False);
plot_grid(mg, grid_size, cell_size, Title='w. central crater');

#### Up to this point is full working + clean code
###################################################

#### Do some LandLab!: 
import numpy as np
import copy
from landlab.components import FlowAccumulator, FastscapeEroder, ChannelProfiler
from landlab import RasterModelGrid, RadialModelGrid, imshow_grid, NodeStatus
import matplotlib as mpl
from matplotlib.pyplot import title, show, figure, plot, subplot, xlabel, ylabel
from landlab.values import random, plane
import random as rdm
import time   

### Define some Variables!
start = time.process_time() #Start time, just for timing purposes. 
xy = 200 #Set the number of nodes in both x and y space
spacing = 100
diameter = 7000 #Set the central crater diameter (m)
radius = diameter/2
slope = 0.10 #Set the slope to be imposed on the model grid domain
zfactor = 200 #set factor to multiply random noise by
rimheight = 100

##Set the colourmap
cmap = "Greys_r"; cmap2 = "bone"; 






#### OLD SCRATCH CODE DOWN HERE
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
# NDs = []
# for D in range(minD, maxD):
#     ND = Kx * D **-delta;
#     NDs.append(ND);

# mg2 = add_craters1(mg, grid_size, cell_size, N_craters, minD, maxD, rim = True)
# plot_grid(mg2, grid_size, cell_size, Title='1.3');
    
# ## Add a central crater:
###########################
# mg = central_crater(mg, 250, grid_size, cell_size, rim = True);
# plot_grid(mg, grid_size, cell_size, Title='1.3');

