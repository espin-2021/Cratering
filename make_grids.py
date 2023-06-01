# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 07:04:30 2023

@author: erb2734
"""
from crater_functions import *
import time
import copy
import numpy as np
from matplotlib.pyplot import title, figure, show, cm, clf
from landlab.io.esri_ascii import read_asc_header, read_esri_ascii, write_esri_ascii
from landlab import RasterModelGrid, imshow_grid
from landlab.plot.imshow import imshow_grid_at_node
from landlab.components import ChannelProfiler, FastscapeEroder, FlowAccumulator, DepressionFinderAndRouter, ErosionDeposition
# import warnings
# warnings.filterwarnings("ignore")

start = time.process_time(); #Start time, just for timing purposes.

## Set file paths for saving
# figdir = 'C:/Users/erb2734/Documents/LANDSCAPE-EVOLUTION-MODELS/code-notebooks-etc/'
griddir = 'C:/Users/erb2734/Documents/LANDSCAPE-EVOLUTION-MODELS/code-notebooks-etc/'
### ASCII file saving info: https://landlab.readthedocs.io/en/master/reference/io/esri_ascii.html#landlab.io.esri_ascii.read_esri_ascii
 
cmap2 = copy.copy(cm.get_cmap("Spectral_r")); #define colour scheme for the topography
cmap1 = copy.copy(cm.get_cmap("Greys_r")); #define colour scheme for the hillshade 
   
### Set up parameters of the grid, See loops below for iterating over different things.  
grid_size = 20000; #m, grid length or width (square)
slope = 0.003; #rise/run (dimensionless)
rf = 100; #m, factor to multiply noise by  
time_interval = [4.2, 4.0]; #in Ga
poisson = True;
cc_d = 10; # diameter of the central crater ('cc'), (km)

################################
####### MAKE THE SURFACE #######
################################
### There are two options:
### (A) With craters on the background, from CSFD
### (B) Without background craters (i.e. just noise + central crater)

### Option (A) Make Grids with craters on the background
for cell_size in [100]: #meters,  ## 100, 50, 20
    cell_size = cell_size; #m, cell size
    xy = int(grid_size / cell_size); ## number of nodes
    
    for size_interval in [[1, 3]]: ##[1, 2], [1, 5], [0.1, 2], [0.1, 5]
        
        ################################
        ####### MAKE THE SURFACE #######
        ################################
        
        #### (1) Make a rough background:
        ###############################
        print ('making surface....');
        # print('grid:{} cells:{} slope:{}'.format(grid_size, cell_size, slope));
        mg = make_noisy_surface(grid_size, cell_size, slope = slope, rf = rf);
        print((time.process_time() - start)/60, 'min : model grid with + noise, slope');
        
        ## (2) Add Craters:
        #############
        print ('adding craters to background....');
        mg = add_craters2(mg, time_interval, size_interval, poisson_intervals = poisson, rim = True);
        
        imshow_grid(mg, 'topographic__elevation', cmap="Spectral_r", colorbar_label = 'Elevation (m)')  # plot the elevation field values
        plt.title('Background (w. craters)'); #add a title
        plt.xlabel("X [m]"); plt.ylabel("Y [m]") #add x and y axis labels
        show();
        
        print((time.process_time() - start)/60, 'min : model grid with noise, slope, + craters');
        
        ### (3) Add a central crater: 
        #########################
        print ('adding central crater ({} m)....'.format(cc_d));
        mg = central_crater(mg, cc_d, rim = False);
        
        imshow_grid(mg, 'topographic__elevation', cmap="Spectral_r", colorbar_label = 'Elevation (m)')  # plot the elevation field values
        plt.title('Topo w. central crater'); #add a title
        plt.xlabel("X [m]"); plt.ylabel("Y [m]") #add x and y axis labels
        show();
        
        print( (time.process_time() - start)/60, 'min : model grid with noise, slope, + central crater (rimless)');
        
        # ### (4) Save the grid! (as an ascii file, using landlab save as ascii function)
        # f = griddir + 'Initial_{}km_{}m_S{}_rf{}_t{}-{}_{}to{}km_cc{}km.asc'.format(grid_size/1000, cell_size, slope, rf, time_interval[0], time_interval[1], size_interval[0], size_interval[1], cc_d)
        # write_esri_ascii(f, mg, clobber = True);   ## Save the raster model grid as an ascii file



# ##### Option (B) Make grids without craters on the background
# for cell_size in [100, 50]: #meters, 
#     cell_size = cell_size; #m, cell size
#     xy = int(grid_size / cell_size); ## number of nodes

#     ################################
#     ####### MAKE THE SURFACE #######
#     ################################
    
#     #### (1) Make a rough background:
#     ###############################
#     print ('making surface....');
#     # print('grid:{} cells:{} slope:{}'.format(grid_size, cell_size, slope));
#     mg = make_noisy_surface(grid_size, cell_size, slope = slope, rf = rf);
#     print((time.process_time() - start)/60, 'min : model grid with + noise, slope');

#     ### (3) Add a central crater: 
#     #########################
#     print ('adding central crater ({} m)....'.format(cc_d));
#     mg = central_crater(mg, cc_d, rim = False);
    
#     imshow_grid(mg, 'topographic__elevation', cmap="Spectral_r", colorbar_label = 'Elevation (m)')  # plot the elevation field values
#     plt.title('Topo w. central crater'); #add a title
#     plt.xlabel("X [m]"); plt.ylabel("Y [m]") #add x and y axis labels
#     show();
    
#     print( (time.process_time() - start)/60, 'min : model grid with noise, slope, + central crater (rimless)');
    
#     ### (4) Save the grid! (as an ascii file, using landlab save as ascii function)
#     f = griddir + 'Initial_{}km_{}m_S{}_rf{}_no-bkgrd-cs_cc{}km.asc'.format(grid_size/1000, cell_size, slope, rf, cc_d)
#     write_esri_ascii(f, mg, clobber = True);   ## Save the raster model grid as an ascii file
