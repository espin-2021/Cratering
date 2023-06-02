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
figdir = 'C:/Users/erb2734/Documents/LANDSCAPE-EVOLUTION-MODELS/code-notebooks-etc/initial-grids/'
griddir = 'C:/Users/erb2734/Documents/LANDSCAPE-EVOLUTION-MODELS/code-notebooks-etc/initial-grids/'
### ASCII file saving info: https://landlab.readthedocs.io/en/master/reference/io/esri_ascii.html#landlab.io.esri_ascii.read_esri_ascii
 
cmap2 = copy.copy(cm.get_cmap("Spectral_r")); #define colour scheme for the topography
cmap1 = copy.copy(cm.get_cmap("Greys_r")); #define colour scheme for the hillshade 
   
### Set up parameters of the grid, See loops below for iterating over different things.  
grid_size = 50000; #m, grid length or width (square)
slope = 0.003; #rise/run (dimensionless)
rf = 100; #m, factor to multiply noise by  
time_interval = [4.2, 3.8]; #in Ga
poisson = True;
cc_d = 10; # diameter of the central crater ('cc'), (km)

################################
####### MAKE THE SURFACE #######
################################
### There are two options:
### (A) With craters on the background, from CSFD
### (B) Without background craters (i.e. just noise + central crater)

## Notes:
## It would actually be quicker to do this without the functions, 
## e.g. generate an array of locations plus the CSFD for [1/0.5, 2/5] and use it for all the simulations

### Option (A) Make Grids with craters on the background
for cell_size in [100, 50, 20, 10]: #meters,  ## 100, 50, 20
    cell_size = cell_size; #m, cell size
    xy = int(grid_size / cell_size); ## number of nodes
    
    for size_interval in [[1, 2], [1, 5], [0.5, 2], [0.5, 5]]: 
        
        for rims in [True, False]:
                
            ################################
            ####### MAKE THE SURFACE #######
            ################################
            
            #### (1) Make a rough background:
            ###############################
            print ('\n making surface....');
            # print('grid:{} cells:{} slope:{}'.format(grid_size, cell_size, slope));
            mg = make_noisy_surface(grid_size, cell_size, slope = slope, rf = rf);
            print((time.process_time() - start)/60, 'min : model grid with + noise, slope');
            
            ## (2) Add Craters:
            #############
            print ('adding craters to background....');
            mg = add_craters2(mg, time_interval, size_interval, poisson_intervals = poisson, rim = rims);
            
            imshow_grid(mg, 'topographic__elevation', cmap="Spectral_r", colorbar_label = 'Elevation (m)')  # plot the elevation field values
            plt.title('Background (w. craters)'); #add a title
            plt.xlabel("X [m]"); plt.ylabel("Y [m]") #add x and y axis labels
            f1 = 'initial_{}km_{}m_S{}_rf{}_t{}-{}_{}to{}km_rims{}_nocenter'.format(grid_size/1000, cell_size, slope, rf, time_interval[0], time_interval[1], size_interval[0], size_interval[1], rims)
            plt.savefig(figdir+'plot-'+f1+'.jpeg');
            plt.show();
            plt.clf();
            
            ### (2b) Save the grid! (as an ascii file, using landlab save as ascii function)
            write_esri_ascii(figdir+'grid-'+f1+'.asc', mg, clobber = True); 
            
            print((time.process_time() - start)/60, 'min : model grid with noise, slope, + craters');
            
            ### (3) Add a central crater: 
            #########################
            print ('adding central crater ({} km)....'.format(cc_d));
            mg = central_crater(mg, cc_d, rim = False);
            
            imshow_grid(mg, 'topographic__elevation', cmap="Spectral_r", colorbar_label = 'Elevation (m)')  # plot the elevation field values
            plt.title('Topo w. central crater \n {}km_{}m_S{}_rf{}_t{}-{}_{}to{}km_rims{}_cc{}km'.format(grid_size/1000, cell_size, slope, rf, time_interval[0], time_interval[1], size_interval[0], size_interval[1], rims, cc_d)); #add a title
            plt.xlabel("X [m]"); plt.ylabel("Y [m]") #add x and y axis labels
            f = 'initial_{}km_{}m_S{}_rf{}_t{}-{}_{}to{}km_rims{}_cc{}km'.format(grid_size/1000, cell_size, slope, rf, time_interval[0], time_interval[1], size_interval[0], size_interval[1], rims, cc_d)
            plt.savefig(figdir+'plot-'+f+'.jpeg');
            plt.show();
            plt.clf();
            
            
            ### (4) Save the grid! (as an ascii file, using landlab save as ascii function)
            write_esri_ascii(figdir+'grid-'+f+'.asc', mg, clobber = True);   ## Save the raster model grid as an ascii file
    
            print( (time.process_time() - start)/60, 'min : model grid with noise, slope, + central crater (rimless) [saved]');

start2 = time.process_time(); #Start time, just for timing purposes.
##### Option (B) Make grids without craters on the background
for cell_size in [100, 50, 20, 10]: #meters, 
    cell_size = cell_size; #m, cell size
    xy = int(grid_size / cell_size); ## number of nodes

    ################################
    ####### MAKE THE SURFACE #######
    ################################
    
    #### (1) Make a rough background:
    ###############################
    print ('\n (1) making surface....');
    # print('grid:{} cells:{} slope:{}'.format(grid_size, cell_size, slope));
    mg = make_noisy_surface(grid_size, cell_size, slope = slope, rf = rf);
    print((time.process_time() - start2)/60, 'min : model grid with + noise, slope');

    ### (3) Add a central crater: 
    #########################
    print (' (2) adding central crater ({} m)....'.format(cc_d));
    mg = central_crater(mg, cc_d, rim = False);
    
    imshow_grid(mg, 'topographic__elevation', cmap="Spectral_r", colorbar_label = 'Elevation (m)')  # plot the elevation field values
    plt.title('Topo w. central crater \n {}km_{}m_S{}_rf{}_no-bkgrd-cs_cc{}km.asc'.format(grid_size/1000, cell_size, slope, rf, cc_d)); #add a title
    plt.xlabel("X [m]"); plt.ylabel("Y [m]") #add x and y axis labels
    f = 'initial_{}km_{}m_S{}_rf{}_no-bkgrd-cs_cc{}km'.format(grid_size/1000, cell_size, slope, rf, cc_d)
    plt.savefig(figdir+'plot-'+f+'.jpeg')
    show();
    plt.clf();
    
    ### (4) Save the grid! (as an ascii file, using landlab save as ascii function)
    write_esri_ascii(griddir+'grid-'+f+'.asc', mg, clobber = True);   ## Save the raster model grid as an ascii file
    print( (time.process_time() - start2)/60, 'min : model grid with noise, slope, + central crater (rimless) [saved]');

print( (time.process_time() - start)/60, 'min : total time for all grid generation');
