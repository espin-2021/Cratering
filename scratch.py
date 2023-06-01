# -*- coding: utf-8 -*-
"""
Created on Thu May 11 07:28:28 2023

@author: erb2734
"""
from crater_functions import *

import numpy as np
import copy
import random as rdm
import time
from matplotlib.pyplot import title, show, figure, plot, subplot, xlabel, ylabel, cm, clf
from landlab import RasterModelGrid, RadialModelGrid, imshow_grid, NodeStatus
from landlab.components import FlowAccumulator, FastscapeEroder, ChannelProfiler, DepressionFinderAndRouter, ErosionDeposition
from landlab.plot.imshow import imshow_grid_at_node
from landlab.values import random, plane
from landlab.io.esri_ascii import read_asc_header, read_esri_ascii, write_esri_ascii
# import warnings
# warnings.filterwarnings("ignore")

### Set up parameters of the grid
grid_size = 10000 #grid length distance in m
cell_size = 100 #size of a cell, m

slope = 0.003
time_interval = [4.0, 3.9]
size_interval = [1, 5]
poisson = True;

cc_d = 8 # diameter of the central crater ('cc'), (km)

cmap2 = copy.copy(cm.get_cmap("Spectral_r")); #define colour scheme for the topography
cmap1 = copy.copy(cm.get_cmap("Greys_r")); #define colour scheme for the hillshade 

figdir = 'C:/Users/erb2734/Documents/LANDSCAPE-EVOLUTION-MODELS/code-notebooks-etc/figs/'
griddir = 'C:/Users/erb2734/Documents/LANDSCAPE-EVOLUTION-MODELS/code-notebooks-etc/grids/'
### ASCII file saving info: https://landlab.readthedocs.io/en/master/reference/io/esri_ascii.html#landlab.io.esri_ascii.read_esri_ascii

#### (1) Plot the rough background:
###############################
mg = make_noisy_surface(grid_size, cell_size, slope = slope, rf=1);

## Plot
# imshow_grid(mg, 'topographic__elevation', cmap="Spectral_r", colorbar_label = 'Elevation (m)')  # plot the elevation field values
# plt.title('(1) Background (noisy, sloping)'); #add a title
# plt.xlabel("X [m]"); plt.ylabel("Y [m]") #add x and y axis labels
# show();

### (2) Add Craters to the background
###############################
mg = add_craters2(mg, time_interval, size_interval, int(grid_size/1000), (cell_size/1000), poisson_intervals=poisson, rim = True);

## Plot
imshow_grid(mg, 'topographic__elevation', cmap="Spectral_r", colorbar_label = 'Elevation (m)')  # plot the elevation field values
plt.title('(2) Background (w. craters)'); #add a title
xlabel("X [m]"); ylabel("Y [m]") #add x and y axis labels
show();

#### (3) Add a Central Crater (with NO RIM):
############################################
mg = central_crater(mg, cc_d, grid_size, cell_size, rim = False);

imshow_grid(mg, 'topographic__elevation', cmap="Spectral_r", colorbar_label = 'Elevation (m)')  # plot the elevation field values
plt.title('(3) Grid with a Central Crater'); #add a title
plt.xlabel("X [m]"); plt.ylabel("Y [m]") #add x and y axis labels
show();


# #### Save the model grid, for ease of access later:
# f = griddir + 'temp-grid.asc'
# write_esri_ascii(f, mg, clobber = True);   ## Save the raster model grid as an ascii file
 

# #### LANDLAB:
# #########################################

# ### Set up the Channel Profiler Tool:
# def runChannelProfiler(Threshold, Title = None):
#     Threshold = 10**Threshold;
#     pf = ChannelProfiler(mg, 
#                           channel_definition_field='drainage_area',
#                           number_of_watersheds = 15, 
#                           minimum_outlet_threshold=0,
#                           main_channel_only=False, 
#                           outlet_nodes=None, 
#                           minimum_channel_threshold=0, 
#                           cmap=cmap1); #number_of_watersheds=5, #number_of_watersheds=15
#     pf.run_one_step();
#     plt.figure();
#     pf.plot_profiles_in_map_view(colorbar_label='Elevation [m]', plot_name = Title, shrink=0.65, color_for_closed=None);
#     plt.show();
#     return pf

# #### Set up parameters of the grid & landlab models
# kt = 0.0 ; #"threshold_sp"; 
# min_channel_thresh = 10; #for channel profile plotting ('channel profiler')
# dt = 1. ; ## timestep length
# steps = 30; ## number of steps

# print ('setting up model....');
# count = 0
# for runoff_rate in [1.0]: ##5.0, 10.0
    
#     for m in [1.0]: ##0.3, 0.5, 0.7, 1.0; ##prev: 1.0 worked
        
#         for n in [0.2]: ##0.1, 0.2, 0.5, 1.0; ##prev: 0.2 worked
            
#             for K_sp in [0.001]: ##0.0001, 0.001, 0.01; ##prev: 0.001 worked
#                 count += 1
    
#                 ## Load a fresh grid
#                 f = griddir + 'temp-initial-crater.asc'
#                 (mg, data) = read_esri_ascii(f, name = 'topographic__elevation');

#                 ## Set the boundary conditions
#                 z = mg.at_node['topographic__elevation'];
#                 np.isclose(z[mg.core_nodes].sum(), 1.); #I don't know what this does?
#                 mg.set_closed_boundaries_at_grid_edges(True, False, True, True) #True, False, True, False
                
#                 ## Set up Calculators for Flow accuulation & Steam Power
#                 df = DepressionFinderAndRouter(mg, routing = "D8", reroute_flow=False); #initiate the depression finder and router
#                 fa = FlowAccumulator(mg, 'topographic__elevation', flow_director="D8", depression_finder = df, runoff_rate = runoff_rate); ## Flow Accumulation  
#                 fa.run_one_step();
#                 sp = FastscapeEroder(mg, K_sp=K_sp, m_sp=m, n_sp=n, threshold_sp=kt, discharge_field="surface_water__discharge") ## Stream Power; #erode_flooded_nodes = False
#                 sp.run_one_step(dt=dt);
                
#                 ## Run the model over a few steps
#                 print ('running model....');
                
#                 steps = steps;
#                 for i in range(steps+1): 
#                     fa.run_one_step()
#                     sp.run_one_step(dt=dt)
                    
#                     if (i in np.arange(10,steps+1,10)):
#                         ###PLOT 1
#                         figure() #show the model grid, values of elevation
#                         imshow_grid(mg, 'topographic__elevation',
#                                     plot_name= 'Topography after {} Model Runs'.format(i), 
#                                     symmetric_cbar=False, cmap=cmap2,
#                                     allow_colorbar=True, colorbar_label='Elevation [m]',
#                                     color_for_closed='black')  # plot the elevation field values
                        
#                         f = figdir + 'temp-Topo_Eroded_ro{}_m{}_n{}_K{}_{}of{}.jpeg'.format(runoff_rate, m, n, K_sp, i, steps)
#                         plt.savefig(f);
#                         show();
                
#                         ###PLOT 2
#                         figure() #show the model grid, values of ???
#                         imshow_grid(mg, 'surface_water__discharge', cmap="bone")  # plot the elevation field values
#                         show();
                
#                         ###PLOT 4       
#                         try: 
#                             profiler = runChannelProfiler(np.log10(min_channel_thresh))
#                         except:
#                             print('No Channel Profiler')
                            
#                         ####SAVE THE GRID AS ASCII FILE
#                         f = griddir + 'temp-Eroded_ro{}_m{}_n{}_K{}_{}of{}.asc'.format(runoff_rate, m, n, K_sp, i, steps)
#                         write_esri_ascii(f, mg, clobber = True); 
                
                            
#                 y = mg.field_values('node', 'topographic__elevation').reshape((xy, xy))[int(xy/2)]; #reshape the topography so it's easy to find the midpoint row (index = xy/2)
#                 x = (np.arange(xy) + 1); #create the x values (from 1 to xy)
#                 plt.plot(x, y, 'k-'); #create the plot, with a black line for topography
#                 plt.title('Final Topo Profile'); #add a title
#                 plt.show() #display the figure
                
#                 print( (time.process_time() - start)/60, 'min : Done model run {} of 1'.format(count) );
#                 ## number of model runs = #runoff options x #m options x #n options x #K options








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

