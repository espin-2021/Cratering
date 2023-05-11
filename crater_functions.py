"""Cratering Functions."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib
import mpl_toolkits.axes_grid1 as axtk
import os
import shutil
import sys
import yaml
import warnings
from packaging import version
import craterstats as cst
from landlab import RasterModelGrid
from landlab import imshow_grid


rn_gen = np.random.default_rng(seed=3);

def weighted_choice_sub(weights):
    ''' randomly generate a number and see which weight number in the input list it falls under,
    return the index of that weight '''
    
    rnd = rn_gen.random(1) * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i

def make_noisy_surface(grid_size, cell_size, rf=1):
    ''' Generate a surface with random topography.
     Parameters
    ----------
    grid size : integer
        Size of the domain, in km 
        (note: domain is square so length = width)
    
    cell_size : integer
        Size of each cell in km.
    
    rf : int, float
        The multiplier (in km) to add to increase or decrease randomness by factor rf. 
        (randomness factor, default = 1)
        
    Returns
    ----------
    mg : landlab.grid.raster.RasterModelGrid
        Landlab raster model grid of the landscape
    '''
    rn_gen = np.random.default_rng(seed=1);
    
    xy = int(grid_size / cell_size) ## The number of cells along each axis of the domain (i.e. length in # cells)
    
    mg = RasterModelGrid((xy,xy), xy_spacing = cell_size); #initiate surface; see above for variables
    z = mg.add_zeros('topographic__elevation', at='node') #create an array of zeros for each node of the model grid  
    
    z += rn_gen.random(mg.number_of_nodes) # Add random elevation values at each node 
    
    mg.at_node["topographic__elevation"] *= rf  # make the noise large enough relative to crater
    return mg


def crater_depth(mg, d, diameter, d_ref=7, rim = True):
    """
    Define in and out of crater changes to topography.

    Parameters
    ----------
    mg : landlab.grid.raster.RasterModelGrid
        Landlab raster model grid of the landscape
   
    d : np.ndarray
        Array of distances of nodes to the center of the crater.
        From `landlab.grid.raster.RasterModelGrid.calc_distances_of_nodes_to_point`
    
    diameter : int, float
        Diameter of the crater

    d_ref: int, float
        Diameter at which craters transition from "simple" to "complex" structural type
        Default = 7 km (relevant to Mars)
        
    rim : boolean, default = True
        whether the crater generated has a crater rim or not. 

    Returns
    -------
    mg : landlab.grid.raster.RasterModelGrid
        Landlab raster model grid after crater has modified the topography
    """
    radius = diameter/2
    diameter *= 1000
    d_ref *= 1000

    if diameter <= d_ref:
        H1 = 2.54*diameter**0.67
        H2 = 1.93*diameter**0.52
        m = 0.73*diameter**0.11  # value: 2 to 3

    elif diameter > d_ref:
        H1 = 12.20*diameter**0.49
        H2 = 0.79*diameter**0.6
        m = 0.64*diameter**0.13  # value: 2 to 3

    # Howard et al, 2007:
    # "The exponent n is constrained such that volume deposited on the rim
    # equals the volume excavated from the bowl and ranges from a value of
    # about 3 for a 7 km crater to 3.5 for a 250 km crater.""
    n = 3
    H2H1 = H2 - H1

    if rim == True: ## (DEFAULT) If the user wants the craters to have a rim
        incrater = d[d <= radius]
        in_idx = np.where(d <= radius)[0]
        # equation for inside the crater
        inDepth = (H2H1 + H1*((2*(incrater*1000))/(diameter))**m)/1000 #In KM (hence /1000)
        mg.at_node['topographic__elevation'][in_idx] += inDepth
    
        outcrater = d[d > radius]
        out_idx = np.where(d > radius)[0]
        # equation for outside the crater (ejecta!)
        outDepth = (H2*((2*(outcrater*1000))/(diameter))**-n)/1000 # IN KM (hence /km)
        mg.at_node['topographic__elevation'][out_idx] += outDepth
    
    elif rim == False: ## If the user doesn't want the crater to have any rims
        incrater = d[d <= radius * 0.9] #Only excavate the crater for 90% of the crater radius
        in_idx = np.where(d <= radius * 0.9)[0] #Only excavate the crater for the first 90% of the crater radius
        ##90% of the crater radius ensures there's no rim on the crater for craters up to about 500 km in diameter
        ## For much smaller craters (< 250 km), a smaller value than 90% could be used, but it still makes a reasonable crater
        
        # equation for inside the crater
        inDepth = (H2H1 + H1*((2*(incrater*1000))/(diameter))**m)/1000 #In KM (hence /1000)
        mg.at_node['topographic__elevation'][in_idx] += inDepth
    
     
    return mg

def do_cratering(mg, grid_size, cell_size, Ncraters, NDs, minD, maxD, rim = True):
    """
    Add craters to some landlab raster model.

    Parameters
    ----------
    mg : landlab.grid.raster.RasterModelGrid
        Landlab raster model grid of the landscape
        
    grid size : integer
        Size of the domain, in km 
        (note: domain is square so length = width)
    
    cell_size : integer
        Size of each cell in km.
        
    Ncraters : int
        Number of craters that impact

    NDs : list
        List of weights for random sampling

    minD : int
        Minimum crater diameter, km

    maxD : int
        Maximum crater diameter, km
        
    rim : boolean, default = True
        whether the crater generated has a crater rim or not. 

    Returns
    -------
    mg : landlab.grid.raster.RasterModelGrid
        Landlab raster model grid after craters have modified the topography

    """
    rn_gen = np.random.default_rng(seed=1); ## create an instance of the Generator class
    
    xy = int(grid_size / cell_size);

    for i in range(Ncraters):  # For N number of craters
        a = weighted_choice_sub(NDs)
        diameter = list(range(minD, maxD))[a]
        cratercenter = (rn_gen.integers(1, grid_size, endpoint=True), rn_gen.integers(1, grid_size, endpoint = True))
        d = mg.calc_distances_of_nodes_to_point(cratercenter)

        crater_depth(mg, d, diameter, d_ref=7, rim = rim)

    return mg

def central_crater(mg, diameter, grid_size, cell_size, rim = True):
	"""
	Add a central crater to a Landlab raster model grid
	
	Parameters
	----------
	mg : Landlab.grid.raster.RasterModelGrid
	Landlab raster model grid of the landscape
	
	diameter: int
	Diameter of crater to be added at the central node
	
    grid size : integer
        Size of the domain, in km 
        (note: domain is square so length = width)
    
    cell_size : integer
        Size of each cell in km.
    
    rim : boolean, default = True
	whether or not the central crater has a rim (True) or no rim (False). 
    This argument is passed to the function crater_depth
    
	Returns
	-------
	mg : landlab.grid.raster.RasterModelGrid
	Landlab raster model grid after a central crater has modified the topography

	"""
	xy = int( grid_size / cell_size);

	d = mg.calc_distances_of_nodes_to_point( (int((grid_size)/2), int((grid_size)/2)) );

	crater_depth(mg, d, diameter, d_ref=7, rim = rim);
	
	return(mg)

def plot_topo_profile(mg, grid_size, cell_size, Title = 'Title'):
    """
    Parameters
	----------
	mg : Landlab.grid.raster.RasterModelGrid
	Landlab raster model grid of the landscape
    
    grid size : integer
        Size of the domain, in km 
        (note: domain is square so length = width)
    
    cell_size : integer
        Size of each cell in km.
        
    Title : string
    title to add to the plot
    
    Returns
    ---------
    Plot of profile from left to right across the whole model grid width
    (through the center of the grid).
    
    """
    xy = int(grid_size / cell_size);
    y = mg.field_values('node', 'topographic__elevation').reshape((xy, xy))[int(xy/2)]; #reshape the topography so it's easy to find the midpoint row (index = xy/2)
    x = (np.arange(xy) + 1); #create the x values (from 1 to xy)
    plt.plot(x, y, 'k-'); #create the plot, with a black line for topography
    plt.title(Title); #add a title
    plt.show() #display the figure
    

def plot_grid(mg, grid_size, cell_size, Title='Title'):
    """
    Parameters
    ----------
    mg : Landlab.grid.raster.RasterModelGrid
        Landlab raster model grid of the landscape
     
           grid size : integer
        Size of the domain, in km 
        (note: domain is square so length = width)
    
    cell_size : integer
        Size of each cell in km.
    
    Title: string
        the title for the plot

    Returns
    -------
    None.

    """
    xy = int(grid_size / cell_size);
    
    cmap2 = mpl.cm.get_cmap("Spectral_r"); #define colour scheme for the topography
    cmap1 = mpl.cm.get_cmap("Greys_r"); #define colour scheme for the hillshade
    
    hs = mg.calc_hillshade_at_node(elevs='topographic__elevation') #create hillshade file
    topo = mg.field_values('node', 'topographic__elevation').reshape((xy, xy)) #reshape the topography to the right (square) shape for display
    hill = np.reshape(hs, (xy, xy)); #reshape the hillshade to the right (square) shape for display

    fig, ax = plt.subplots() #initiate figure
    img1 = plt.imshow(hill, cmap=cmap1, alpha=1, extent = [0,grid_size, 0, grid_size]) #plot hillshade
    img2 = plt.imshow(topo, cmap=cmap2, alpha=0.6, extent = [0,grid_size, 0, grid_size]) #plot topograpy
    fig.colorbar(img2,ax=ax, label="Elevation [km]") #add & label the colorbar
    plt.title(Title); #add a title
    plt.xlabel("X [km]"); plt.ylabel("Y [km]") #add x and y axis labels
    plt.show() #show the figure