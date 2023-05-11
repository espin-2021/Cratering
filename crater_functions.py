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

np.random.seed(1); ## set the numpy random seed to be consistent
rn_gen = np.random.default_rng() ## create an instance of the Generator class

def weighted_choice_sub(weights):
    ''' randomly generate a number and see which weight number in the input list it falls under,
    return the index of that weight '''
    rnd = rn_gen.random() * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i

def make_noisy_surface(xy,spacing,rf=1):
    ''' Generate a surface with random topography.
     Parameters
    ----------
    length : integer
        Number of cells for each axis (i.e., 200 for a 200x200 grid)
    spacing : integer
        Size of each cell in km.
    rf : int, float
        The multiplier (in km) to add to increase or decrease randomness by factor rf. 
        (randomness factor, default = 1)
    '''
    mg = RasterModelGrid((xy,xy), xy_spacing = spacing); #initiate surface; see above for variables
    z = mg.add_zeros('topographic__elevation', at='node') #create an array of zeros for each node of the model grid  
    
    z += rn_gen.random(mg.number_of_nodes) # Add random elevation values at each node 
    
    mg.at_node["topographic__elevation"] *= rf  # make the noise large enough relative to crater
    return mg

def crater_depth(d, diameter, mg, d_ref=7):
    """
    Define in and out of crater changes to topography.

    Parameters
    ----------
    d : np.ndarray
        Array of distances of nodes to the center of the crater.
        From `landlab.grid.raster.RasterModelGrid.calc_distances_of_nodes_to_point`
    diameter : int, float
        Diameter of the crater
    mg : landlab.grid.raster.RasterModelGrid
        Landlab raster model grid of the landscape

    d_ref: int, float
        Diameter at which craters transition from "simple" to "complex" structural type
        Default = 7 km (relevant to Mars)

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
    return mg


def do_cratering(Ncraters, NDs, minD, maxD, xy, mg, spacing):
    """
    Add craters to some landlab raster model.

    Parameters
    ----------
    Ncraters : int
        Number of craters that impact

    NDs : list
        List of weights for random sampling

    minD : int
        Minimum crater diameter

    maxD : int
        Maximum crater diameter

    xy : int
        Domain size (# cells) in one direction (domain is square)

    mg : landlab.grid.raster.RasterModelGrid
        Landlab raster model grid of the landscape

    Returns
    -------
    mg : landlab.grid.raster.RasterModelGrid
        Landlab raster model grid after craters have modified the topography

    """
    grid_size = xy * spacing
    for i in range(Ncraters):  # For N number of craters
        a = weighted_choice_sub(NDs)
        diameter = list(range(minD, maxD))[a]
        cratercenter = (rn_gen.integers(1, grid_size, endpoint=True), rn_gen.integers(1, grid_size, endpoint = True))
        d = mg.calc_distances_of_nodes_to_point(cratercenter)

        crater_depth(d, diameter, mg, d_ref=7)

    return mg

def central_crater(mg, diameter, xy, spacing):
	"""
	Add a central crater to a Landlab raster model grid
	
	Parameters
	----------
	mg : Landlab.grid.raster.RasterModelGrid
	Landlab raster model grid of the landscape
	
	diameter: int
	Diameter of crater to be added at the central node
	
	xy : int
	domain size (# cells) in one direction (domain is square)
	
	spacing: int
	size of spacing between nodes
	
	Returns
	-------
	mg : landlab.grid.raster.RasterModelGrid
	Landlab raster model grid after a central crater has modified the topography

	"""
	
	## centerpoint = ( (int((xy*spacing)/2), int((xy*spacing)/2)) )
	## print(centerpoint)
	d = mg.calc_distances_of_nodes_to_point( (int((xy*spacing)/2), int((xy*spacing)/2)) );

	crater_depth(d, diameter, mg, d_ref=7);
	
	return(mg)

def plot_topo_profile(mg, xy, Title = 'Title'):
    """
    Parameters
	----------
	mg : Landlab.grid.raster.RasterModelGrid
	Landlab raster model grid of the landscape
    
    xy: int
    number of nodes along one length of the model grid
    
    Title : string
    title to add to the plot
    
    Returns
    ---------
    Plot of profile from left to right across the whole model grid width
    (through the center of the grid).
    
    """
    y = mg.field_values('node', 'topographic__elevation').reshape((xy, xy))[int(xy/2)]; #reshape the topography so it's easy to find the midpoint row (index = xy/2)
    x = (np.arange(xy) + 1); #create the x values (from 1 to xy)
    plt.plot(x, y, 'k-'); #create the plot, with a black line for topography
    plt.title(Title); #add a title
    plt.show() #display the figure
    

def plot_grid(mg, xy, grid_size, Title='Title'):
    """
    Parameters
    ----------
    mg : Landlab.grid.raster.RasterModelGrid
        Landlab raster model grid of the landscape
        
    xy : int
        number of nodes along one length of the model grid
        
    grid_size : int
        length of one side of the model grid
    
    Title: string
        the Title for the plot

    Returns
    -------
    None.

    """
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