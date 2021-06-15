### Code to randomly generatate a number of craters
######################################################
## Some Resources??
## https://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/
## beta, a = 1, b = 10, from here? https://cmdlinetips.com/2018/03/probability-distributions-in-python/#:~:text=There%20are%20at%20least%20two%20ways%20to%20draw,9%20most%20commonly%20used%20probability%20distributions%20using%20SciPy.stats.

import numpy as np
from landlab import RasterModelGrid, imshow_grid, NodeStatus, values
import random as rdm
import matplotlib as mpl

##Set the colourmap
#cmap = "Greys_r"; cmap2 = "bone";
cmap =  mpl.cm.get_cmap("Greys_r").copy()
cmap2 = mpl.cm.get_cmap("bone").copy()

### Define some Variables!
xy = 200 #Set the number of nodes in both x and y space
spacing = 1
zfactor = 10 #set factor to multiply random noise by
 
## Set up a LandLab model grid object, with random topographic noise on the order of "zfactor"  
mg = RasterModelGrid((xy,xy), xy_spacing = spacing); #see above for variables
z = mg.add_zeros('topographic__elevation', at='node') #create an array of zeros for each node of the model grid
#### Add noise to the surface
rf = 0.1 #randomness factor (how rough initial terrain is)
np.random.seed(30) # Keep this constant (e.g., at 30) so the initial randomness it always the same
noise = (np.random.rand(mg.number_of_nodes)) #between 0 and 1
z += (noise) * rf #if you need a multiplier
#### Increase the zfactor to 10 to get realistic values of noisy elevation relative to the crater
mg.at_node["topographic__elevation"] *= zfactor 

mpl.pyplot.figure() #show the model grid, values of elevation
imshow_grid(mg, 'topographic__elevation', 
            plot_name= 'Initial Topography', 
            symmetric_cbar=False, cmap=cmap,
            allow_colorbar=True, colorbar_label='Elevation [m]',
            shrink=1.,
            color_for_closed='black',
            color_for_background=None,            
            show_elements=False)  


## The quickest weighted choice method:
def weighted_choice_sub(weights):
    ''' randomly generate a number and see which weight number in the input list it falls under,
    return the index of that weight '''
    rdm.seed(25)
    rnd = rdm.random() * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i          
            
## Assigning crater shape 
def crater_depth(d, diameter, mg, d_ref = 7):
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
    radius = diameter / 2;
    
    if diameter <= d_ref:
        H1 = 2.54*diameter**0.67;
        H2 = 1.93*diameter**0.52;
        m = 0.73*diameter**0.11; ## value: 2 to 3
    
    elif diameter > d_ref:
        H1 = 12.20*diameter**0.49;
        H2 = 0.79*diameter**0.6;
        m = 0.64*diameter**0.13; ## value: 2 to 3
        
    n = 3 ## Howard et al, 2007: "The exponent n is constrained such that volume deposited on the rim equals the volume excavated from the bowl and ranges from a value of about 3 for a 7 km crater to 3.5 for a 250 km crater.""
    H2H1 = H2 - H1
    
    incrater = d[d<=radius]
    in_idx = np.where(d <= radius)[0]
    inDepth = H2H1 + H1*((2*incrater)/(diameter))**m #equation for inside the crater
    mg.at_node['topographic__elevation'][in_idx] += inDepth
    
    outcrater = d[d>radius]
    out_idx = np.where(d > radius)[0]
    outDepth = H2*((2*outcrater)/(diameter))**-n #equation for outside the crater (ejecta!)
    mg.at_node['topographic__elevation'][out_idx] += outDepth
    return mg        
    

## Set the "distribution" as weights, and choose a value from those weights
##################################
###  ND = Kx * D ** -delta
Kx = 1.0    #Scaling coefficient
delta = 2.0 #km, scaling exponent
Ncraters = 1001 #Number of craters to add
minD = 1 #int(spacing*3) # min diameter (m), xy_spacing times 3 i.e. more tha 3 cells in diameter
maxD = 100 #int((mg.number_of_nodes)/4) #max diameter (m), a quarter of the total domain width
NDs = []
for D in range(minD, maxD):
    ND = Kx * D **-delta
    NDs.append(ND)
    
    
count = 0
rdm.seed(50) #Chose random seed number 50 (this ensures crater locations are same every time)
for i in range(Ncraters): #For N number of craters
    count += 1
    a = weighted_choice_sub(NDs); 
    diameter = list(range(minD, maxD))[a]
    #print(diameter)
    cratercenter = (rdm.randint(1, xy), rdm.randint(1, xy))
    d = mg.calc_distances_of_nodes_to_point(cratercenter)
    
    crater_depth(d, diameter, mg, d_ref = 7)

    if count%500 == 0:
        hs = mg.calc_hillshade_at_node(elevs='topographic__elevation') #create hillshade file
        mpl.pyplot.figure() #show the model grid, values of elevation
        imshow_grid(mg,hs,cmap=cmap,allow_colorbar=False) # Easier to vizualize for now?
        #imshow_grid(mg, 'topographic__elevation', # this is for visualizing topography (no hillshade; maybe can add on top as transparent layer?)
                    #plot_name= 'Cratered Topography', 
                    #symmetric_cbar=False, cmap=cmap,
                    #allow_colorbar=True, colorbar_label='Elevation [m]',
                    #shrink=1.,
                    #color_for_closed='black',
                    #color_for_background=None,            
                    #show_elements=False)  
