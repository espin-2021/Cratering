### Code to randomly generatate a number of craters
######################################################
## Some Resources??
## https://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/
## beta, a = 1, b = 10, from here? https://cmdlinetips.com/2018/03/probability-distributions-in-python/#:~:text=There%20are%20at%20least%20two%20ways%20to%20draw,9%20most%20commonly%20used%20probability%20distributions%20using%20SciPy.stats.

import numpy as np
from landlab import RasterModelGrid, imshow_grid
import random as rdm
import matplotlib.pyplot as plt
import matplotlib as mpl
from crater_functions import weighted_choice_sub, crater_depth

##Set the colourmap
cmap1 =  mpl.cm.get_cmap("Greys_r").copy().reversed()
cmap2 = mpl.cm.get_cmap("Spectral").copy().reversed()

### Define some variables for our grid!
##################################
xy = 200 #Set the number of nodes in both x and y space
spacing = 1
zfactor = 10 #set factor to multiply random noise by

## Set up a LandLab model grid object, with random topographic noise on the order of "zfactor"
##################################
mg = RasterModelGrid((xy,xy), xy_spacing = spacing); #see above for variables
z = mg.add_zeros('topographic__elevation', at='node') #create an array of zeros for each node of the model grid
#### Add noise to the surface
rf = 0.1 #randomness factor (how rough initial terrain is)
np.random.seed(30) # Keep this constant (e.g., at 30) so the initial randomness it always the same
noise = (np.random.rand(mg.number_of_nodes)) #between 0 and 1
z += (noise) * rf # make the noise large enough relative to crater
mg.at_node["topographic__elevation"] *= zfactor 

# Show the initial, noisy topography:
##################################
plt.figure() #show the model grid, values of elevation
imshow_grid(mg, 'topographic__elevation', plot_name= 'Initial topographic surface', cmap=cmap2, colorbar_label='Elevation [m]')

## Set the "distribution" as of crater size as weights, and choose a value from those weights
##################################
###  These parameters describe the population frequency for crater diameters:
Kx = 1.0    #Scaling coefficient (Howard, 2007)
delta = 2.0 #km, scaling exponent (Howard, 2007)
minD = spacing * 3 # want it to be greater than spacing, three times the smallest cell
maxD = int(xy/4) #max diameter (m), a quarter of the total domain width
NDs = []
for D in range(minD, maxD):
    ND = Kx * D **-delta
    NDs.append(ND)

## Blast the surface with craters!
##################################
#Ncraters = 201 #Number of craters to add
Ncraters = int(input("How many craters do you want to add? ")) # see console to add your own number!
count = 0
rdm.seed(50) #Chose random seed number 50 (this ensures crater locations are same every time)
for i in range(Ncraters): #For N number of craters
    count += 1
    a = weighted_choice_sub(NDs); 
    diameter = list(range(minD, maxD))[a]
    cratercenter = (rdm.randint(1, xy), rdm.randint(1, xy))
    d = mg.calc_distances_of_nodes_to_point(cratercenter)
    crater_depth(d, diameter, mg, d_ref = 7)

## Make figure of final cratered surface:
##################################
hs = mg.calc_hillshade_at_node(elevs='topographic__elevation') #create hillshade file
topo = mg.field_values('node', 'topographic__elevation').reshape((xy, xy))
hill = np.reshape(hs, (xy, xy))

fig, ax = plt.subplots() #initiate figure
img1 = plt.imshow(hill, cmap=cmap1, alpha=1)
img2 = plt.imshow(topo, cmap=cmap2, alpha=0.6)
fig.colorbar(img2,ax=ax, label="Elevation [m]")
plt.title("Final cratered surface")
plt.show()

    