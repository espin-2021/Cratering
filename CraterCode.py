### Code to randomly generatate a number of craters on a planetary surface
## by Emily Bamber & Gaia Stucky de Quay (June, 2021)
######################################################
## Some Resources??
## https://eli.thegreenplace.net/2010/01/22/weighted-random-generation-in-python/
## beta, a = 1, b = 10, from here? https://cmdlinetips.com/2018/03/probability-distributions-in-python/#:~:text=There%20are%20at%20least%20two%20ways%20to%20draw,9%20most%20commonly%20used%20probability%20distributions%20using%20SciPy.stats.

import numpy as np
from landlab import imshow_grid
import random as rdm
import matplotlib.pyplot as plt
import matplotlib as mpl
from crater_functions import do_cratering, make_noisy_surface, central_crater

### Define some variables for our grid!
##################################
size = 1000 #size of study area in km (x or y length)
spacing = 5 #size of cells in km
zfactor = 1 #set factor to multiply random noise by

## Set up a LandLab model grid object, with random topographic noise on the order of "zfactor"
##################################
xy = int(size/spacing) # number of nodes
mg = make_noisy_surface(xy,spacing,rf=1) # rf is randomness factor
cmap2 = mpl.cm.get_cmap("Spectral").copy().reversed() # pastel blue-red for topography
fig, ax = plt.subplots() #initiate figure
imshow_grid(mg, 'topographic__elevation', cmap=cmap2, colorbar_label='Elevation [km]')
plt.title("Initial topographic surface")
plt.xlabel("X [km]"); plt.ylabel("Y [km]")
plt.show()

## Set the "distribution" as of crater size as weights, and choose a value from those weights
##################################
###  These parameters describe the population frequency for crater diameters:
Kx = 1.0    #Scaling coefficient (Howard, 2007)
delta = 2.0 #km, scaling exponent (Howard, 2007)
minD = int(spacing) * 3 # want it to be greater than spacing, three times the smallest cell
print("Minimum possible crater size:", minD, "km")
maxD = int((size)/4) #max diameter (m), a quarter of the total domain width
print("Maximum possible crater size:", maxD, "km")
NDs = []
for D in range(minD, maxD):
    ND = Kx * D **-delta
    NDs.append(ND)

## Blast the surface with craters!
##################################
Ncraters = 200 #Number of craters to add
#Ncraters = int(input("How many craters do you want to add? ")) # see console to add your own number!
count = 0
rdm.seed(50) #Chose random seed number 50 (this ensures crater locations are same every time)
mg = do_cratering(Ncraters, NDs, minD, maxD, xy, mg, spacing)

## Make figure of final cratered surface:
##################################
##Set the colourmap
cmap1 =  mpl.cm.get_cmap("Greys_r").copy().reversed() # grey for hillshade

#Edit data for improved plotting
hs = mg.calc_hillshade_at_node(elevs='topographic__elevation') #create hillshade file
topo = mg.field_values('node', 'topographic__elevation').reshape((xy, xy))
hill = np.reshape(hs, (xy, xy))

fig, ax = plt.subplots() #initiate figure
img1 = plt.imshow(hill, cmap=cmap1, alpha=1, extent = [0,size, 0, size])
img2 = plt.imshow(topo, cmap=cmap2, alpha=0.6, extent = [0,size, 0, size])
fig.colorbar(img2,ax=ax, label="Elevation [km]")
plt.title("Final cratered surface")
plt.xlabel("X [km]"); plt.ylabel("Y [km]")
plt.show()

## Add one central crater!
##########################
mg = central_crater(mg, 40, xy, spacing)

## Plot:
hs = mg.calc_hillshade_at_node(elevs='topographic__elevation') #create hillshade file
topo = mg.field_values('node', 'topographic__elevation').reshape((xy, xy))
hill = np.reshape(hs, (xy, xy))

fig, ax = plt.subplots() #initiate figure
img1 = plt.imshow(hill, cmap=cmap1, alpha=1, extent = [0,size, 0, size])
img2 = plt.imshow(topo, cmap=cmap2, alpha=0.6, extent = [0,size, 0, size])
fig.colorbar(img2,ax=ax, label="Elevation [km]")
plt.title("Final cratered surface with central crater")
plt.xlabel("X [km]"); plt.ylabel("Y [km]")
plt.show()
