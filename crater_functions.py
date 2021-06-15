"""Cratering Functions."""
import numpy as np
from landlab import RasterModelGrid


def weighted_choice_sub(weights):
    ''' randomly generate a number and see which weight number in the input list it falls under,
    return the index of that weight '''
    rnd = np.random.random() * sum(weights)
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
        The multiplier to add to increase or decrease randomness by factor rf. 
        (randomness factor, default = 1)
    '''
    mg = RasterModelGrid((xy,xy), xy_spacing = spacing); #initiate surface; see above for variables
    z = mg.add_zeros('topographic__elevation', at='node') #create an array of zeros for each node of the model grid
    np.random.seed(30) # Keep this constant (e.g., at 30) so the initial randomness it always the same
    z += np.random.rand(mg.number_of_nodes)  # make the noise large enough relative to crater
    mg.at_node["topographic__elevation"] *= rf 
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
    radius = diameter / 2

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
    inDepth = H2H1 + H1*((2*incrater)/(diameter))**m
    mg.at_node['topographic__elevation'][in_idx] += inDepth

    outcrater = d[d > radius]
    out_idx = np.where(d > radius)[0]
    # equation for outside the crater (ejecta!)
    outDepth = H2*((2*outcrater)/(diameter))**-n
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

    for i in range(Ncraters):  # For N number of craters
        a = weighted_choice_sub(NDs)
        diameter = list(range(minD, maxD))[a]
        cratercenter = (np.random.randint(1, xy*spacing), np.random.randint(1, xy*spacing))
        d = mg.calc_distances_of_nodes_to_point(cratercenter)

        crater_depth(d, diameter, mg, d_ref=7)

    return mg
