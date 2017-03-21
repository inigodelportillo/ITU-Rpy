import numpy as np
from scipy.interpolate import griddata, RegularGridInterpolator, interp2d

"""
Interpolation methods for the geophysical properties used to compute 
propagation effects. These methods are based on those in Recommendation
ITU-R P.1144-7.
 
References
----------
[1] Guide to the application of the propagation methods of Radiocommunication 
Study Group 3: https://www.itu.int/rec/R-REC-P.1144/en
"""

def is_regular_grid(lats_o, lons_o):
    '''
    Returns whether the grids in lats_o and lons_o are both regular grids or
    not. A grid is regular if the difference (column-wise or row-wise) 
    between consecutive values is constant across the grid.
    
    Parameters:
    -----------
    lats_o : numpy.ndarray
        Grid of latitude coordinates
    lons_o : numpy.ndarray
        Grid of longitude coordinates
        
    Returns:
    --------
        is_regular: boolean
    '''
    Delta_lons = np.unique(np.diff(lons_o, axis=1))
    Delta_lats = np.unique(np.diff(lats_o, axis=0))
    
    return (Delta_lons.size == 1 and (Delta_lons !=0).all() and\
        Delta_lats.size == 1 and (Delta_lats !=0).all())

###############################################################################    
####                    Nearest Neighbour Interpolation                    ####
###############################################################################

def nearest_2D_interpolator(lats_o, lons_o, values):
    '''
    Produces a 2D interpolator function using the nearest value interpolation 
    method. If the grids are regular grids, uses the 
    scipy.interpolate.RegularGridInterpolator, 
    otherwise, scipy.intepolate.griddata
    
    Values can be interpolated from the returned function as follows:
       f = nearest_2D_interpolator(lat_origin, lon_origin, values_origin)
       interp_values = f(lat_interp, lon_interp)
       
    Parameters:
    -----------
    lats_o: numpy.ndarray
        Latitude coordinates of the values usde by the interpolator
    lons_o: numpy.ndarray      
        Longitude coordinates of the values usde by the interpolator
    values: numpy.ndarray
        Values usde by the interpolator
    
    Returns:
    --------
    interpolator: function
        Nearest neighbour interpolator function
    '''
    # Determine if we are dealing with a regular grid
    if is_regular_grid(lats_o, lons_o):
        return _nearest_2D_interpolator_reg(lats_o, lons_o, values)
    else:
        return _nearest_2D_interpolator_arb(lats_o, lons_o, values)
        
def _nearest_2D_interpolator_reg(lats_o, lons_o, values, lats_d, lons_d):
    f = RegularGridInterpolator((np.flipud(lats_o[:,0]), lons_o[0,:]), 
                               np.flipud(values), method='nearest',
                               bounds_error=False)
    return f

def _nearest_2D_interpolator_arb(lats_o, lons_o, values):
    return lambda x: griddata((lats_o.ravel(), lons_o.ravel()), values.ravel(),
             (x[:,0], x[:,1]), 'nearest')


###############################################################################    
####                         Bilinear Interpolation                        ####
###############################################################################

def bilinear_2D_interpolator(lats_o, lons_o, values):
    '''
    Produces a 2D interpolator function using the bilinear interpolation 
    method. If the grids are regular grids, uses the 
    scipy.interpolate.RegularGridInterpolator, 
    otherwise, scipy.intepolate.griddata
    
    Values can be interpolated from the returned function as follows:
       f = nearest_2D_interpolator(lat_origin, lon_origin, values_origin)
       interp_values = f(lat_interp, lon_interp)
       
    Parameters:
    -----------
    lats_o: numpy.ndarray
        Latitude coordinates of the values usde by the interpolator
    lons_o: numpy.ndarray      
        Longitude coordinates of the values usde by the interpolator
    values: numpy.ndarray
        Values usde by the interpolator
    
    Returns:
    --------
    interpolator: function
        Bilinear interpolator function
    '''
    if is_regular_grid(lats_o, lons_o):
        return _bilinear_2D_interpolator_reg(lats_o, lons_o, values)
    else:
        return _bilinear_2D_interpolator_arb(lats_o, lons_o, values)
        
def _bilinear_2D_interpolator_reg(lats_o, lons_o, values):
    f = RegularGridInterpolator((np.flipud(lats_o[:,0]), lons_o[0,:]),
                               np.flipud(values), method='linear',
                               bounds_error=False)
    return f

def _bilinear_2D_interpolator_arb(lats_o, lons_o, values):
    return lambda x: griddata((lats_o.ravel(), lons_o.ravel()), values.ravel(),
             (x[:,0], x[:,1]), 'linear')

###############################################################################    
####                         Bicubic Interpolation                         ####
###############################################################################
    
def bicubic_2D_interpolator(lats_o, lons_o, values):
    '''
    Produces a 2D interpolator function using the bicubic interpolation 
    method. Uses the scipy.intepolate.griddata method.
    
    Values can be interpolated from the returned function as follows:
       f = nearest_2D_interpolator(lat_origin, lon_origin, values_origin)
       interp_values = f(lat_interp, lon_interp)
       
    Parameters:
    -----------
    lats_o: numpy.ndarray
        Latitude coordinates of the values usde by the interpolator
    lons_o: numpy.ndarray      
        Longitude coordinates of the values usde by the interpolator
    values: numpy.ndarray
        Values usde by the interpolator
    
    Returns:
    --------
    interpolator: function
        Bicubic interpolator function
    '''
    return lambda x:  griddata((lats_o.ravel(), lons_o.ravel()), values.ravel(),
             (x[:,0], x[:,1]), 'cubic')
    