import numpy as np
import os
import numbers
import csv

from tempfile import mkdtemp
from joblib import Memory

from io import StringIO
from astropy import units as u

dir_path = os.path.dirname(os.path.realpath(__file__))
dataset_dir = os.path.join(dir_path, '../../data/')

# Create a memory cache to memoize results of some functions
cachedir = mkdtemp()
memory = Memory(cachedir=cachedir, verbose=0)


def load_data(path, is_text=False, **kwargs):
    if is_text:
        data = np.loadtxt(path, dtype=np.string_, delimiter=',', **kwargs)
    else:
        data = np.genfromtxt(path, dtype=float, delimiter=',', **kwargs)
    return data


def prepare_input_array(array):
    if isinstance(array, numbers.Number):
        return np.array([[array]])
    elif isinstance(array, list):
        if isinstance(array[0], list):
            return np.array(array)
        else:
            return np.array([array])
    elif isinstance(array, np.ndarray):
        if array.ndim == 1:
            return np.array([array])
        else:
            return array


def prepare_output_array(array, type_input=None):

    if isinstance(array, u.Quantity):
        value = array.value
        unit = array.unit
    else:
        value = array
        unit = None

    if type_input in [numbers.Number, int, float, complex]:
        value = float(value)

    elif type_input is list:
        value = list(value[0])
    else:
        value = value

    if unit is not None:
        return value * unit
    else:
        return value


def prepare_quantity(value, units=None, name_val=None):
    if value is None:
        return None

    if isinstance(value, u.Quantity):
        if units in [u.K, u.deg_C, u.Kelvin, u.Celsius, u.imperial.deg_F]:
            return value.to(units, equivalencies=u.temperature()).value
        else:
            return value.to(units).value

    elif isinstance(value, numbers.Number) and units is not None:
        return value
    elif isinstance(value, np.ndarray) and units is not None:
        return value
    else:
        raise ValueError('%s has not the correct format. It must be a value,'
                         'sequence, array, or a Quantity with %s units' %
                         (name_val, str(units)))


def format_text_file_ITU(path_values, path_lon, path_lat):
    lon = read_raw_file_ITU(path_lon)
    lat = read_raw_file_ITU(path_lat)
    if lon[0:-1] == 360:
        lon = lon[:, :-1]

#    mode = stats.mode(lon, axis=0)[0]
#    idx = np.where(mode == 180)[1]
#    lon[lon > 180] = lon[lon>180] - 360
    idx = 0
    lon2 = np.roll(lon, idx, axis=1)
    lat2 = np.roll(lat, idx, axis=1)

    np.savetxt(path_lon, lon2, '%.2f', delimiter=',')
    np.savetxt(path_lat, lat2, '%.2f', delimiter=',')

    for path in path_values:
        values = read_raw_file_ITU(path)
        values2 = np.roll(values, idx, axis=1)
        np.savetxt(path, values2, '%.6g', delimiter=',')


def compute_distance_earth_to_earth(lat_p, lon_p, lat_grid, lon_grid):
    '''
    Compute the distance between a point (P) in (lat_s, lon_s) and a matrix of
    latitude and longitudes (lat_grid, lon_grid)

    Parameters
    ----------
        lat_p : number
            latitude projection of the point P (degrees)
        lon_p : number
            longitude projection of the point P (degrees)
        lat_grid : number, sequence of np.ndarray
            Grid of latitude points to which compute the distance (degrees)
        lon_grid : number, sequence of np.ndarray
            Grid of longitude points to which compute the distance (degrees)

    Returns
    -------
        d : numpy array
        Distance between the point P and each point in (lat_grid, lon_grid)
        (km)

    References:
    This is based on the Haversine formula
    '''
    RE = 6371.0  # Radius of the Earth, km

    lat1 = np.deg2rad(lat_grid)
    lat2 = np.deg2rad(lat_p)
    lon1 = np.deg2rad(lon_grid)
    lon2 = np.deg2rad(lon_p)

    dlat = lat2-lat1
    dlon = lon2-lon1

    # Compute the distance
    a = np.clip((np.sin(dlat / 2.0))**2 + np.cos(lat1) * np.cos(lat2) *
                (np.sin(dlon / 2))**2, -1, 1)
    c = 2 * np.arcsin(np.sqrt(a))
    d = RE * c
    return d


def regular_lat_lon_grid(resolution_lat=1, resolution_lon=1, lon_start_0=False,
                         lat_min=-90, lat_max=90, lon_min=-180, lon_max=180):
    '''
    Build latitude and longitude coordinate matrix with resolution
    resolution_lat, resolution_lon

    Parameters
    ----------
        resolution_lat, number
            Resolution for the latitude axis (deg)
        resolution_lon, number
            Resolution for the longitude axis (deg)

    Returns
    -------
        lat, numpy.array :
            Grid of coordinates of the latitude point
        lon, numpy.array :
            Grid of coordinates of the latitude point
    '''
    if lon_start_0:
        lon, lat = np.meshgrid(np.arange(lon_min + 180.0, lon_max + 180.0,
                                         resolution_lon),
                               np.arange(lat_max, lat_min, - resolution_lat))
    else:
        lon, lat = np.meshgrid(np.arange(lon_min, lon_max, resolution_lon),
                               np.arange(lat_max, lat_min, - resolution_lat))

    return lat, lon


def elevation_angle(h, lat_s, lon_s, lat_grid, lon_grid):
    '''
    Compute the elevation angle between a satellite located in an orbit
    at height h and located above coordinates (lat_s, lon_s) and a matrix of
    latitude and longitudes (lat_grid, lon_grid)

    Args:
        h : Altitude of the satellite (km)
        lat_s : latitude of the projection of the satellite (degrees)
        lon_s : longitude of the projection of the satellite (degrees)
        lat_grid : Grid of latitude points to which compute the elevation angle (degrees)
        lon_grid : Grid of longitude points to which compute the elevation angle (degrees)

    Returns:
        elevation : Elevation angle between the satellite and each point in (lat_grid, lon_grid) (degrees)

    References:
    [1] http://www.propagation.gatech.edu/ECE6390/notes/ASD5.pdf - Slides 3 and 4
    '''
    RE = 6371.0     # Radius of the Earth (km)
    rs = RE + h

    # Transform latitude_longitude values to radians
    lat1 = np.deg2rad(lat_grid)
    lat2 = np.deg2rad(lat_s)
    lon1 = np.deg2rad(lon_grid)
    lon2 = np.deg2rad(lon_s)

    # Compute the elevation angle as described in
    gamma = np.arccos(np.clip(np.sin(lat2)*np.sin(lat1)+np.cos(lat1)*np.cos(lat2)*np.cos(lon2-lon1), -1, 1))
    elevation = np.arccos(np.sin(gamma)/np.sqrt(1 + (RE/rs)**2-2*(RE/rs)*np.cos(gamma))) # In radians

    return np.rad2deg(elevation)


def read_raw_file_ITU(path):
    print(path)
    with open(path, 'r') as f:
        lines = f.readlines()
        values = []
        for line in lines:
            line = line.replace('                      NaN', ',NaN')
            line = line.replace('             NaN', ',NaN')
            line = line.replace(',,,,,, NaN', ',NaN')
            line = line.replace(',,,,, NaN', ',NaN')
            line = line.replace(',,,, NaN', ',NaN')
            line = line.replace(',,, NaN', ',NaN')
            line = line.replace('     ', ',')
            line = line.replace('    ', ',')
            line = line.replace('   ', ',')
            line = line.replace('  ', ',')
            line = line.replace(' ', ',')
            line = line.replace('\t', ',')
            line = line.replace(' \n', '')
            line = line.replace(',\n', '')
            if line[0] == ',':
                line = line[1:]

            f = StringIO(line)
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                values.append(row)
        return np.array(values, dtype=float)


if __name__ == '__main__':
    pass
    import glob

    folders = ['453/', '530/', '836/', '837/', '839/', '840/', '1510/',
               '1511/']
    for f in folders:
        for doc in glob.iglob(dataset_dir + f + '*.*'):
            vals = read_raw_file_ITU(doc)
            if 'lat' in doc.lower() or 'lon' in doc.lower():
                np.savetxt(doc, vals, '%.5f', delimiter=',')
            else:
                np.savetxt(doc, vals, '%.14g', delimiter=',')
