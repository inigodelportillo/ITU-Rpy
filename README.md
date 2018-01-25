# ITU-Rpy [![GitHub license](https://img.shields.io/badge/license-MIT-lightgrey.svg)](https://raw.githubusercontent.com/Carthage/Carthage/master/LICENSE.md)

A python implementation of the ITU-P R.XXXX Recommendations to compute atmospheric attenuation in slant and horizontal paths.

The propagation loss on an Earth-space path and a horizontal-path, relative to the free-space loss, is the sum of different contributions, namely:  attenuation by atmospheric gases; attenuation by rain, other precipitation and clouds; scintillation and multipath effects; attenuation by sand and dust storms. Each of these contributions has its own characteristics as a function of frequency, geographic location and elevation angle. This package allows for fast, vectorial computation of the different contributions to the atmospheric attenuation. 

## Installation
ITU-Rpy has the followind dependencies: `numpy`, `scipy`, `joblib`, and `astropy`

Using pip, you can install all of them by running:
```
pip install itur
```

and using conda using:
```
conda install itur
```

## ITU-P Recommendations implemented:
The following ITU-P Recommendations are implemented in ITU-Rpy
* **ITU-P R.453-12:** The radio refractive index: its formula and refractivity data
* **ITU-P R.618-12:** Propagation data and prediction methods required for the design of Earth-space telecommunication systems
* **ITU-P R.676-11:** Attenuation by atmospheric gases
* **ITU-P R.835-12:** Reference Standard Atmospheres
* **ITU-P R.836-5:** Water vapour: surface density and total columnar content
* **ITU-P R.837-6:** Characteristics of precipitation for propagation modelling
* **ITU-P R.838-3:** Specific attenuation model for rain for use in prediction methods
* **ITU-P R.839-4:** Rain height model for prediction methods.
* **ITU-P R.840-6:** Attenuation due to clouds and fog 
* **ITU-P R.1144-7:** Interpolation methods for the geophysical properties used to compute propagation effects 
* **ITU-P R.1511-1:** Topography for Earth-to-space propagation modelling
* **ITU-P R.1853-1:** Tropospheric attenuation time series synthesis

## Usage

The following code examples shows the usage of ITU-Rpy
```python
import itur
from astropy import units as u

f = 86 * u.GHz    # Link frequency
D = 1 * u.m       # Size of the receiver antenna
el = 60           # Elevation angle constant of 60 degrees
p = 3             # Percentage of time that attenuation values are exceeded.
	
# Generate a regular grid latitude and longitude points with 0.1 degrees resolution	
lat, lon = regular_lat_lon_grid() 

# Comute the atmospheric attenuation
Att = itur.atmospheric_attenuation_slant_path(lat, lon, f, el, p, D) 
```

The individual models can be accessed in `itur.models`.
Examples for other use cases can be found in the `examples` folder.