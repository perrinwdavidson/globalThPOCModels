# -------------------------------------
#    Thorium Dynamics in the Ocean
# -------------------------------------
#         Make 2D Model Grid
# -------------------------------------
#         perrin w. davidson
#          pwd@uchicago.edu
# -------------------------------------
# Import packages ---------------------
import numpy as np
import netCDF4 as nc
from functions.interp2d import interp2d
import matplotlib.pyplot as plt

# Define grids ------------------------
# Longitude:
longitude = np.arange(
    -179.5,
    180.5,
    1
)

# Latitude:
latitude = np.arange(
    -89.5,
    90.5,
    1
)

# -------- Working Code Below ---------
# Set paths ---------------------------
input_path = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/inputs/'
output_path = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/outputs/'

# Read data ---------------------------
# Set filename:
filename = input_path + 'bathymetry/GEBCO_2014_2D.nc'

# Read data:
bathymetryfile = nc.Dataset(filename)

# Get data:
bathlon = bathymetryfile['lon'][:].data
bathlat = bathymetryfile['lat'][:].data
bathymetry = bathymetryfile['elevation'][:].data

# Make inputs -------------------------
# Sample arrays:
x = bathlat
y = bathlon
z = bathymetry

# Sample indices:
nx = len(x)
ny = len(y)

# Query arrays:
yq, xq = np.meshgrid(
    longitude,
    latitude
)
xq = xq.flatten()
yq = yq.flatten()

# Query indices:
nxq = len(latitude)
nyq = len(longitude)
nq = len(xq)

# Fill value:
fillvalue = -9999

# Interpolate data --------------------
modelbath = interp2d(
    x,
    y,
    z,
    xq,
    yq,
    fillvalue,
    [nx, ny, nq]
)

# Clean data --------------------------
# Reshape:
modelbath = np.reshape(
    modelbath,
    (nxq, nyq),
)
modelbath[modelbath == fillvalue] = np.nan
modelbath = np.transpose(np.flipud(np.rot90(modelbath)))

# Remove land masses:
modelbath[modelbath > 0] = 0

# Make bathymetry positive:
modelbath = np.abs(modelbath)

# Change longitude convention ---------
# Make longitude to [0, 360]:
n = int(len(longitude) / 2)
longitude[longitude < 0] = longitude[longitude < 0] + 360
longitude = np.concatenate(
    (longitude[n:], longitude[:n]),
    axis=0
)

# Update bathymetry:
modelbath = np.concatenate(
    (modelbath[n:, :], modelbath[:n, :]),
    axis=0
)

# Make mask ---------------------------
# Make:
mask = np.array(
    modelbath,
    dtype=bool
)

# Replace zeros:
modelbath[modelbath == 0] = np.nan

# Transpose ---------------------------
modelbath = np.transpose(modelbath)
mask = np.transpose(mask)

# Save data ---------------------------
# Set filename:
filename = output_path + 'make2dgrid/modelgrid2d.nc'

# Open .nc file:
ncfile = nc.Dataset(
    filename,
    mode='w',
    format='NETCDF4_CLASSIC'
)

# Set dimensions:
lon_dim = ncfile.createDimension('lon', len(longitude))
lat_dim = ncfile.createDimension('lat', len(latitude))

# Set title and subtitle:
ncfile.title = 'THOR Model Grid Data'
ncfile.subtitle = 'THorium and ORganic carbon flux Model'

# Create variables - longitude:
lonnc = ncfile.createVariable('lon', np.float64, 'lon')
lonnc.units = 'degrees_east'
lonnc.long_name = 'longitude'

# Create variables - latitude:
latnc = ncfile.createVariable('lat', np.float64, 'lat')
latnc.units = 'degrees_north'
latnc.long_name = 'latitude'

# Create variables - bathymetry:
bathnc = ncfile.createVariable('bathymetry', np.float64, ('lon', 'lat'))
bathnc.units = 'meters'
bathnc.standard_name = 'bathymetry'

# Create variables - mld variance:
masknc = ncfile.createVariable('mask', np.float64, ('lon', 'lat'))
masknc.units = 'boolean'
masknc.standard_name = 'grid_mask'

# Assign data:
lonnc[:] = longitude
latnc[:] = latitude
bathnc[:, :] = modelbath
masknc[:, :] = mask

# Close file:
ncfile.close()

# Plot data ---------------------------
plt.rcParams['figure.figsize'] = (10, 5)
plt.contourf(np.transpose(modelbath))
plt.colorbar(label='meters')
plt.title('Bathymetry')
plt.show()
plt.close()

# End program -------------------------
