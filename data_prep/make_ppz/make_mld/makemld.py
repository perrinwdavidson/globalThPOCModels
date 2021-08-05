# -------------------------------------
#    Thorium Dynamics in the Ocean
# -------------------------------------
#              Make MLDs
# -------------------------------------
#         perrin w. davidson
#          pwd@uchicago.edu
# -------------------------------------
# Import ------------------------------
import numpy as np
import wget
from netCDF4 import Dataset
from pykrige import OrdinaryKriging
import pykrige.kriging_tools as kt
import matplotlib.pyplot as plt
import gstools as gs
from pykrige.rk import Krige
from sklearn.model_selection import GridSearchCV

# Configure environment ---------------
# What are your base paths?
inputs_basepath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/inputs/'
outputs_basepath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/outputs/'

# Do you want to download the dataset? (yes or no)
download_mld = 'no'

# Do you want to CV the variogram? (yes or no)
cv_variogram = 'yes'

# Download MLD data -------------------
# Set local filename:
savefilename = inputs_basepath + 'mld/mld_dto2.nc'

# Download:
if download_mld == 'yes':

    url = 'http://www.ifremer.fr/cerweb/deboyer/data/mld_DT02_c1m_reg2.0.nc'
    wget.download(
        url=url,
        out=savefilename
    )

# Read mld data -----------------------
# Set reading file:
file2read = Dataset(
    savefilename,
    'r'
)

# Read data:
mld = file2read['mld_smth'][:].data
mlderr = file2read['krig_std_dev'][:].data

# Read coordinates:
lon = file2read['lon'][:].data
lat = file2read['lat'][:].data
time = file2read['time'][:].data

# Read mask:
mask = np.array(
    file2read['mask'][:].data,
    dtype=bool
)

# Get fill value:
fillvaluesample = file2read['mld_smth'].missing_value
fillvalueerr = file2read['krig_std_dev'].missing_value

# Get mask value:
maskvaluesample = file2read['mld_smth'].mask_value
maskvalueeerr = file2read['krig_std_dev'].mask_value

# Read grid data ----------------------
# Set filename:
filename = outputs_basepath + 'make2dgrid/modelmask2d.asc'

# Read data:
maskq, latq, lonq, cellsize, fillvaluequery = kt.read_asc_grid(filename)

# Process mld data --------------------
# Make sample coordinates:
latgrid, longrid = np.meshgrid(lat, lon)
latvec = latgrid.flatten()
lonvec = longrid.flatten()

# Update mask with monthly fill values:
mld[mld == fillvaluesample] = 0
mlderr[mlderr == fillvalueerr] = 0

# Remove masked values:
mld[mld == maskvaluesample] = 0
mlderr[mlderr == maskvalueeerr] = 0

# Transpose:
mld = mld.transpose((2, 1, 0))
mlderr = mlderr.transpose((2, 1, 0))

# Make mask:
mask = np.array(mld, dtype=bool)

# Replace 0s with nans:
mld[mld == 0] = np.nan

# Process grid data -------------------
# Make mask:
maskq = np.array(
    maskq,
    dtype=bool
)

# Make grids:
latqgrid, lonqgrid = np.meshgrid(latq, lonq)

# Make vectors:
lonqvec = lonqgrid.flatten()
latqvec = latqgrid.flatten()
maskqvec = maskq.flatten()

# Mask.
lonqvec = lonqvec[maskqvec]
latqvec = latqvec[maskqvec]

# Flip to match PyKrige definition:
maskq = ~maskq

# Loop through all times --------------
# Preallocate arrays:
mldok = np.empty((len(lonq), len(latq), 12))
mldvar = np.empty((len(lonq), len(latq), 12))
mldok[:] = np.nan
mldvar[:] = np.nan

# Specify 5-fold CV PyKrige variograms:
variomod = {'variogram_models': ['power',
                                 'power',
                                 'gaussian',
                                 'linear',
                                 'gaussian',
                                 'gaussian',
                                 'gaussian',
                                 'power',
                                 'power',
                                 'linear',
                                 'power',
                                 'linear']
            }

# Loop:
for imonth in range(mld.shape[2]):

    # Get month data:
    mldtime = mld[:, :, imonth]
    masktime = ~mask[:, :, imonth]

    # Make masked array:
    mldmasked = np.ma.masked_array(mldtime, masktime)

    # Make vectors:
    mldvectime = mldtime.flatten()
    maskvectime = masktime.flatten()

    # Mask:
    lonvectime = lonvec[~maskvectime]
    latvectime = latvec[~maskvectime]
    mldvectime = mldvectime[~maskvectime]

    # Calculate number of lags:
    numlags = int(np.round(1 + np.log2(len(mldvectime)), 0))

    # Normalize data:
    mldtimenormraw = gs.normalizer.remove_trend_norm_mean(
        (lon, lat),
        mldmasked,
        normalizer=gs.normalizer.BoxCox(),
        mesh_type="structured",
        fit_normalizer=True
    )
    mldtimenorm = mldtimenormraw[0]
    lmbda = mldtimenormraw[1].lmbda
    mldtimenormvec = mldtimenormraw[0].flatten()
    mldtimenormvec = mldtimenormvec[~maskvectime]

    # Choose best variogram model:
    if cv_variogram == 'yes':

        # Set parameters to test:
        param_dict = {
            "variogram_model": ["linear",
                                "exponential",
                                "gaussian",
                                "spherical",
                                "power"],
            "nlags": [numlags],
            "weight": [True],
            "coordinates_type": ["geographic"],
            "exact_values": [True]
        }

        # Make estimator:
        estimator = GridSearchCV(
            Krige(),
            param_dict,
            verbose=True,
            return_train_score=True
        )

        # Set coordinate and data arrays:
        X = np.array([lonvectime, latvectime]).T
        y = mldtimenormvec

        # Estimate:
        estimator.fit(X=X, y=y)

        # Specify best model:
        best_model = estimator.best_params_['variogram_model']

        # Save best model:
        variomod['variogram_models'][imonth] = best_model

    else:

        best_model = variomod['variogram_models'][imonth]

    # Create PyKrige object:
    OK = OrdinaryKriging(
        lonvectime,
        latvectime,
        mldtimenormvec,
        variogram_model=best_model,
        nlags=numlags,
        weight=True,
        verbose=True,
        enable_plotting=True,
        coordinates_type="geographic",
        exact_values=True
    )

    # PyKrige to masked grid:
    fieldnorm, k_varnorm = OK.execute(
        "masked",
        lonq,
        latq,
        maskq,
        backend="C",
        n_closest_points=500
    )

    # De-normalize, retrend, and add back in mean:
    field = gs.normalizer.apply_mean_norm_trend(
        (latq, lonq),
        fieldnorm,
        normalizer=gs.normalizer.BoxCox(lmbda=lmbda),
        mesh_type="structured"
    )
    k_var = gs.normalizer.apply_mean_norm_trend(
        (latq, lonq),
        k_varnorm,
        normalizer=gs.normalizer.BoxCox(lmbda=lmbda),
        mesh_type="structured"
    )

    # Remove out of bounds field data:
    field[np.transpose(maskq)] = np.nan
    k_var[np.transpose(maskq)] = np.nan

    # Plot:
    plt.figure()
    im = plt.imshow(np.flipud(field))
    plt.colorbar(
        im,
        fraction=0.046,
        pad=0.04,
        shrink=0.61,
        label='Mixed Layer Depths [m]'
    )
    plt.title('Month ' + str(imonth + 1))
    plt.show()
    plt.close()

    # Put into array:
    mldok[:, :, imonth] = np.transpose(field)
    mldvar[:, :, imonth] = np.transpose(k_var)

    # Print out time:
    print('Done kriging MLDs for month', str(imonth + 1))

# Save data ---------------------------
# Set filename:
filename = outputs_basepath + 'makemld/mldok.nc'

# Open .nc file:
ncfile = Dataset(
    filename,
    mode='w',
    format='NETCDF4_CLASSIC'
)

# Set dimensions:
lon_dim = ncfile.createDimension('lon', len(lonq))
lat_dim = ncfile.createDimension('lat', len(latq))
time_dim = ncfile.createDimension('time', len(time))

# Set title and subtitle:
ncfile.title = 'THOR Model MLD Data'
ncfile.subtitle = 'THorium and ORganic carbon flux Model'

# Create variables - longitude:
lonnc = ncfile.createVariable('lon', np.float64, 'lon')
lonnc.units = 'degrees_east'
lonnc.long_name = 'longitude'

# Create variables - latitude:
latnc = ncfile.createVariable('lat', np.float64, 'lat')
latnc.units = 'degrees_north'
latnc.long_name = 'latitude'

# Create variables - time:
timenc = ncfile.createVariable('time', np.float64, 'time')
timenc.units = 'days since january 1'
timenc.long_name = 'date_of_year'

# Create variables - mld:
mldnc = ncfile.createVariable('mld', np.float64, ('lon', 'lat', 'time'))
mldnc.units = 'meters'
mldnc.standard_name = 'mixed_layer_depths'

# Create variables - mld variance:
mldvarnc = ncfile.createVariable('variance', np.float64, ('lon', 'lat', 'time'))
mldvarnc.units = 'meters'
mldvarnc.standard_name = 'mixed_layer_depths_variance'

# Create variables - mld mask:
masknc = ncfile.createVariable('mask', np.float64, ('lon', 'lat'))
masknc.units = 'bool'
masknc.standard_name = 'land_mask'

# Assign data:
lonnc[:] = lonq
latnc[:] = latq
timenc[:] = time
mldnc[:, :, :] = mldok
mldvarnc[:, :, :] = mldvar
masknc[:, :] = maskq

# Close file:
ncfile.close()

# End program -------------------------
