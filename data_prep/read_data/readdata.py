# -------------------------------------
#    Thorium Dynamics in the Ocean
# -------------------------------------
#              Read Data
# -------------------------------------
#         perrin w. davidson
#          pwd@uchicago.edu
# -------------------------------------
# Import packages ---------------------
import numpy as np
import pandas as pd
from netCDF4 import Dataset

# Configure environment ---------------
# What are your base paths?
inputs_basepath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/inputs/'
outputs_basepath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/outputs/'

# Do you want to remove out-of-bounds data? (set to yes if you do)
boundsCoords = 'no'

# Do you want to find and print undated data? (set to yes if you do)
getnodate = 'no'

# For both of the above, these should really be done only once.
# After this, the data needed is then read in below, so it is
# cyclical.

# Read in data ------------------------
# Read in CGODB data:
filename = inputs_basepath + 'radionuclide/cgodb/cgodb_220720.xlsx'
sheet = 'metadata&data'
cgodb_rawData = pd.read_excel(
                    filename,
                    sheet
)

# Read in GP15 data:
filename = inputs_basepath + 'radionuclide/gp15/gp15_150621.xlsx'
sheet = 'Sheet1'
gp15_rawData = pd.read_excel(
                   filename,
                   sheet
)

# Read in EXPORTS data:
filename = inputs_basepath + 'radionuclide/exports/exports_101220.xlsx'
sheet = 'Table S1'
exports_rawData = pd.read_excel(
                      filename,
                      sheet
)

# Read in column names:
filename = inputs_basepath + 'radionuclide/format/columnnames.dat'
columnnames = pd.read_csv(
                      filename,
                      sep=','
)
columnlist = columnnames.columns

# Read in new dates:
filename = inputs_basepath + 'radionuclide/date/adddates.dat'
adddates = pd.read_csv(
                      filename,
                      sep=','
)

# Collate data into one array ---------
frames = [
    cgodb_rawData,
    exports_rawData,
    gp15_rawData
]
ocean = pd.concat(
            frames,
            axis=0
)

# Remove unnamed data -----------------
ocean.drop(
    ocean.columns[ocean.columns.str.contains('unnamed', case=False)],
    axis=1,
    inplace=True
)
ocean.reset_index(
    drop=True,
    inplace=True
)

# Remove uncoordinated data -----------
# Make blank coordinates nans:
ocean["lat_decimal_degrees"].replace(
    "",
    float("nan"),
    inplace=True
)
ocean["lon_decimal_degrees"].replace(
    "",
    float("nan"),
    inplace=True
)
ocean["depth(m)"].replace(
    "",
    float("nan"),
    inplace=True
)

# Remove:
ocean.dropna(
    subset=["lat_decimal_degrees"],
    inplace=True
)
ocean.dropna(
    subset=["lon_decimal_degrees"],
    inplace=True
)
ocean.dropna(
    subset=["depth(m)"],
    inplace=True
)

# Give new names ----------------------
ocean.columns = columnlist

# Add in dates ------------------------
# Get no date data:
if getnodate == 'yes':

    # Get data:
    nodatedata = ocean[ocean['Month'].isna()]

    # Print:
    filename = outputs_basepath + 'readdata/ocean_nodate.xlsx'
    sheetname = 'Complete Dataset'
    nodatedata.to_excel(
        filename,
        sheetname,
        index=False
    )

# Add in new dates:
ocean.loc[ocean['Month'].isna(), 'Month'] = adddates['Month']

# Remove undated data:
ocean["Month"].replace(
    "",
    float("nan"),
    inplace=True
)
ocean.dropna(
    subset=["Month"],
    inplace=True
)

# Subset data -------------------------
# Make mask:
mask = ['Dataset',
        'Dataset_Serial_Number',
        'Longitude',
        'Latitude',
        'Depth(m)',
        'Month',
        '238U(dpm/L)',
        'Uncert_238U(dpm/L)',
        'Total_234Th(dpm/L)',
        'Uncert_Total_234Th(dpm/L)',
        'Dissolved_234Th(dpm/L)',
        'Uncert_Dissolved_234Th(dpm/L)',
        'Particulate_234Th_Small(dpm/L)',
        'Uncert_Particulate_234Th_Small(dpm/L)',
        'POC_Small(umol/L)',
        'Uncert_POC_Small(umol/L)',
        'POC:234Th_Ratio_Small(umol/dpm)',
        'Uncert_POC:234Th_Ratio_Small(umol/dpm)',
        'Particulate_234Th_Large(dpm/L)',
        'Uncert_Particulate_234Th_Large(dpm/L)',
        'POC_Large(umol/L)',
        'Uncert_POC_Large(umol/L)',
        'POC:234Th_Ratio_Large(umol/dpm)',
        'Uncert_POC:234Th_Ratio_Large(umol/dpm)']
ocean_data = ocean[mask]
ocean_data = ocean_data.applymap(lambda x: x.strip() if isinstance(x, str) else x)
ocean_data = ocean_data.apply(pd.to_numeric)

# Remove out-of-bounds coordinates:
if boundsCoords == 'yes':
    ocean_data.drop(
        ocean_data[ocean_data['lat_decimal_degrees'] < -90].index,
        inplace=True
    )
    ocean_data.drop(
        ocean_data[ocean_data['lat_decimal_degrees'] > 90].index,
        inplace=True
    )
    ocean_data.drop(
        ocean_data[ocean_data['lon_decimal_degrees'] < -180].index,
        inplace=True
    )
    ocean_data.drop(
        ocean_data[ocean_data['lon_decimal_degrees'] > 360].index,
        inplace=True
    )

# Make all longitude [0, 360]:
lon_index = ocean_data['Longitude'].index[ocean_data['Longitude'] < 0]
ocean_data.loc[lon_index, 'Longitude'] = ocean_data.loc[lon_index, 'Longitude'] + 360

# Add in seasonality ------------------
# Preallocate:
ocean_data[['Season']] = np.nan

# Winter:
ocean_data.loc[(ocean_data['Month'] == 12)
               | (ocean_data['Month'] == 1)
               | (ocean_data['Month'] == 2)
               & (ocean_data['Latitude'] > 0), 'Season'] = 1
ocean_data.loc[(ocean_data['Month'] == 6)
               | (ocean_data['Month'] == 7)
               | (ocean_data['Month'] == 8)
               & (ocean_data['Latitude'] < 0), 'Season'] = 1

# Spring:
ocean_data.loc[(ocean_data['Month'] == 3)
               | (ocean_data['Month'] == 4)
               | (ocean_data['Month'] == 5)
               & (ocean_data['Latitude'] > 0), 'Season'] = 2
ocean_data.loc[(ocean_data['Month'] == 9)
               | (ocean_data['Month'] == 10)
               | (ocean_data['Month'] == 11)
               & (ocean_data['Latitude'] < 0), 'Season'] = 2

# Summer:
ocean_data.loc[(ocean_data['Month'] == 6)
               | (ocean_data['Month'] == 7)
               | (ocean_data['Month'] == 8)
               & (ocean_data['Latitude'] > 0), 'Season'] = 3
ocean_data.loc[(ocean_data['Month'] == 12)
               | (ocean_data['Month'] == 1)
               | (ocean_data['Month'] == 2)
               & (ocean_data['Latitude'] < 0), 'Season'] = 3

# Fall:
ocean_data.loc[(ocean_data['Month'] == 9)
               | (ocean_data['Month'] == 10)
               | (ocean_data['Month'] == 11)
               & (ocean_data['Latitude'] > 0), 'Season'] = 4
ocean_data.loc[(ocean_data['Month'] == 3)
               | (ocean_data['Month'] == 4)
               | (ocean_data['Month'] == 5)
               & (ocean_data['Latitude'] < 0), 'Season'] = 4

# Put season before month:
idx = ocean_data.columns.get_loc('Month')
season = ocean_data['Season']
ocean_data.drop(
    'Season',
    axis=1,
    inplace=True
)
ocean_data.insert(
    idx,
    'Season',
    season
)

# Give final serial number ------------
ocean_data.insert(
    0,
    'Internal_Serial_Number',
    np.arange(
        0,
        ocean_data.shape[0],
        1
    ) + 1
)

# Write data to netCDF ----------------
# Set filename:
filename = outputs_basepath + 'readdata/ocean.nc'

# Open .nc file:
ncfile = Dataset(
    filename,
    mode='w',
    format='NETCDF4_CLASSIC'
)

# Set dimensions:
serial_number_dim = ncfile.createDimension(
    'index',
    ocean_data.shape[0]
)

# Set title and subtitle:
ncfile.title = 'THOR Model Observational Data'
ncfile.subtitle = 'THorium and ORganic carbon flux Model'

# Create variables - serial number:
snnc = ncfile.createVariable(
    'serial_number',
    np.float64,
    'index'
)
snnc.units = 'index'
snnc[:] = ocean_data['Internal_Serial_Number'].to_numpy()

# Create variables - dataset:
dsnc = ncfile.createVariable(
    'dataset',
    np.float64,
    'index'
)
dsnc.units = 'index'
dsnc[:] = ocean_data['Dataset'].to_numpy()

# Create variables - dataset serial number:
dssnnc = ncfile.createVariable(
    'dataset_serial_number',
    np.float64,
    'index'
)
dssnnc.units = 'index'
dssnnc[:] = ocean_data['Dataset_Serial_Number'].to_numpy()

# Create variables - longitude:
lonnc = ncfile.createVariable(
    'longitude',
    np.float64,
    'index'
)
lonnc.units = 'degrees east'
lonnc[:] = ocean_data['Longitude'].to_numpy()

# Create variables - latitude:
latnc = ncfile.createVariable(
    'latitude',
    np.float64,
    'index'
)
latnc.units = 'degrees north'
latnc[:] = ocean_data['Latitude'].to_numpy()

# Create variables - depth:
dpthnc = ncfile.createVariable(
    'depth',
    np.float64,
    'index'
)
dpthnc.units = 'meters'
dpthnc[:] = ocean_data['Depth(m)'].to_numpy()

# Create variables - season:
seasonnc = ncfile.createVariable(
    'season',
    np.float64,
    'index'
)
seasonnc.units = 'season index'
seasonnc[:] = ocean_data['Season'].to_numpy()

# Create variables - month:
mnthnc = ncfile.createVariable(
    'month',
    np.float64,
    'index'
)
mnthnc.units = 'month of year'
mnthnc[:] = ocean_data['Month'].to_numpy()

# Create variables - u238:
u238nc = ncfile.createVariable(
    'total_238u',
    np.float64,
    'index'
)
u238nc.units = 'dpm L-1'
u238nc[:] = ocean_data['238U(dpm/L)'].to_numpy()

# Create variables - u238 uncertainty:
u238errnc = ncfile.createVariable(
    'total_238u_uncert',
    np.float64,
    'index'
)
u238errnc.units = 'dpm L-1'
u238errnc[:] = ocean_data['Uncert_238U(dpm/L)'].to_numpy()

# Create variables - th234:
th234nc = ncfile.createVariable(
    'total_234th',
    np.float64,
    'index'
)
th234nc.units = 'dpm L-1'
th234nc[:] = ocean_data['Total_234Th(dpm/L)'].to_numpy()

# Create variables - th234 uncertainty:
th234errnc = ncfile.createVariable(
    'total_234th_uncert',
    np.float64,
    'index'
)
th234errnc.units = 'dpm L-1'
th234errnc[:] = ocean_data['Uncert_Total_234Th(dpm/L)'].to_numpy()

# Create variables - dissolved 234th:
dissnc = ncfile.createVariable(
    'dissolved_234th',
    np.float64,
    'index'
)
dissnc.units = 'dpm L-1'
dissnc[:] = ocean_data['Dissolved_234Th(dpm/L)'].to_numpy()

# Create variables - dissolved 234th uncertainty:
disserrnc = ncfile.createVariable(
    'dissolved_234th_uncert',
    np.float64,
    'index'
)
disserrnc.units = 'dpm L-1'
disserrnc[:] = ocean_data['Uncert_Dissolved_234Th(dpm/L)'].to_numpy()

# Create variables - small particulate 234th:
partsnc = ncfile.createVariable(
    'small_particulate_234th',
    np.float64,
    'index'
)
partsnc.units = 'dpm L-1'
partsnc[:] = ocean_data['Particulate_234Th_Small(dpm/L)'].to_numpy()

# Create variables - small particulate 234th uncertainty:
partserrnc = ncfile.createVariable(
    'small_particulate_234th_uncert',
    np.float64,
    'index'
)
partserrnc.units = 'dpm L-1'
partserrnc[:] = ocean_data['Uncert_Particulate_234Th_Small(dpm/L)'].to_numpy()

# Create variables - poc small:
pocsnc = ncfile.createVariable(
    'small_poc',
    np.float64,
    'index'
)
pocsnc.units = 'umol L-1'
pocsnc[:] = ocean_data['POC_Small(umol/L)'].to_numpy()

# Create variables - poc small uncertainty:
pocserrnc = ncfile.createVariable(
    'small_poc_uncert',
    np.float64,
    'index'
)
pocserrnc.units = 'umol L-1'
pocserrnc[:] = ocean_data['Uncert_POC_Small(umol/L)'].to_numpy()

# Create variables - small ratio:
ratsnc = ncfile.createVariable(
    'small_th234_poc_ratio',
    np.float64,
    'index'
)
ratsnc.units = 'umol dpm-1'
ratsnc[:] = ocean_data['POC:234Th_Ratio_Small(umol/dpm)'].to_numpy()

# Create variables - small ratio uncert:
ratserrnc = ncfile.createVariable(
    'small_th234_poc_ratio_uncert',
    np.float64,
    'index'
)
ratserrnc.units = 'umol dpm-1'
ratserrnc[:] = ocean_data['Uncert_POC:234Th_Ratio_Small(umol/dpm)'].to_numpy()

# Create variables - large particulate 234th:
partlnc = ncfile.createVariable(
    'large_particulate_234th',
    np.float64,
    'index'
)
partlnc.units = 'dpm L-1'
partlnc[:] = ocean_data['Particulate_234Th_Large(dpm/L)'].to_numpy()

# Create variables - large particulate 234th uncertainty:
partlerrnc = ncfile.createVariable(
    'large_particulate_234th_uncert',
    np.float64,
    'index'
)
partlerrnc.units = 'dpm L-1'
partlerrnc[:] = ocean_data['Uncert_Particulate_234Th_Large(dpm/L)'].to_numpy()

# Create variables - poc large:
poclnc = ncfile.createVariable(
    'large_poc',
    np.float64,
    'index'
)
poclnc.units = 'umol L-1'
poclnc[:] = ocean_data['POC_Large(umol/L)'].to_numpy()

# Create variables - poc large uncertainty:
poclnc = ncfile.createVariable(
    'large_poc_uncert',
    np.float64,
    'index'
)
poclnc.units = 'umol L-1'
poclnc[:] = ocean_data['Uncert_POC_Large(umol/L)'].to_numpy()

# Create variables - large ratio:
ratlnc = ncfile.createVariable(
    'large_th234_poc_ratio',
    np.float64,
    'index'
)
ratlnc.units = 'umol dpm-1'
ratlnc[:] = ocean_data['POC:234Th_Ratio_Large(umol/dpm)'].to_numpy()

# Create variables - large ratio uncert:
ratlerrnc = ncfile.createVariable(
    'large_th234_poc_ratio_uncert',
    np.float64,
    'index'
)
ratlerrnc.units = 'umol dpm-1'
ratlerrnc[:] = ocean_data['Uncert_POC:234Th_Ratio_Large(umol/dpm)'].to_numpy()

# Close file:
ncfile.close()

# End routine.
