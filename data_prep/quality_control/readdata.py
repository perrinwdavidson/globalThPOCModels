# -------------------------------------
#    Thorium Dynamics in the Ocean
# -------------------------------------
#              Read Data
# -------------------------------------
#         perrin w. davidson
#          pwd@uchicago.edu
# -------------------------------------
# Import packages:
import numpy as np
import pandas as pd

# Configure environment ---------------
# What are your base paths?
inputs_basepath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/inputs/'
outputs_basepath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/outputs/'

# Do you want to remove out-of-bounds data? (set to yes if you do)
boundsCoords = 'no'

# Do you want to find and print undated data?
getnodate = 'no'

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
                   sheet,
)

# Read in EXPORTS18 data:
filename = inputs_basepath + 'radionuclide/exports/exports_v6.xlsx'
sheet = 'Table S1'
exports_rawData = pd.read_excel(
                      filename,
                      sheet,
)

# Read in column names:
filename = inputs_basepath + 'radionuclide/columnnames.csv'
columnnames = pd.read_csv(
                      filename,
                      sep=',',
)
columnlist = columnnames.columns

# Read in new dates:
filename = inputs_basepath + 'radionuclide/adddates.xlsx'
sheet = 'Sheet1'
adddates = pd.read_excel(
                      filename,
                      sheet,
)

# Give serial numbers -----------------
# CGODB:
cgodb_dataset = 1
dataset = np.tile(
              np.array(cgodb_dataset),
              [cgodb_rawData.shape[0], 1]
)
serialNumbers = np.linspace(
                    1,
                    cgodb_rawData.shape[0],
                    cgodb_rawData.shape[0],
                    dtype='int'
)
cgodb_rawData.loc[:, 'dataset'] = dataset
cgodb_rawData.loc[:, 'serialNumber'] = serialNumbers

# GP15:
gp15_rawData.sort_values(
    'station_ID',
    axis=0,
    ascending=True
)
gp15_dataset = 2
dataset = np.tile(
              np.array(gp15_dataset),
              [gp15_rawData.shape[0], 1]
)
serialNumbers = np.linspace(
                    1,
                    gp15_rawData.shape[0],
                    gp15_rawData.shape[0],
                    dtype='int'
)
gp15_rawData.loc[:, 'dataset'] = dataset
gp15_rawData.loc[:, 'serialNumber'] = serialNumbers

# EXPORTS18:
exports_dataset = 3
dataset = np.tile(
              np.array(exports_dataset),
              [exports_rawData.shape[0], 1]
)
serialNumbers = np.linspace(
                    1,
                    exports_rawData.shape[0],
                    exports_rawData.shape[0],
                    dtype='int'
)
exports_rawData.loc[:, 'dataset'] = dataset
exports_rawData.loc[:, 'serialNumber'] = serialNumbers

# Collate data into one array ---------
frames = [cgodb_rawData, gp15_rawData, exports_rawData]
ocean = pd.concat(
            frames,
            axis=0
)
ocean.drop(
    ocean.columns[ocean.columns.str.contains('unnamed', case=False)],
    axis=1,
    inplace=True
)
ocean.reset_index(drop=True, inplace=True)

# Remove uncoordinated data -----------
# Make blank coordinates nans:
ocean["lat_decimal_degrees"].replace("", float("nan"), inplace=True)
ocean["lon_decimal_degrees"].replace("", float("nan"), inplace=True)
ocean["depth(m)"].replace("", float("nan"), inplace=True)

# Remove:
ocean.dropna(subset=["lat_decimal_degrees"], inplace=True)
ocean.dropna(subset=["lon_decimal_degrees"], inplace=True)
ocean.dropna(subset=["depth(m)"], inplace=True)

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
ocean["Month"].replace("", float("nan"), inplace=True)
ocean.dropna(subset=["Month"], inplace=True)

# Subset data -------------------------
# Full ocean:
ocean_full = ocean

# Info ocean:
ocean_id = ocean.iloc[:, :43]
ocean_id[['Dataset', 'Serial_Number']] = ocean[['Dataset', 'Serial_Number']]

# Data ocean:
mask = ['Dataset',
        'Serial_Number',
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
    ocean_data.drop(ocean_data[ocean_data['lat_decimal_degrees'] < -90].index, inplace=True)
    ocean_data.drop(ocean_data[ocean_data['lat_decimal_degrees'] > 90].index, inplace=True)
    ocean_data.drop(ocean_data[ocean_data['lon_decimal_degrees'] < -180].index, inplace=True)
    ocean_data.drop(ocean_data[ocean_data['lon_decimal_degrees'] > 360].index, inplace=True)

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

# Write data --------------------------
# Full ocean:
filename = outputs_basepath + 'readdata/ocean_full.xlsx'
sheetname = 'CGODB'
ocean_full.to_excel(
    filename,
    sheetname,
    index=False
)

# Info ocean:
filename = outputs_basepath + 'readdata/ocean_info.xlsx'
sheetname = 'CGODB'
ocean_id.to_excel(
    filename,
    sheetname,
    index=False
)

# Data ocean:
filename = outputs_basepath + 'readdata/ocean_data.xlsx'
sheetname = 'CGODB'
ocean_data.to_excel(
    filename,
    sheetname,
    index=False
)

# End routine.
