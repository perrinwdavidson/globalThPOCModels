# -------------------------------------
#    Thorium Dynamics in the Ocean
# -------------------------------------
#           Quality Control
# -------------------------------------
#         perrin w. davidson
#          pwd@uchicago.edu
# -------------------------------------
# Import packages:
import numpy as np
import pandas as pd

# Configure environment:
# Do you want to print?
printing = 1

# What are your base paths?
inputs_basepath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/inputs/'
outputs_basepath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/outputs/'

# Read in CGODB data:
filename = inputs_basepath + 'radionuclide/cgodb/cgodb_220720.xlsx'
sheet = 'metadata&data'
cgodb_rawData = pd.read_excel(
                    filename,
                    sheet,
                    na_values='nan'
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

# Give serial numbers to each dataset:
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

# Collate data into one array:
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

# Remove nan coordinates:
ocean.dropna(subset=["lat_decimal_degrees"], inplace=True)
ocean.dropna(subset=["lon_decimal_degrees"], inplace=True)
ocean.dropna(subset=["depth(m)"], inplace=True)

# Subset data:
# Full ocean:
ocean_full = ocean

# Info ocean:
ocean_id = ocean.iloc[:, :43]
ocean_id[['dataset', 'serialNumber']] = ocean[['dataset', 'serialNumber']]

# Data ocean:
mask = ['serialNumber',
        'lon_decimal_degrees',
        'lat_decimal_degrees',
        'total_234Th(dpm/L)',
        'uncert_total234Th',
        'POC/Th_large(umol/dpm)',
        'uncert_POC/Th_large']
ocean_data = ocean[mask]
ocean_data = ocean_data.applymap(lambda x: x.strip() if isinstance(x, str) else x)
ocean_data = ocean_data.apply(pd.to_numeric)

# Remove out-of-bounds coordinates (commented as there are no such data points):
# ocean_data.drop(ocean_data[ocean_data['lat_decimal_degrees'] < -90].index, inplace=True)
# ocean_data.drop(ocean_data[ocean_data['lat_decimal_degrees'] > 90].index, inplace=True)
# ocean_data.drop(ocean_data[ocean_data['lon_decimal_degrees'] < -180].index, inplace=True)
# ocean_data.drop(ocean_data[ocean_data['lon_decimal_degrees'] > 360].index, inplace=True)

# Make all longitude [0, 360]:
lon_index = ocean_data['lon_decimal_degrees'].index[ocean_data['lon_decimal_degrees'] < 0]
ocean_data.loc[lon_index, 'lon_decimal_degrees'] = ocean_data.loc[lon_index, 'lon_decimal_degrees'] + 360

# Write data:
# Full ocean:
filename = outputs_basepath + 'readdata/ocean_full.xlsx'
sheetname = 'CGODB'
ocean_full.to_excel(filename, sheetname)

# Info ocean:
filename = outputs_basepath + 'readdata/ocean_info.xlsx'
sheetname = 'CGODB'
ocean_id.to_excel(filename, sheetname)

# Data ocean:
filename = outputs_basepath + 'readdata/ocean_data.xlsx'
sheetname = 'CGODB'
ocean_data.to_excel(filename, sheetname)

# End routine.
