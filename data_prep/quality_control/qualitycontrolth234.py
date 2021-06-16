# -------------------------------------
#    Thorium Dynamics in the Ocean
# -------------------------------------
#         Quality Control Th-234
# -------------------------------------
#         perrin w. davidson
#          pwd@uchicago.edu
# -------------------------------------
# Import packages:
import pandas as pd
import wget
import numpy as np
import netCDF4

# Configure environment:
# What are your base paths?
inputs_basepath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/inputs/'
outputs_basepath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/outputs/'

# Do you want to download WOA18 data?
downloadWOA18 = 'yes'

# Read in ocean info:
filename = outputs_basepath + 'readdata/ocean_info.xlsx'
sheet = 'CGODB'
ocean_info = pd.read_excel(
                filename,
                sheet
)

# Read in ocean data:
filename = outputs_basepath + 'readdata/ocean_data.xlsx'
sheet = 'CGODB'
ocean_data = pd.read_excel(
                filename,
                sheet
)

# Get salinity data:
timeDim = 12
basefilename = 'https://www.ncei.noaa.gov/data/oceans/ncei/woa/salinity/decav/1.00/'
for iMonth in range(timeDim):
    # Get filenames:
    if iMonth <= 8:
        url = basefilename + 'woa18_decav_s0' + str(iMonth+1) + '_01.nc'
        savefilename = inputs_basepath + 'woa18/woa18_sal_0' + str(iMonth+1) + '.nc'
    else:
        url = basefilename + 'woa18_decav_s' + str(iMonth+1) + '_01.nc'
        savefilename = inputs_basepath + 'woa18/woa18_sal_' + str(iMonth+1) + '.nc'

    # Download data:
    if downloadWOA18 == 'yes':
        wget.download(
            url=url,
            out=savefilename
        )

    # Read salinity data:
    file2read = netCDF4.Dataset(savefilename, 'r')
    if iMonth == 0:
        lonDim = file2read.dimensions['lon'].size
        latDim = file2read.dimensions['lat'].size
        depthDim = file2read.dimensions['depth'].size
        salinity = np.ndarray(shape=(lonDim, latDim, depthDim, timeDim))




# End program.
