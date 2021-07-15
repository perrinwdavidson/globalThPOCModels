# -------------------------------------
#    Thorium Dynamics in the Ocean
# -------------------------------------
#       Quality Control Th-234
# -------------------------------------
#         perrin w. davidson
#          pwd@uchicago.edu
# -------------------------------------
# Import packages:
import pandas as pd
import wget
import numpy as np
from numpy import sin, cos, arctan2, sqrt, pi
import netCDF4 as nc
import os
import matlab.engine
from scipy.interpolate import RBFInterpolator
from scipy.interpolate import NearestNDInterpolator


# Begin local functions ---------------
def haversine(loc1, loc2):

    # Get coordinates:
    lat1 = loc1[1]
    lon1 = loc1[0]
    lat2 = loc2.iloc[:, 1]
    lon2 = loc2.iloc[:, 0]

    # Convert to radians:
    lon1 = lon1 * pi / 180.0
    lon2 = lon2 * pi / 180.0
    lat1 = lat1 * pi / 180.0
    lat2 = lat2 * pi / 180.0

    # Calculate using haversine:
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (sin(dlat / 2)) ** 2 + cos(lat1) * cos(lat2) * (sin(dlon / 2.0)) ** 2
    c = 2.0 * arctan2(sqrt(a), sqrt(1.0 - a))
    m = 6371.0 * c * 1000

    # Set output:
    return m

# End functions -----------------------

# Configure environment ---------------
# What are your base paths?
inputs_basepath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/inputs/'
outputs_basepath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/outputs/'

# Do you want to download WOA18 data?
downloadWOA18 = 'no'

# Set 238U interpolation type (rbf or nearest):
interptype = 'nearest'

# If rbf, set interp type (of scipy rbf kind):
rbftype = 'thin_plate_spline'

# Start MATLAB engine -----------------
eng = matlab.engine.start_matlab()

# Read 234th data ---------------------
# Read in ocean info:
filename = outputs_basepath + 'readdata/ocean_info.xlsx'
sheet = 'CGODB'
oceaninfo = pd.read_excel(
    filename,
    sheet
)

# Read in ocean data:
filename = outputs_basepath + 'readdata/ocean_data.xlsx'
sheet = 'CGODB'
oceandata = pd.read_excel(
    filename,
    sheet
)

# Get salinity data -------------------
# Set time dimension:
timeDim = 12

# Set fill value:
fillval = 9969209968386869046778552952102584320.00000

# Set base filename for WOA18:
basefilename = 'https://www.ncei.noaa.gov/data/oceans/ncei/woa/salinity/decav/1.00/'

# Loop through each time dimension:
for iMonth in range(timeDim):

    # Get filenames:
    if iMonth <= 8:

        url = basefilename + 'woa18_decav_s0' + str(iMonth + 1) + '_01.nc'
        savefilename = inputs_basepath + 'woa18/woa18_sal_0' + str(iMonth + 1) + '.nc'

    else:

        url = basefilename + 'woa18_decav_s' + str(iMonth + 1) + '_01.nc'
        savefilename = inputs_basepath + 'woa18/woa18_sal_' + str(iMonth + 1) + '.nc'

    # Download data:
    if downloadWOA18 == 'yes':

        # Remove file if exists:
        if os.path.exists(savefilename):
            os.remove(savefilename)

        # Download:
        wget.download(
            url=url,
            out=savefilename
        )

    # Set reading file:
    file2read = nc.Dataset(savefilename, 'r')

    # Read salinity data:
    if iMonth == 0:
        # Make array:
        lonDim = file2read.dimensions['lon'].size
        latDim = file2read.dimensions['lat'].size
        depthDim = file2read.dimensions['depth'].size
        salinity = np.ndarray(shape=(
            timeDim,
            depthDim,
            latDim,
            lonDim)
        )

        # Get coordinate arrays:
        sallon = file2read['lon'][:].data
        sallon[sallon < 0] = sallon[sallon < 0] + 360
        sallat = file2read['lat'][:].data
        saldepth = file2read['depth'][:].data

        # Convert to meters:
        sallatmat = matlab.double([sallat.tolist()])
        sallonmat = matlab.double([sallon.tolist()])
        sallon = np.transpose(np.array(eng.deg2m(
            sallonmat,
            nargout=1))
        )
        sallat = np.transpose(np.array(eng.deg2m(
            sallatmat,
            nargout=1))
        )

        # Make time array:
        saltime = np.linspace(1, 12, 12)

        # Make coordinate meshgrids:
        zv, tv, yv, xv = np.meshgrid(
            saldepth,
            saltime,
            sallat,
            sallon
        )

    # Put salinity data into array:
    sal = file2read['s_an'][:].data
    sal[sal == fillval] = np.nan
    salinity[iMonth, :, :, :] = sal

    # Print out month:
    print(
        'Done downloading salinity from month',
        str(iMonth + 1)
    )

# Remove 0 from salinity --------------
salinity[salinity == 0] = np.nan

# Remove negative 234Th ---------------
oceandata.loc[oceandata['238U(dpm/L)'] < 0,
              '238U(dpm/L)'] = np.nan
oceandata.loc[oceandata['Total_234Th(dpm/L)'] < 0,
              'Total_234Th(dpm/L)'] = np.nan
oceandata.loc[oceandata['Dissolved_234Th(dpm/L)'] < 0,
              'Dissolved_234Th(dpm/L)'] = np.nan
oceandata.loc[oceandata['Particulate_234Th_Small(dpm/L)'] < 0,
              'Particulate_234Th_Small(dpm/L)'] = np.nan
oceandata.loc[oceandata['Particulate_234Th_Large(dpm/L)'] < 0,
              'Particulate_234Th_Large(dpm/L)'] = np.nan

# Correct large values ----------------
oceandata.loc[oceandata['238U(dpm/L)'] > 1000,
              '238U(dpm/L)'] = oceandata.loc[oceandata['238U(dpm/L)'] > 1000,
                                             '238U(dpm/L)'] / 1000
oceandata.loc[oceandata['Total_234Th(dpm/L)'] > 1000,
              'Total_234Th(dpm/L)'] = oceandata.loc[oceandata['Total_234Th(dpm/L)'] > 1000,
                                                    'Total_234Th(dpm/L)'] / 1000
oceandata.loc[oceandata['Dissolved_234Th(dpm/L)'] > 1000,
              'Dissolved_234Th(dpm/L)'] = oceandata.loc[oceandata['Dissolved_234Th(dpm/L)'] > 1000,
                                                        'Dissolved_234Th(dpm/L)'] / 1000
oceandata.loc[oceandata['Particulate_234Th_Small(dpm/L)'] > 10,
              'Particulate_234Th_Small(dpm/L)'] = oceandata.loc[oceandata['Particulate_234Th_Small(dpm/L)'] > 10,
                                                                'Particulate_234Th_Small(dpm/L)'] / 1000
oceandata.loc[oceandata['Particulate_234Th_Large(dpm/L)'] > 10,
              'Particulate_234Th_Large(dpm/L)'] = oceandata.loc[oceandata['Particulate_234Th_Large(dpm/L)'] > 10,
                                                                'Particulate_234Th_Large(dpm/L)'] / 1000

# Make total 234Th --------------------
# Find indices:
idxtot = oceandata['Total_234Th(dpm/L)'].index[oceandata['Total_234Th(dpm/L)'].apply(np.isnan)].values.tolist()
idxdiss = oceandata['Dissolved_234Th(dpm/L)'].index[
    oceandata['Dissolved_234Th(dpm/L)'].apply(np.isfinite)].values.tolist()
idxpartsmall = oceandata['Particulate_234Th_Small(dpm/L)'].index[
    oceandata['Particulate_234Th_Small(dpm/L)'].apply(np.isfinite)].values.tolist()
idxpartlarge = oceandata['Particulate_234Th_Large(dpm/L)'].index[
    oceandata['Particulate_234Th_Large(dpm/L)'].apply(np.isfinite)].values.tolist()

# Concatenate:
idx = [
    idxdiss,
    idxpartsmall,
    idxpartlarge
]

# Find intersection:
makeidx = list(set(idxtot).intersection(*idx))

# Make new 234Th values:
oceandata.loc[makeidx, 'Total_234Th(dpm/L)'] = oceandata.loc[makeidx, 'Dissolved_234Th(dpm/L)'] \
                                               + oceandata.loc[makeidx, 'Particulate_234Th_Small(dpm/L)'] \
                                               + oceandata.loc[makeidx, 'Particulate_234Th_Large(dpm/L)']

# Collect 234Th values ----------------
oceandata.dropna(subset=["Total_234Th(dpm/L)"], inplace=True)
th234 = oceandata

# Exclude large values ----------------
idxlarge = oceandata['Total_234Th(dpm/L)'].index[oceandata['Total_234Th(dpm/L)'] > 10].values.tolist()
if len(idxlarge) != 0:
    th234 = th234.drop(idxlarge)

# Make 238U values from salinity ------
u238 = (0.0786 * salinity) - 0.315

# Fill 238U nan values ----------------
# Find nan indices:
idxu = th234['238U(dpm/L)'].index[th234['238U(dpm/L)'].apply(np.isnan)].values.tolist()

# Replace values:
uLat = th234.loc[idxu, 'Latitude']
uLon = th234.loc[idxu, 'Longitude']
uDepth = th234.loc[idxu, 'Depth(m)']
uNan = th234.loc[idxu, '238U(dpm/L)']

# Get 238U time:
uTime = th234.loc[idxu, 'Month']
uTime = uTime.fillna(0)

# Convert degrees to meters:
uLatmat = matlab.double([uLat.tolist()])
uLonmat = matlab.double([uLon.tolist()])
uLon = np.transpose(np.array(eng.deg2m(uLonmat, nargout=1)))
uLat = np.transpose(np.array(eng.deg2m(uLatmat, nargout=1)))

# Make coordinate query arrays:
xq = np.transpose(np.array(
    [
        uNan,
        uLon[:, 0],
        uLat[:, 0],
        uDepth,
        uTime
    ])
)

# Interpolate 238U data ---------------
for itime in range(timeDim):

    # Get month sample data:
    samples = np.transpose(np.array(
        [
            u238[itime, :, :, :].flatten(),
            xv[itime, :, :, :].flatten(),
            yv[itime, :, :, :].flatten(),
            zv[itime, :, :, :].flatten()
        ])
    )
    nan_array = np.isnan(samples[:, 0])
    not_nan_array = ~ nan_array
    samples = samples[not_nan_array, :]
    d = samples[:, 0]
    y = samples[:, 1:4]

    # Get month query data:
    x = xq[xq[:, 4] == itime + 1, 1:4]

    # Interpolate:
    if interptype == 'rbf':

        vq = RBFInterpolator(
            y,
            d,
            neighbors=1000,
            kernel=rbftype,
        )(x)

    elif interptype == 'nearest':

        vq = NearestNDInterpolator(
            y,
            d,
        )(x)

    # Put into array:
    xq[xq[:, 4] == itime + 1, 0] = np.transpose(vq)

    # Print out month:
    print(
        'Done interpolating 238U from month',
        str(itime + 1)
    )

# Put new 238U values into 234Th dataframe:
th234.loc[idxu, '238U(dpm/L)'] = pd.Series(xq[:, 0], index=idxu)
th234.loc[idxu, 'Uncert_238U(dpm/L)'] = 0.047

# Fill 234Th errors -------------------
# Find uncertainties indices:
idxth = th234['Uncert_Total_234Th(dpm/L)'].index[th234['Uncert_Total_234Th(dpm/L)'].apply(np.isfinite)].values.tolist()

# Get data:
th234datauncert = th234.loc[idxth, :]

# Find nan uncertainties indices:
idxth = th234['Uncert_Total_234Th(dpm/L)'].index[th234['Uncert_Total_234Th(dpm/L)'].apply(np.isnan)].values.tolist()

# Get data:
th234uncert = th234.loc[idxth, :]

# Make sample arrays:
y = th234datauncert[['Latitude', 'Longitude', 'Depth(m)', 'Month']]
d = th234datauncert[['Uncert_Total_234Th(dpm/L)']]
n = y.shape[0]

# Make query arrays:
x = th234uncert[['Latitude', 'Longitude', 'Depth(m)', 'Month']]

# Set radii:
rcoord = eng.deg2m(15.0, nargout=1)
rdepth = 10
rtime = 3

# Loop through data:
numprof = 1
numnnint = 1
for iprof in idxth:

    # Get coords:
    ilat = th234uncert.loc[iprof, 'Latitude']
    ilon = th234uncert.loc[iprof, 'Longitude']
    idepth = th234uncert.loc[iprof, 'Depth(m)']
    imon = th234uncert.loc[iprof, 'Month']

    # Make query array:
    xcoord = np.array([ilon, ilat])

    # Make coordinate arrays:
    loc1 = xcoord
    loc2 = y[['Longitude', 'Latitude']]

    # Calculate coordinate distances:
    coorddists = haversine(loc1, loc2)
    depthdists = np.array(np.abs(idepth - y['Depth(m)']))
    timedists = np.array(np.abs(imon - y['Month']))

    # Specify data:
    iuncert = d.loc[
        (coorddists <= rcoord)
        & (depthdists <= rdepth)
        & (timedists <= rtime), 'Uncert_Total_234Th(dpm/L)'
    ]

    # NNint if none within radius:
    if len(iuncert) == 0:

        # Make query array:
        x = np.array([ilat, ilon, idepth, imon])

        # Interpolate:
        vq = NearestNDInterpolator(
            y,
            d,
        )(x)

        # Extract value:
        vq = vq[0][0]

        # Count:
        numnnint += 1

    else:

        # Average:
        vq = np.nanmean(iuncert)

    # Put into array:
    th234uncert.loc[iprof, 'Uncert_Total_234Th(dpm/L)'] = vq

    # Print out:
    if numprof % 500 == 0:
        print('Done with 234Th error profile', str(numprof))
    numprof += 1

    # Clear:
    vq = np.nan

# Put back into arrays:
th234.loc[idxth, 'Uncert_Total_234Th(dpm/L)'] = th234uncert[['Uncert_Total_234Th(dpm/L)']]

# Fill 238U errors --------------------
# Make zero errors into nans:
th234['Uncert_238U(dpm/L)'].replace(0, float("nan"), inplace=True)

# Find uncertainties indices:
idxu = th234['Uncert_238U(dpm/L)'].index[th234['Uncert_238U(dpm/L)'].apply(np.isfinite)].values.tolist()

# Get data:
u238datauncert = th234.loc[idxu, :]

# Find nan uncertainties indices:
idxu = th234['Uncert_238U(dpm/L)'].index[th234['Uncert_238U(dpm/L)'].apply(np.isnan)].values.tolist()

# Get data:
u238uncert = th234.loc[idxu, :]

# Make sample arrays:
y = u238datauncert[['Latitude', 'Longitude', 'Depth(m)', 'Month']]
d = u238datauncert[['Uncert_238U(dpm/L)']]
n = y.shape[0]

# Make query arrays:
x = u238uncert[['Latitude', 'Longitude', 'Depth(m)', 'Month']]

# Set radii:
rcoord = eng.deg2m(15.0, nargout=1)
rdepth = 10
rtime = 3

# Loop through data:
numprof = 1
numnnint = 1
for iprof in idxu:

    # Get coords:
    ilat = u238uncert.loc[iprof, 'Latitude']
    ilon = u238uncert.loc[iprof, 'Longitude']
    idepth = u238uncert.loc[iprof, 'Depth(m)']
    imon = u238uncert.loc[iprof, 'Month']

    # Make query array:
    xcoord = np.array([ilon, ilat])

    # Make coordinate arrays:
    loc1 = xcoord
    loc2 = y[['Longitude', 'Latitude']]

    # Calculate coordinate distances:
    coorddists = haversine(loc1, loc2)
    depthdists = np.array(np.abs(idepth - y['Depth(m)']))
    timedists = np.array(np.abs(imon - y['Month']))

    # Specify data:
    iuncert = d.loc[
        (coorddists <= rcoord)
        & (depthdists <= rdepth)
        & (timedists <= rtime), 'Uncert_238U(dpm/L)'
    ]

    # NNint if none within radius:
    if len(iuncert) == 0:

        # Make query array:
        x = np.array([ilat, ilon, idepth, imon])

        # Interpolate:
        vq = NearestNDInterpolator(
            y,
            d,
        )(x)

        # Extract value:
        vq = vq[0][0]

        # Count:
        numnnint += 1

    else:

        # Average:
        vq = np.nanmean(iuncert)

    # Put into array:
    u238uncert.loc[iprof, 'Uncert_238U(dpm/L)'] = vq

    # Print out:
    if numprof % 500 == 0:
        print('Done with 238U error profile', str(numprof))
    numprof += 1

    # Clear:
    vq = np.nan

# Put back into arrays:
th234.loc[idxu, 'Uncert_238U(dpm/L)'] = u238uncert[['Uncert_238U(dpm/L)']]

# Exclude 234Th w/ 238U and depth -----
# Make observation arrays:
th234obs = th234['Total_234Th(dpm/L)'] + th234['Uncert_Total_234Th(dpm/L)'].abs()
u238obs = th234['238U(dpm/L)'] + th234['Uncert_238U(dpm/L)'].abs()

# Exclude if <50 & >U238:
th234.drop(th234[(th234obs >= u238obs)
                 & (th234['Depth(m)'] <= 50)].index, inplace=True)

# Exclude if 50><300 & >3:
th234obs = th234obs[th234.index]
u238obs = u238obs[th234.index]
th234.drop(th234[(th234obs >= 3)
                 & (th234['Depth(m)'] > 50)
                 & (th234['Depth(m)'] <= 300)].index, inplace=True)

# Exclude if >300 & >238U:
th234obs = th234obs[th234.index]
u238obs = u238obs[th234.index]
th234.drop(th234[(th234obs >= u238obs)
                 & (th234['Depth(m)'] > 300)].index, inplace=True)

# Exclude 234Th w/ counter ------------
# Set limit:
betacounter = 0.3

# Exclude:
th234obs = th234obs[th234.index]
th234.drop(th234[th234obs <= betacounter].index, inplace=True)

# Exclude 234Th w/ plotting -----------
# Set sns:
removesn = [
    1447,
    3398,
    7727,
    18577,
    21665
]

# Exclude:
for isr in removesn:
    th234 = th234.loc[~(th234['Serial_Number'] == isr), :]

# Write out ---------------------------
filename = outputs_basepath + 'qualitycontrolth234/th234.xlsx'
sheetname = '234Th'
th234.to_excel(
    filename,
    sheetname,
    index=False
)

# End program -------------------------
