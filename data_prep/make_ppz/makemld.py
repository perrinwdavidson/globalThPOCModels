# -------------------------------------
#    Thorium Dynamics in the Ocean
# -------------------------------------
#              Make MLDs
# -------------------------------------
#         perrin w. davidson
#          pwd@uchicago.edu
# -------------------------------------
# Import libraries --------------------
import numpy as np
import wget
import netCDF4 as nc
import matplotlib.pyplot as plt
import seaborn as sns

# Configure environment ---------------
# What are your base paths?
inputs_basepath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/inputs/'
outputs_basepath = '/Users/perrindavidson/Research/whoi/current/globalThPOCModels/outputs/'

# Do you want to download the dataset? (yes or no)
download_mld = 'yes'

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
file2read = nc.Dataset(
    savefilename,
    'r'
)

# Read data:
mldfull = file2read['mld_smth'][:].data
mlderrfull = file2read['krig_std_dev'][:].data

# Read coordinates:
lonmld = file2read['lon'][:].data
latmld = file2read['lat'][:].data
timemld = file2read['time'][:].data

# Get fill value:
fillvaluesample = file2read['mld_smth'].missing_value
fillvalueerr = file2read['krig_std_dev'].missing_value

# Get mask value:
maskvaluesample = file2read['mld_smth'].mask_value
maskvalueeerr = file2read['krig_std_dev'].mask_value

# Read grid data ----------------------
# Set filename:
filename = outputs_basepath + 'make2dgrid/modelgrid2d.nc'

# Set reading file:
file2read = nc.Dataset(
    filename,
    'r'
)

# Read coordinates:
longrid = file2read['lon'][:].data
latgrid = file2read['lat'][:].data

# Read mask:
maskgrid = np.array(
    file2read['mask'][:].data,
    dtype=bool
)

# Process mld data --------------------
# Update mask with monthly fill values:
mldfull[mldfull == fillvaluesample] = 0
mlderrfull[mlderrfull == fillvalueerr] = 0

# Remove masked values:
mldfull[mldfull == maskvaluesample] = 0
mlderrfull[mlderrfull == maskvalueeerr] = 0

# Make mask:
maskmldfull = np.array(
    mldfull,
    dtype=bool
)

# Correct mld arrays:
mldfull[mldfull == 0] = np.nan
mlderrfull[mlderrfull == 0] = np.nan

# Transpose:
mldfull = mldfull.transpose((2, 1, 0))
mlderrfull = mlderrfull.transpose((2, 1, 0))
maskmldfull = np.transpose(maskmldfull)

# Make coordinate grids:
latmldgrid, lonmldgrid = np.meshgrid(
    latmld,
    lonmld
)

# Make coordinate vectors:
lonmldvec = lonmldgrid.flatten()
latmldvec = latmldgrid.flatten()

# Loop through months -----------------
# Get number of months:
nummonths = len(timemld)

# Loop through months:
imonth = 0  # for imonth in range(nummonths):

# Data management -----------------
# Get month arrays:
mld = mldfull[:, :, imonth]
mlderr = mlderrfull[:, :, imonth]
mldmask = maskmldfull[:, :, imonth]

# Make vectors:
mldvec = mld.flatten()
maskmldvec = mldmask.flatten()
mldvec = mldvec[maskmldvec]
lonmldvecmonth = lonmldvec[maskmldvec]
latmldvecmonth = latmldvec[maskmldvec]

# Write out data ------------------
# Make one vector:
mldout = np.array(
    [mldvec, lonmldvecmonth, latmldvecmonth]
).T

# Write:
np.savetxt(
    'calcvariogram/inputs/mld' + str(imonth + 1) + '.txt',
    mldout,
    fmt='%4.2f',
    delimiter=' ',
)

# ---------- work below ----------

# Normalize data ------------------
# Plot distribution (histogram):
plt.figure(figsize=(5, 5))
sns.histplot(mldvec)
plt.xlabel('MLD (m)')
plt.ylabel('Count')
plt.title('Monthly MLD Distribution')
plt.show()
plt.close()

# Define likelihood function inputs:
data = mldvec
ax = -3
bx = 0
cx = 3
tol = 1.48E-8
n = len(data)

# Minimize function:
lmbdahat = brent(
    data,
    ax,
    bx,
    cx,
    llf,
    tol,
    n
)

# Normalize and subtract mean:
mldvecnormmean, mldmean = normalize0mean(
    lmbdahat,
    data,
    n
)

# Plot normalized distribution (histogram):
plt.figure(figsize=(6, 6))
sns.histplot(mldvecnormmean)
plt.xlabel('MLD (m)')
plt.ylabel('Count')
plt.title('Normalized, 0 Mean Monthly MLD Distribution')
plt.show()
plt.close()

# Bin data ------------------------
# Set conversion:
deg2rad = np.pi / 180

# Set distances:
x = lonmldvecmonth * deg2rad
y = latmldvecmonth * deg2rad
z = mldvecnormmean
m = (n * (n - 1)) / 2

# Calculate pairwise distances:
dist, vals, maxdist, numbins = pdist(
    x,
    y,
    z,
    m,
    n
)

# Make new maxdist, user defined:
maxdist = greatcircle(
    0,
    0,
    0,
    30
)

# Make bins:
binedges, bins = makebins(
    numbins,
    numbins + 1,
    maxdist
)

# Write to binary to be read by FORTRAN:
# Here, we are using single-precision (NumPy type float32)
# with big-endian byte ordering:
dist.astype('>f4').tofile('expvario/inputs/dist.bin')
vals.astype('>f4').tofile('expvario/inputs/vals.bin')
binedges.astype('>f4').tofile('expvario/inputs/binedges.bin')

# End routine.
