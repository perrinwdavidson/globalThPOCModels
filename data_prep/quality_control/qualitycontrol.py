# -------------------------------------
#    Thorium Dynamics in the Ocean
# -------------------------------------
#           Quality Control
# -------------------------------------
#         perrin w. davidson
#          pwd@uchicago.edu
# -------------------------------------
# Import packages:
import pandas as pd
# import numpy as np
# import scipy as sp

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
                    0
)

# Read in GP15 data:
filename = inputs_basepath + 'radionuclide/gp15/gp15_150621.xlsx'
sheet = 'Sheet1'
gp15_rawData = pd.read_excel(
                   filename,
                   sheet,
                   0
)

# Read in 

# End routine.
