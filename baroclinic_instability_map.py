# This script creates a plot for Baroclinic Instabilty for the range of values.
# Definition: BI = 0.31 * ((f) / (SQRT( (g/PTm) * ((PTu - PTl) / (GHu - GHl)) ))) * ((Vu - Vl) / (GHu - GHl))

# PT = T (P0 / P) ^ (R/Cp); P0 = 1000mb, R = 287, Cp = 1004

# Import plotting, number, pylab tools
import matplotlib.pyplot as plt
import numpy as np
import pylab as py
import os
import csv

# Create our data folder if we need to.
currentFilePath = os.path.realpath(__file__)
currentDir = os.path.dirname(currentFilePath)
trgDir = currentDir + '/Data/'
if not os.path.exists(trgDir):
    os.makedirs(trgDir)	

# Set up vectors representing 2.5 deg lat, lon values for NCEP/NCAR reanalysis data:
lat=np.linspace(5,85,33)
lon=np.linspace(-177.5,180,144)

# Step 1: Define Constants, Easily Calculated Fields (I.E: Coriolis Parameter, Wind Velocity)

# Step 2: Calculate Potential Temperature for each point at 1000mb (Lower), 700mb (Middle), and 500mb (Upper)

# Step 3: Calculate Baroclinic Instability

# Step 4: Plot...