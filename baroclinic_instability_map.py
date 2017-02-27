# This script creates a plot for Baroclinic Instabilty for the range of values.
# NC Files Can be Obtained From: ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/
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
	
# Load in Ze Data
os.chdir(trgDir) # Insert directory here or comment out if running script in directory where files are saved
uwndnc=Dataset('uwnd.2010.nc') # u-wind data
vwndnc=Dataset('vwnd.2010.nc') # v-wind data
hgtnc=Dataset('hgt.2010.nc') # Geopotential height data
landnc=Dataset('land.nc') # Landmask

u=np.squeeze(uwndnc.variables['uwnd'][1194,:,2:35,:]) # u-wind on 26 Oct 2010 1200 UTC (Python index value 1194; index 0 is 00Z 1 Oct. 2010)
v=np.squeeze(vwndnc.variables['vwnd'][1194,:,2:35,:]) # v-wind on 26 Oct 2010 1200 UTC (Python index value 1194; index 0 is 00Z 1 Oct. 2010)
geo_hght=np.squeeze(hgtnc.variables['hgt'][1194,:,2:35,:]) # Geo. Hght on 26 Oct 2010 1200 UTC (Python index value 1194; index 0 is 00Z 1 Oct. 2010)

# Set up vectors representing 2.5 deg lat, lon values for NCEP/NCAR reanalysis data:
lat=np.linspace(5,85,33)
lon=np.linspace(-177.5,180,144)
p=uwndnc.variables['level'][:] # Vector of isobaric levels (17)

# Step 1: Define Constants, Easily Calculated Fields (I.E: Coriolis Parameter, Wind Velocity)
## Coriolis Parameter (f = 2OM * sin(phi))
Omega = 7.29 * math.pow(10, -5)
TwoOmega = 2 * Omega
CorPar = TwoOmega * np.sin(lat)

## Fetch the landmask for the plot region
land=np.squeeze(nc_file_land.variables['land'][:,:])
land=land[2:35,:]

## Wind Velocity (Magnitude)
wVel = (u**2+v**2)**0.5

# Step 2: Calculate Potential Temperature for each point at 1000mb (Lower), 700mb (Middle), and 500mb (Upper)

# Step 3: Calculate Baroclinic Instability

# Step 4: Plot...