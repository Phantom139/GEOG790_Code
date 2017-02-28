# This script creates a plot for Baroclinic Instabilty for the range of values.
# NC Files Can be Obtained From: ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/
# You can "inventory" a NetCDF file by using 'ncdump -b c "infile.nc" > "outfile.cdl"', which can be loaded in notepad, etc...
# Definition: BI = 0.31 * ((f) / (SQRT( (g/PTm) * ((PTu - PTl) / (GHu - GHl)) ))) * ((Vu - Vl) / (GHu - GHl))

# Import plotting, number, pylab tools
import matplotlib.pyplot as plt
import numpy as np
import pylab as py
import os
import csv
from netCDF4 import Dataset # This is important for reading in netCDF4 files below
import math

# Create our data folder if we need to.
currentFilePath = os.path.realpath(__file__)
currentDir = os.path.dirname(currentFilePath)
trgDir = currentDir + '/Data/'
if not os.path.exists(trgDir):
    os.makedirs(trgDir)	
	
# Load in Ze Data
os.chdir(trgDir) # Insert directory here or comment out if running script in directory where files are saved
atempnc=Dataset('air.2010.nc') # air temperature data (In Units of Kelvin)
uwndnc=Dataset('uwnd.2010.nc') # u-wind data
vwndnc=Dataset('vwnd.2010.nc') # v-wind data
hgtnc=Dataset('hgt.2010.nc') # Geopotential height data
landnc=Dataset('land.nc') # Landmask

T=np.squeeze(atempnc.variables['air'][1194,:,2:35,:]) 
u=np.squeeze(uwndnc.variables['uwnd'][1194,:,2:35,:]) # u-wind on 26 Oct 2010 1200 UTC (Python index value 1194; index 0 is 00Z 1 Oct. 2010)
v=np.squeeze(vwndnc.variables['vwnd'][1194,:,2:35,:]) # v-wind on 26 Oct 2010 1200 UTC (Python index value 1194; index 0 is 00Z 1 Oct. 2010)
geo_hght=np.squeeze(hgtnc.variables['hgt'][1194,:,2:35,:]) # Geo. Hght on 26 Oct 2010 1200 UTC (Python index value 1194; index 0 is 00Z 1 Oct. 2010)

# Set up vectors representing 2.5 deg lat, lon values for NCEP/NCAR reanalysis data:
lat=np.linspace(5,85,33)
lon=np.linspace(-177.5,180,144)
p=uwndnc.variables['level'][:] # Vector of isobaric levels (17)

# Step 1: Define Constants, Easily Calculated Fields (I.E: Coriolis Parameter, Wind Velocity)
## Coriolis Parameter (f = 2OM * sin(phi))
Omega = 7.29*math.pow(10,-5)
TwoOmega = 2*Omega

## Fetch the landmask for the plot region
land=np.squeeze(landnc.variables['land'][:,:])
land=land[2:35,:]

## Wind Velocity (Magnitude)
wVel = (u**2+v**2)**0.5
wVelLow = wVel[np.where(p==1000),:,:]
wVelHigh = wVel[np.where(p==500),:,:]

geoHgtLow=geo_hght[np.where(p==1000),:,:]
geoHgtHigh=geo_hght[np.where(p==500),:,:]

# Step 2: Calculate Potential Temperature for each point at 1000mb (Lower), 700mb (Middle), and 500mb (Upper)
## PT = T (P0 / P) ^ (R/Cp); P0 = 1000mb, R = 287, Cp = 1004

ROverCp = 287./1004

PHi = float(500)
PMed = float(700)
PLow = float(1000)

PotTempLow = T[np.where(p==1000),:,:]*(1000/PLow)**ROverCp
PotTempMid = T[np.where(p==700),:,:]*(1000/PMed)**ROverCp
PotTempHigh = T[np.where(p==500),:,:]*(1000/PHi)**ROverCp

# Step 3: Calculate Baroclinic Instability
## BI = 0.31 * ((f) / (SQRT( (g/PTm) * ((PTu - PTl) / (GHu - GHl)) ))) * ((Vu - Vl) / (GHu - GHl))
gConst=9.81

## Indexing
ilat=len(lat) # 'ilat' represents the total number of latitude points for each longitude
ilon=len(lon) # 'ilon' represents the total number of longitude points for each latitude
iday=31 # 'iday' represents the total number of calendar days in October
ihr=4 # 'ihr' represents the number of 6-hourly periods per day (4 total: 00Z, 06Z, 12Z and 18Z)

BI = np.zeros((ilat, ilon))

landFlip=np.flipud(land)

for i in range(0,ilat):
	corPar = TwoOmega*math.sin(lat[i])
	for j in range(0,ilon):	
		PTDif = PotTempHigh[0,0,i,j] - PotTempLow[0,0,i,j]	
		GeoDif = geoHgtHigh[0,0,i,j] - geoHgtLow[0,0,i,j]
		Root = math.sqrt((gConst/PotTempMid[0,0,i,j]) * (PTDif / GeoDif))
		Outer = ((wVelHigh[0,0,i,j]-wVelLow[0,0,i,j])/(geoHgtHigh[0,0,i,j]-geoHgtLow[0,0,i,j]))
		
		BI[i,j] = 0.31 * (corPar / Root) * (Outer)

# Step 4: Plot...
BIContours = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
fig2_plt=plt.contour(lon,lat,landFlip,[0,1],colors=[(0.6,0.6,0.6)])
fig2_plt=plt.contour(lon,lat,BI,BIContours,colors=[(255./255,204./255,204./255), (255./255, 102./255, 102./255), (255./255, 0, 0), (204./255, 0, 0), (153./255, 0, 0)])
plt.xlabel('Longitude (degrees)') # x-axis label
plt.ylabel('Latitude (degrees)') # y-axis label
plt.suptitle('Baroclinic Instability')

#plt.clabel(fig2_plt, fig2_plt.levels, inline=False, fmt='%r', fontsize=4)

plt.show()