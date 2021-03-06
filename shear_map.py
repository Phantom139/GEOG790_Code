# This script creates a plot for Baroclinic Instabilty for the range of values.
# NC Files Can be Obtained From: ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/
# You can "inventory" a NetCDF file by using 'ncdump -b c "infile.nc" > "outfile.cdl"', which can be loaded in notepad, etc...
# Definition: N = SQRT( (g/PTm) * ((PTu - PTl) / (GHu - GHl)) )

# Import plotting, number, pylab tools
import matplotlib.pyplot as plt
import numpy as np
import pylab as py
import os
import csv
from netCDF4 import Dataset # This is important for reading in netCDF4 files below
import math
from datetime import date, datetime, time, timedelta
from mpl_toolkits.basemap import Basemap

# Create our data folder if we need to.
currentFilePath = os.path.realpath(__file__)
currentDir = os.path.dirname(currentFilePath)
trgDir = currentDir + '/Data/'
if not os.path.exists(trgDir):
    os.makedirs(trgDir)	
saveDir = currentDir + '/Figures/'
if not os.path.exists(saveDir):
    os.makedirs(saveDir)		
	
# Load in Ze Data
os.chdir(trgDir) # Insert directory here or comment out if running script in directory where files are saved
atempnc=Dataset('air.2010.nc') # air temperature data (In Units of Kelvin)
uwndnc=Dataset('uwnd.2010.nc') # u-wind data
vwndnc=Dataset('vwnd.2010.nc') # v-wind data
hgtnc=Dataset('hgt.2010.nc') # Geopotential height data

startDay = date(2010, 1, 1) 
startHour = time(0, 0, 0)
startingTimeIndex=1184
endingTimeIndex=1208

for time in range(startingTimeIndex,endingTimeIndex):

	validTime = datetime.combine(startDay, startHour) + timedelta(hours=(6 * time))
	figName = saveDir + "shear_map" + str(time)

	u=np.squeeze(uwndnc.variables['uwnd'][time,:,2:35,:]) # u-wind on 26 Oct 2010 1200 UTC (Python index value 1194; index 0 is 00Z 1 Oct. 2010)
	v=np.squeeze(vwndnc.variables['vwnd'][time,:,2:35,:]) # v-wind on 26 Oct 2010 1200 UTC (Python index value 1194; index 0 is 00Z 1 Oct. 2010)
	geo_hght=np.squeeze(hgtnc.variables['hgt'][time,:,2:35,:]) # Geo. Hght on 26 Oct 2010 1200 UTC (Python index value 1194; index 0 is 00Z 1 Oct. 2010)

	# Set up vectors representing 2.5 deg lat, lon values for NCEP/NCAR reanalysis data:
	lat=np.linspace(5,85,33)
	lon=np.linspace(-177.5,180,144)
	p=uwndnc.variables['level'][:] # Vector of isobaric levels (17)

	## Wind Velocity (Magnitude)
	wVel = (u**2+v**2)**0.5
	wVelLow = wVel[np.where(p==1000),:,:]
	wVelHigh = wVel[np.where(p==500),:,:]

	geoHgtLow=geo_hght[np.where(p==1000),:,:]
	geoHgtHigh=geo_hght[np.where(p==500),:,:]

	## Indexing
	ilat=len(lat) # 'ilat' represents the total number of latitude points for each longitude
	ilon=len(lon) # 'ilon' represents the total number of longitude points for each latitude
	iday=31 # 'iday' represents the total number of calendar days in October
	ihr=4 # 'ihr' represents the number of 6-hourly periods per day (4 total: 00Z, 06Z, 12Z and 18Z)

	Shear = np.zeros((ilat, ilon))

	for i in range(0,ilat):
		for j in range(0,ilon):	
			Outer = ((wVelHigh[0,0,i,j]-wVelLow[0,0,i,j])/(geoHgtHigh[0,0,i,j]-geoHgtLow[0,0,i,j]))			
			Shear[i,j] = Outer * (86400)

	# Step 4: Plot...
	#m = Basemap(llcrnrlon=0,llcrnrlat=5,urcrnrlon=360,urcrnrlat=85,projection='mill')	
	m = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
	m.drawcoastlines()
	m.drawstates()
	m.drawcountries()
	m.drawmapboundary()
	pRange = np.linspace(1, 1000, 15, endpoint=True)

	ny = Shear.shape[0] 
	nx = Shear.shape[1]
	lons, lats = m.makegrid(nx, ny)
	x, y = m(lons, lats)
	cs = m.contourf(x, y, Shear, pRange, cmap=plt.cm.jet)
	cbar = m.colorbar(cs,location='bottom',pad="5%")
	cbar.set_label('m/day')
	
	timeFormat = "%a %b %d %Y %H:%M"
	title = "BC Shear Term (" + validTime.strftime(timeFormat) + ") (Normalized)"
	plt.suptitle(title)

	#plt.clabel(fig2_plt, fig2_plt.levels, inline=False, fmt='%r', fontsize=4)

	plt.savefig(figName)
	plt.close()