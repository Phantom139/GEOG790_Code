# This script creates a plot for Baroclinic Instabilty for the range of values.
# NC Files Can be Obtained From: ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/
# You can "inventory" a NetCDF file by using 'ncdump -b c "infile.nc" > "outfile.cdl"', which can be loaded in notepad, etc...
# Definition: BI = 0.31 * ((f) / (SQRT( (g/PTm) * ((PTu - PTl) / (GHu - GHl)) ))) * ((Vu - Vl) / (GHu - GHl))

# Import plotting, number, pylab tools
import matplotlib
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
	fileTimeStr = validTime.strftime("%m-%d-%Y-%HZ")
	figName = saveDir + "bci_triplot_" + str(fileTimeStr)

	T=np.squeeze(atempnc.variables['air'][time,:,2:35,:]) 
	u=np.squeeze(uwndnc.variables['uwnd'][time,:,2:35,:]) # u-wind on 26 Oct 2010 1200 UTC (Python index value 1194; index 0 is 00Z 1 Oct. 2010)
	v=np.squeeze(vwndnc.variables['vwnd'][time,:,2:35,:]) # v-wind on 26 Oct 2010 1200 UTC (Python index value 1194; index 0 is 00Z 1 Oct. 2010)
	geo_hght=np.squeeze(hgtnc.variables['hgt'][time,:,2:35,:]) # Geo. Hght on 26 Oct 2010 1200 UTC (Python index value 1194; index 0 is 00Z 1 Oct. 2010)

	# Set up vectors representing 2.5 deg lat, lon values for NCEP/NCAR reanalysis data:
	lat=np.linspace(5,85,33)
	lon=np.linspace(-177.5,180,144)
	p=uwndnc.variables['level'][:] # Vector of isobaric levels (17)

	# Step 1: Define Constants, Easily Calculated Fields (I.E: Coriolis Parameter, Wind Velocity)
	## Coriolis Parameter (f = 2OM * sin(phi))
	Omega = 7.29*math.pow(10,-5)
	TwoOmega = 2*Omega

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
	Shear = np.zeros((ilat, ilon))
	N = np.zeros((ilat, ilon))
	
	matplotlib.rcParams.update({'font.size': 8})

	for i in range(0,ilat):
		corPar = TwoOmega*math.sin((lat[i]) * (np.pi / 180))
		for j in range(0,ilon):	
			PTDif = PotTempHigh[0,0,i,j] - PotTempLow[0,0,i,j]	
			GeoDif = geoHgtHigh[0,0,i,j] - geoHgtLow[0,0,i,j]
			Root = math.sqrt((gConst/PotTempMid[0,0,i,j]) * (PTDif / GeoDif))
			Outer = ((wVelHigh[0,0,i,j]-wVelLow[0,0,i,j])/(geoHgtHigh[0,0,i,j]-geoHgtLow[0,0,i,j]))
			
			BI[i,j] = 0.31 * (corPar / Root) * (Outer) * 86400
			Shear[i,j] = Outer * (86400)
			N[i,j] = Root * 86400

	# Step 4: Plot...
	fig = plt.figure()
	
	m1 = plt.subplot2grid((2,2), (0,0), colspan=2)
	m1 = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
	m1.drawcoastlines()
	m1.drawstates()
	m1.drawcountries()
	m1.drawmapboundary()
	pRange1 = np.linspace(0.25, 2.5, 15, endpoint=True)
	ny1 = BI.shape[0] 
	nx1 = BI.shape[1]
	lons1, lats1 = m1.makegrid(nx1, ny1)
	x1, y1 = m1(lons1, lats1)
	cs1 = m1.contourf(x1, y1, BI, pRange1, cmap=plt.cm.jet)
	cbar1 = m1.colorbar(cs1,location='bottom',pad="5%")
	cbar1.set_label('day^-1')

	m2 = plt.subplot2grid((2,2), (1,0))	
	m2 = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
	m2.drawcoastlines()
	m2.drawstates()
	m2.drawcountries()
	m2.drawmapboundary()
	pRange2 = np.linspace(500, 2000, 15, endpoint=True)
	ny2 = N.shape[0] 
	nx2 = N.shape[1]
	lons2, lats2 = m2.makegrid(nx2, ny2)
	x2, y2 = m2(lons2, lats2)
	cs2 = m2.contourf(x2, y2, N, pRange2, cmap=plt.cm.jet)
	cbar2 = m2.colorbar(cs2,location='bottom',pad="5%")
	cbar2.set_label('day^-1')	
	
	m3 = plt.subplot2grid((2,2), (1,1))	
	m3 = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
	m3.drawcoastlines()
	m3.drawstates()
	m3.drawcountries()
	m3.drawmapboundary()
	pRange3 = np.linspace(1, 1000, 15, endpoint=True)
	ny3 = Shear.shape[0] 
	nx3 = Shear.shape[1]
	lons3, lats3 = m3.makegrid(nx3, ny3)
	x3, y3 = m2(lons3, lats3)
	cs3 = m2.contourf(x3, y3, Shear, pRange3, cmap=plt.cm.jet)
	cbar3 = m3.colorbar(cs3,location='bottom',pad="5%")
	cbar3.set_label('m/day')		
	
	# Finalize the plot
	timeFormat = "%a %b %d %Y %H:%M"
	title = "Baroclinic Instability (" + validTime.strftime(timeFormat) + ")"
	plt.suptitle(title)

	#plt.clabel(fig2_plt, fig2_plt.levels, inline=False, fmt='%r', fontsize=4)

	plt.savefig(figName)
	plt.close()