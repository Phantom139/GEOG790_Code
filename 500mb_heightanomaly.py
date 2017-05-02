# This script reads-in CMIP5 data in netcdf format and is used to explore various output data.

# Import relevant packages; many of these come with Anaconda Python 2.7 version, but you will probably have to install the netCDF4 package.  This is installed on the met lab computers in Davis Hall.
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

matplotlib.rcParams.update({'savefig.dpi': 300})

################
# READ IN DATA #
################
startYear = 1981
endYear = 2010
fList = []
for y in range(startYear, endYear):
	fName = "hgt." + str(y) + ".nc"
	fList.append(Dataset(fName))

# Load Data (2.5 deg horizontal resolution; 17 isobaric levels)
lon = fList[0].variables['lon'][:]   # Vector of longitude values (144 longitudes)
lat = fList[0].variables['lat'][:]   # Vector of latitude values; only selecting latitude values ranging from 5N to 85N (Python index values 2:35), since that is the latitude range of the jet ID data (33 latitudes selected)
p   = fList[0].variables['level'][:] # Vector of isobaric levels (17)

###########################################################
# DECLARING ARRAYS AND VARIABLES, PERFORMING CALCULATIONS #
###########################################################

startDay = date(2010, 1, 1) 
startHour = time(0, 0, 0)
startingTimeIndex=1168
endingTimeIndex=1201

geo_hght = []

for time in range(startingTimeIndex,endingTimeIndex):   
	validTime = datetime.combine(startDay, startHour) + timedelta(hours=(6 * time))
	figName = saveDir + "500hPa_heightanomaly" + str(time)
    # Set 850hPa level by letting k = 2
	k = 5 #k[5] = 500 hPa
	i = 0
	for y in range(startYear, endYear):
		geo_hght.append(np.squeeze(fList[i].variables['hgt'][:,k,:,:]))
		i = i + 1
	
    # Indexing
	ilat=len(lat) # 'ilat' represents the total number of latitude points for each longitude
	ilon=len(lon) # 'ilon' represents the total number of longitude points for each latitude
	ilvl=len(p)   # 'ilvl' represents the total number of pressure level points for each latitude
	iday=31       # 'iday' represents the total number of calendar days in October
	ihr=4         # 'ihr' represents the number of 6-hourly periods per day (4 total: 00Z, 06Z, 12Z and 18Z)
	
	#############################################
	# CALCULATING HEIGHT AVERAGES AND ANOMALIES #
	#############################################

	# Setting increments for averaging
	startingaveragingTimeIndex = 1092
	endingaveragingTimeIndex   = 1216
	inputcount = endingaveragingTimeIndex - startingaveragingTimeIndex

	# Matrices full of "-99999" values to be filled with real data below in the for loop
	heighttotals  = np.ones([inputcount,ilat,ilon])*0 # Geopot. height totals for 30 yrs of Oct. before taking avg (m)
	heightaverage = np.ones([ilat,ilon])*0 # Geopotential height average for 30 Octobers (m)
	heightanomaly = np.ones([ilat,ilon])*0 # Geopotential height anomaly (m)
	yearcount     = 0
	cYR           = 0
	sum           = 0
	for averagingtime in range(startingaveragingTimeIndex,endingaveragingTimeIndex):
		for j in range(1,len(lat)-1):
			for i in range(1,len(lon)-1):
				for x in range(startYear, endYear):
					if (((cYR+1) % 4) == 0):
						sum = sum + geo_hght[x][averagingtime+4,j,i]
					else:
						sum = sum + geo_hght[x][averagingtime,j,i]			
				heighttotals[yearcount,j,i] = sum				
        yearcount = yearcount + 1

	for j in range(1,len(lat)-1):
		for i in range(1,len(lon)-1):
			heightaverage[j,i] = sum(heighttotals[:,j,i]) / (inputcount * 30)

	for j in range(1,len(lat)-1):
		for i in range(1,len(lon)-1):
			heightanomaly[j,i] = geo_hght[length(geo_hght)][time,j,i] - heightaverage[j,i]



	#####################
	# PLOTTING THE DATA #
	#####################

	m = Basemap(llcrnrlon=120,llcrnrlat=20,urcrnrlon=300,urcrnrlat=70,projection='mill')
	m.drawcoastlines()
	m.drawcountries()
	m.drawmapboundary()
	m.drawstates(linestyle = ':')
	heights = np.linspace(-450,450,37, endpoint=True)
	lons, lats = np.meshgrid(lon,lat)
	x, y = m(lons, lats)
	cs = m.contourf(x, y, heightanomaly, heights, cmap=plt.cm.bwr)
	cbar = m.colorbar(cs,location='bottom',pad="5%")
	cbar.set_label('Based on a 30-year geopotential height mean \n 1981-2010 (meters)')
	labelrange = [-400,-350,-300,-250,-200,-150,-100,-50,50,100,150,200,250,300,350,400]
	cs = m.contour(x, y, heightanomaly, labelrange, linewidths = 0, colors = [(0,0,0)])
	plt.clabel(cs, fmt="%1.0f", fontsize=6)

	# Title
	timeFormat  = "%H00z %d %b %Y"
	title = "500hPa Height Anomalies (" + validTime.strftime(timeFormat) + ")"
	plt.suptitle(title)

	#plt.clabel(fig2_plt, fig2_plt.levels, inline=False, fmt='%r', fontsize=4)   
	plt.savefig(figName)
	plt.close()
