# This script reads-in CMIP5 data in netcdf format and is used to explore various output data.

# Import relevant packages; many of these come with Anaconda Python 2.7 version, but you will probably have to install the netCDF4 package.  This is installed on the met lab computers in Davis Hall.
import matplotlib.pyplot as plt
import numpy as np
import pylab as py
import os
import csv
from netCDF4 import Dataset # This is important for reading in netCDF4 files below
import math

# Change directory to where data is stored
os.chdir('D:/Robert Docs/College/NIU/GEOG 790/OCT26_CASE') # Insert directory here or comment out if running script in directory where files are saved

# Assign each nc file to variable to read in data
nc_file_uwnd=Dataset('uwnd.2010.nc') # u-wind data
nc_file_vwnd=Dataset('vwnd.2010.nc') # v-wind data
nc_file_hgt=Dataset('hgt.2010.nc') # Geopotential height data
nc_file_land=Dataset('land.nc') # Landmask

# Load Data (2.5 deg horizontal resolution; 17 isobaric levels)
lon=nc_file_uwnd.variables['lon'][:] # Vector of longitude values (144 longitudes)
lat=nc_file_uwnd.variables['lat'][2:35] # Vector of latitude values; only selecting latitude values ranging from 5N to 85N (Python index values 2:35), since that is the latitude range of the jet ID data (33 latitudes selected)
p=nc_file_uwnd.variables['level'][:] # Vector of isobaric levels (17)

# Load u, v and geopotential height data; dimensions of data are times (1460 6-hrly times starting at 00Z 1 Jan. 2010 and ending at 18Z 31 Dec. 2010) x isobaric levels (17) x latitudes (33 selected; see comment above) x longitudes (144)
u=np.squeeze(nc_file_uwnd.variables['uwnd'][1194,:,2:35,:]) # u-wind on 26 Oct 2010 1200 UTC (Python index value 1194; index 0 is 00Z 1 Oct. 2010)
v=np.squeeze(nc_file_vwnd.variables['vwnd'][1194,:,2:35,:]) # v-wind on 26 Oct 2010 1200 UTC (Python index value 1194; index 0 is 00Z 1 Oct. 2010)
geo_hght=np.squeeze(nc_file_hgt.variables['hgt'][1194,:,2:35,:]) # Geo. Hght on 26 Oct 2010 1200 UTC (Python index value 1194; index 0 is 00Z 1 Oct. 2010)

# Indexing
ilat=len(lat) # 'ilat' represents the total number of latitude points for each longitude
ilon=len(lon) # 'ilon' represents the total number of longitude points for each latitude
iday=31 # 'iday' represents the total number of calendar days in October
ihr=4 # 'ihr' represents the number of 6-hourly periods per day (4 total: 00Z, 06Z, 12Z and 18Z)

# Load landmask
land=np.squeeze(nc_file_land.variables['land'][:,:])
land=land[2:35,:]

# Calculate the magnitude of the wind
mag_wind=(u**2+v**2)**0.5

# Set values here for October calender date and 6-hrly time you wish to look at based on date/time loaded in above with NCEP/NCAR data:
oct_date=26 # Enter calendar date here between 1-31
ioct_date=oct_date-1 # Subtract 1 to get proper Python index value

six_hr_time=12 # Enter '0' for 00Z, '6' for 06Z, '12' for 12Z or '18' for 18Z
i6hr_time=six_hr_time/6 # Divide by 6 so '0' = 00Z, '1' = 06Z, '2' = 12Z and '3' = 18Z

# Load in jet superposition ID Data
os.chdir('D:/Robert Docs/College/NIU/GEOG 790/OCT26_CASE') # Insert directory here or comment out if running script in directory where files are saved
ovrlp_matrix_2d=np.loadtxt('ovrlps_oct_2010_NCEP_python.txt') # Load .txt file of jet superposition ID's
ovrlp_ids=ovrlp_matrix_2d.reshape(ilat,ilon,iday,ihr) # The variable 'ovrlp_ids' will have dimensions lat x lon x # days in October x 4 6-hr periods (00Z, 06Z, 12Z and 18Z)
ovrlp_case=np.squeeze(ovrlp_ids[:,:,ioct_date,i6hr_time]) # This line selects overlap data for all latitude and longitude points for date/time specified earlier

# Load in polar ID Data
os.chdir('D:/Robert Docs/College/NIU/GEOG 790/OCT26_CASE') # Insert directory here or comment out if running script in directory where files are saved
polj_matrix_2d=np.loadtxt('poljs_oct_2010_NCEP_python.txt') # Load .txt file of polar jet ID's
polj_ids=polj_matrix_2d.reshape(ilat,ilon,iday,ihr) # The variable 'ovrlp_ids' will have dimensions lat x lon x # days in October x 4 6-hr periods (00Z, 06Z, 12Z and 18Z)
polj_case=np.squeeze(polj_ids[:,:,ioct_date,i6hr_time]) # This line selects overlap data for all latitude and longitude points for date/time specified earlier

# Load in subtropical jet ID Data
os.chdir('D:/Robert Docs/College/NIU/GEOG 790/OCT26_CASE') # Insert directory here or comment out if running script in directory where files are saved
stj_matrix_2d=np.loadtxt('stjs_oct_2010_NCEP_python.txt') # Load .txt file of subtropical jet ID's
stj_ids=stj_matrix_2d.reshape(ilat,ilon,iday,ihr) # The variable 'ovrlp_ids' will have dimensions lat x lon x # days in October x 4 6-hr periods (00Z, 06Z, 12Z and 18Z)
stj_case=np.squeeze(stj_ids[:,:,ioct_date,i6hr_time]) # This line selects overlap data for all latitude and longitude points for date/time specified earlier

# For loop that shifts all points 180 degrees in the jet ID data to line up with NCEP/NCAR Reanalysis 1 Data for each latitude
for i in range(0,ilat):
    ovrlp_case[i,:]=np.roll(ovrlp_case[i,:],72)
    polj_case[i,:]=np.roll(polj_case[i,:],72)
    stj_case[i,:]=np.roll(stj_case[i,:],72)

# Flip latitude dimension to line up ovrlps with NCEP/NCAR Reanalysis 1 Data
ovrlp_case=np.flipud(ovrlp_case)
polj_case=np.flipud(polj_case)
stj_case=np.flipud(stj_case)

# Plot Data with ID's
contour_lvls_wind=[30,40,50,60,70,80,90,100] # Contour levels for wind at 250 hPa
fig2_plt=plt.contour(lon,lat,land,[0,0],colors=[(0.6,0.6,0.6)]) # Plot basic NCEP/NCAR Reanalysis 1 LandMask (use Basemap tool to plot better one)
fig2_plt=plt.contour(lon,lat,np.squeeze(mag_wind[np.where(p==250),:,:]),contour_lvls_wind,colors=[(0,0,0)],linewidths=2) # Plot wind speed as fill pattern
# plt.colorbar(fig2_plt) # Colorbar; insert as line after any line where you use the "plt.contourf" command
fig2_plt=plt.contour(lon,lat,polj_case,colors=[(0,0.5,0)],linewidths=1) # Plot polar jet ID's
fig2_plt=plt.contour(lon,lat,stj_case,colors=[(1,0,0)],linewidths=1) # Plot subtropical jet ID's
fig2_plt=plt.contour(lon,lat,ovrlp_case,colors=[(0,0,0)],linewidths=1) # Plot jet superposition ID's
plt.xlabel('Longitude (degrees)') # x-axis label
plt.ylabel('Latitude (degrees)') # y-axis label
plt.suptitle('12Z 26 Oct. 2010 250 hPa Wind Speed (m/s) and Jet IDs (polar = green; subtropical = red; superposition = black)') # Plot title
# plt.axis([225,300,20,80]) # Use if you want to specifiy lon/lat boundaries to zoom in on particular region - right now it's set to zoom in on U.S.
plt.show() # Show plot!
