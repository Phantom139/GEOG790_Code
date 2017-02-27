# This script loads one month of subtropical ID .txt data and grids it into an array such that all ID points are marked as a '1' and all non-ID points are marked as a '0'.

# Import plotting, number, pylab tools
import matplotlib.pyplot as plt
import numpy as np
import pylab as py
import os
import csv

# Create our data folder if we need to.
currentFilePath = os.path.realpath(__file__)
currentDir = os.path.dirname(full_path)
trgDir = currentDir + '/Data/'
if not os.path.exists(trgDir):
    os.makedirs(trgDir)	

# Set up vectors representing 2.5 deg lat, lon values for NCEP/NCAR reanalysis data:
lat=np.linspace(5,85,33)
lon=np.linspace(-177.5,180,144)

# Indexing
ilat=len(lat) # 'ilat' represents the total number of latitude points for each longitude
ilon=len(lon) # 'ilon' represents the total number of longitude points for each latitude
iday=31 # 'iday' represents the total number of calendar days being looped through for specified month in 'for' loop below
iyr=1 # 'iyr' represents the total number of years being looped through below
ihr=4 # 'ihr' represents the number of 6-hourly periods per day (4 total: 00Z, 06Z, 12Z and 18Z)

# For loop that reads in each .txt file of subtropical ID data in the Northern Hemisphere and stores data within array 'subtropical_matrix'
for year in range(2010,2011): # Loop through all years of interest (end of range is last year you want plus 1)
    subtropical_matrix=np.zeros([ilat,ilon,iday,ihr]) # Array with all '0' values to use to replace jet ID points with '1' values later
    for month in range(10,11): # Loop through all months of interest (end of range is last month you want plus 1) - To deal with memory issues, I only loop through one month
        for day in range(1,32): # Loop through all days of interest (end of range is last day you want plus 1)
            for hr in range(0,24,6): # Loop through 6-hourly time periods of interest (end of range is last 6-hourly time period you want plus 1)
                if year<2000:
                    if day<10 and hr<10:
                        filename="stj-%d%d0%d0%d.txt" % (year-1900,month,day,hr) # NOTE: IF YOU LOOP THROUGH OCT, NOV, OR DEC MONTH, MAKE SURE THERE IS NOT A '0' IN FRONT OF SECOND '%d' ENTRY ON EVERY LINE THAT SAYS 'FILENAME=....'
                    elif day>=10 and hr<10:
                        filename="stj-%d%d%d0%d.txt" % (year-1900,month,day,hr) # NOTE: IF YOU LOOP THROUGH JAN-AUG MONTHS, MAKE SURE THERE IS A '0' IN FRONT OF SECOND '%d' ENTRY ON EVERY LINE THAT SAYS 'FILENAME=....'
                    elif day<10 and hr>=10:
                        filename="stj-%d%d0%d%d.txt" % (year-1900,month,day,hr)
                    elif day>=10 and hr>=10:
                        filename="stj-%d%d%d%d.txt" % (year-1900,month,day,hr)
                    else:
                        print filename
                elif year>=2000 and year<2010:
                    if day<10 and hr<10:
                        filename="stj-0%d%d0%d0%d.txt" % (year-2000,month,day,hr)
                    elif day>=10 and hr<10:
                        filename="stj-0%d%d%d0%d.txt" % (year-2000,month,day,hr)
                    elif day<10 and hr>=10:
                        filename="stj-0%d%d0%d%d.txt" % (year-2000,month,day,hr)
                    elif day>=10 and hr>=10:
                        filename="stj-0%d%d%d%d.txt" % (year-2000,month,day,hr)
                    else:
                        print filename
                elif year==2010:
                    if day<10 and hr<10:
                        filename="stj-%d%d0%d0%d.txt" % (year-2000,month,day,hr)
                    elif day>=10 and hr<10:
                        filename="stj-%d%d%d0%d.txt" % (year-2000,month,day,hr)
                    elif day<10 and hr>=10:
                        filename="stj-%d%d0%d%d.txt" % (year-2000,month,day,hr)
                    elif day>=10 and hr>=10:
                        filename="stj-%d%d%d%d.txt" % (year-2000,month,day,hr)              
                    else:
                        print filename
                else:
                    print [filename,2]
                os.chdir(trgDir) # Change directory to where .txt files are saved (comment out if not necessary)
                reader=csv.reader(open(filename,"rb"),delimiter=',')
                x=list(reader)
                stj_data=np.array(x).astype('double')
                os.chdir(trgDir) # Change directory to where you want to save newly created .txt files (comment out if not necessary)
                stjs=stj_data[:,4] # Subtropical jet ID's now in vector form
                stj_find=np.array(stjs==10) # Find all elements where subtropical jet ID is present (marked as a "10" in the NCEP/NCAR Reanalysis 1 ID dataset)
                lat_pts=stj_data[stj_find,2] # Vectors of lat, lon and stj data (this line and next two lines)
                lon_pts=stj_data[stj_find,3]
                stj_pts=stj_data[stj_find,4]
                for i in range(0,len(stj_pts)): # For loop that places a '1' in a grid box that has a subtropical jet ID associated with it
                    lat_find=np.array(lat==lat_pts[i])
                    lon_find=np.array(lon==lon_pts[i])
                    subtropical_matrix[lat_find,lon_find,day-1,hr/6]=stj_pts[i]/10 # Divide by 10 because jet ID's assigned a "10" rather than a "1" for some reason.
    print year # This is to make sure the script is still running, since if you loop through many years, it takes a while
    subtropical_matrix_2d=subtropical_matrix.reshape(ilat*ilon,iday*ihr) # Reshape array of ID points to 2D matrix to save as a .txt file; dimensions of .txt file are lat/lon x time dimensions merged together
    filename_save="stjs_oct_%d_NCEP_python.txt" % (year) # What you want to name your .txt file
    np.savetxt(filename_save,subtropical_matrix_2d) # Save .txt file

print year

# Sum total number of subtropical jet ID's for all dates/times; produce a plot below to see if script worked
subtropical_sum=np.sum(np.sum(subtropical_matrix,axis=3),axis=2)

# Plot results
fig1_plt=plt.contourf(lon,lat,subtropical_sum)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.suptitle('Total Number of Subtropical Jet IDs over Northern Hemisphere for October 2010')
cb=plt.colorbar(fig1_plt)
plt.show()
