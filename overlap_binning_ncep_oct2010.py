# This script loads one month of jet superposition ID .txt data and grids it into an array such that all ID points are marked as a '1' and all non-ID points are marked as a '0'.

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

# Indexing
ilat=len(lat) # 'ilat' represents the total number of latitude points for each longitude
ilon=len(lon) # 'ilon' represents the total number of longitude points for each latitude
iday=31 # 'iday' represents the total number of calendar days being looped through for specified month in 'for' loop below
iyr=1 # 'iyr' represents the total number of years being looped through below
ihr=4 # 'ihr' represents the number of 6-hourly periods per day (4 total: 00Z, 06Z, 12Z and 18Z)

# For loop that reads in each .txt file of jet superposition ID data in the Northern Hemisphere and stores data within array 'overlap_matrix'
for year in range(2010,2011): # Loop through all years of interest (end of range is last year you want plus 1)
    overlap_matrix=np.zeros([ilat,ilon,iday,ihr]) # Array with all '0' values to use to replace jet ID points with '1' values later
    for month in range(10,11): # Loop through all months of interest (end of range is last month you want plus 1) - To deal with memory issues, I only loop through one month
        for day in range(1,32): # Loop through all days of interest (end of range is last day you want plus 1)
            for hr in range(0,24,6): # Loop through 6-hourly time periods of interest (end of range is last 6-hourly time period you want plus 1)
                if year<2000:
                    if day<10 and hr<10:
                        filename="ovrlp-%d%d0%d0%d.txt" % (year-1900,month,day,hr) # NOTE: IF YOU LOOP THROUGH OCT, NOV, OR DEC MONTH, MAKE SURE THERE IS NOT A '0' IN FRONT OF SECOND '%d' ENTRY ON EVERY LINE THAT SAYS 'FILENAME=....'
                    elif day>=10 and hr<10:
                        filename="ovrlp-%d%d%d0%d.txt" % (year-1900,month,day,hr) # NOTE: IF YOU LOOP THROUGH JAN-AUG MONTHS, MAKE SURE THERE IS A '0' IN FRONT OF SECOND '%d' ENTRY ON EVERY LINE THAT SAYS 'FILENAME=....'
                    elif day<10 and hr>=10:
                        filename="ovrlp-%d%d0%d%d.txt" % (year-1900,month,day,hr)
                    elif day>=10 and hr>=10:
                        filename="ovrlp-%d%d%d%d.txt" % (year-1900,month,day,hr)
                    else:
                        print filename
                elif year>=2000 and year<2010:
                    if day<10 and hr<10:
                        filename="ovrlp-0%d%d0%d0%d.txt" % (year-2000,month,day,hr)
                    elif day>=10 and hr<10:
                        filename="ovrlp-0%d%d%d0%d.txt" % (year-2000,month,day,hr)
                    elif day<10 and hr>=10:
                        filename="ovrlp-0%d%d0%d%d.txt" % (year-2000,month,day,hr)
                    elif day>=10 and hr>=10:
                        filename="ovrlp-0%d%d%d%d.txt" % (year-2000,month,day,hr)
                    else:
                        print filename
                elif year==2010:
                    if day<10 and hr<10:
                        filename="ovrlp-%d%d0%d0%d.txt" % (year-2000,month,day,hr)
                    elif day>=10 and hr<10:
                        filename="ovrlp-%d%d%d0%d.txt" % (year-2000,month,day,hr)
                    elif day<10 and hr>=10:
                        filename="ovrlp-%d%d0%d%d.txt" % (year-2000,month,day,hr)
                    elif day>=10 and hr>=10:
                        filename="ovrlp-%d%d%d%d.txt" % (year-2000,month,day,hr)              
                    else:
                        print filename
                else:
                    print [filename,2]
                os.chdir(trgDir) # Change directory to where .txt files are saved (comment out if not necessary)
                reader=csv.reader(open(filename,"rb"),delimiter=',')
                x=list(reader)
                ovrlp_data=np.array(x).astype('double')
                os.chdir(trgDir) # Change directory to where you want to save newly created .txt files (comment out if not necessary)
                ovrlps=ovrlp_data[:,4] # Superposition ID's now in vector form
                ovrlp_find=np.array(ovrlps==10) # Find all elements where superposition ID is present (marked as a "10" in the NCEP/NCAR Reanalysis 1 ID dataset)
                lat_pts=ovrlp_data[ovrlp_find,2] # Vectors of lat, lon and ovrlp data (this line and next two lines)
                lon_pts=ovrlp_data[ovrlp_find,3]
                ovrlp_pts=ovrlp_data[ovrlp_find,4]
                for i in range(0,len(ovrlp_pts)): # For loop that places a '1' in a grid box that has a jet superposition ID associated with it
                    lat_find=np.array(lat==lat_pts[i])
                    lon_find=np.array(lon==lon_pts[i])
                    overlap_matrix[lat_find,lon_find,day-1,hr/6]=ovrlp_pts[i]/10 # Divide by 10 because jet ID's assigned a "10" rather than a "1" for some reason.
    print year # This is to make sure the script is still running, since if you loop through many years, it takes a while
    overlap_matrix_2d=overlap_matrix.reshape(ilat*ilon,iday*ihr) # Reshape array of ID points to 2D matrix to save as a .txt file; dimensions of .txt file are lat/lon x time dimensions merged together
    filename_save="ovrlps_oct_%d_NCEP_python.txt" % (year) # What you want to name your .txt file
    np.savetxt(filename_save,overlap_matrix_2d) # Save .txt file

print year

# Sum total number of superposition ID's for all dates/times; produce a plot below to see if script worked
overlap_sum=np.sum(np.sum(overlap_matrix,axis=3),axis=2)

# Plot results
fig1_plt=plt.contourf(lon,lat,overlap_sum)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.suptitle('Total Number of Superposition IDs over Northern Hemisphere for October 2010')
cb=plt.colorbar(fig1_plt)
plt.show()
