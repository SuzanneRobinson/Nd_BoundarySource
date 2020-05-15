#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 11:57:00 2020

Calculating the source of boundary exchange [Nd] per gridbox 
g(Nd)/m3/yr 
Then converted into kg(Nd)/m3/s *(1E18) : the flux specified in the mod


Files needed for this script

1. A(i,j,k) - a file with the sediment area per gridbox (m2) 
Define the depth to which boundary exchange occurs in m
This will act as a mask for the depth of boundary exchange 

2. A_TOTAL- constant 
The total sediment area m2 : again this will vary with depth specified

3. f(bs)/BOUNDARY_SOURCE - constant 
Tuning parameter, this is the force of boundary source (g(Nd)/yr)

4. V(i,j,k)- a file with the volume of the ocean per gridbox (m3)


S(be) = A(i,j,k)/A(tot)  *   boundarySource   * 1/V(i,j,k)
S(be) = AreaRatio_cube * boundarySource * VolumeRatio_cube


@author: suzie
"""

# Load packages
import numpy as np
import matplotlib.pyplot as plt
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import iris.analysis.cartography
from matplotlib import ticker, cm
import cmocean

# **********************************************************************

""" INSERT DEPTH OF OCEAN AND AREA OF THE SEDIMENT BASED ON THE DEPTH HERE

DEPTH OF SEDIMENT SET TO = [ WHOLE OCEAN (6000 M) ]

"""

# Total area of sediment in boundary source (m2)
# This is calculated in area_ocean_sediment.py
A_TOTAL = 118211932610694.69 # 3000 m depth
#A_TOTAL= 891043820092836.6 # full ocean depth 

# *******************************************************************

""" INSERT THE F(BS), THIS IS THE TUNING PARAMETER
REPRESENTS THE FLUX OF THE BOUNDARY SOURCE 
UNITS: GRAMS OF ND PER YEAR
"""
# f(bs), g(Nd)/yr
BOUNDARY_SOURCE = (5.5*10**9)


# *****************************************************************


def area_ratio(sed_cube):
    
    """ Calculates the proportion of area in each grid box 
    to the total BE sediment area (m2)
    A(i,j,k)/A_TOTAL 
    
    
    """
    depth, lat, lon = np.shape(sed_cube.data)
 
    area_ratio_data = np.zeros((depth, lat, lon))
    
    
    for j in range(0,lon):
        for i in range(0,lat):
            for k in range(0,depth):
                area_ratio_data[k,i,j]= 1.0 * (sed_data[k,i,j]/A_TOTAL)
    
    area_ratio_cube = sed_cube.copy(data=area_ratio_data)
    return area_ratio_cube
    



def volume_ratio(vol_cube):
    
    """ Calculates the inverse of the volume of each gridbox (m^-3)
    1/V(i,j,k)
    
    """
    time,depth,lat,lon =np.shape(vol_cube.data)
    
    volume_ratio_data = np.zeros((time, depth, lat, lon))
    
    
    for j in range(0,lon):
        for i in range(0,lat):
            for k in range(0,depth):
                volume_ratio_data[0,k,i,j] = 1.0 * (1.0/vol_data[0,k,i,j])
    volume_ratio_cube= vol_cube.copy(data=volume_ratio_data)
    return volume_ratio_cube




# LOAD FILES

# Load A(i,j,k)
fname1 = '/Users/suzie/Documents/PhD_2ndYear/Self_Isolation_Work/Git/BoundarySource/input_files/sedimentAreaPerGridbox.nc'
# Load V(i,j,k)
fname2 = '/Users/suzie/Documents/PhD_2ndYear/Self_Isolation_Work/Git/BoundarySource/input_files/oceanVolumePerGridbox.nc'

cube=iris.load_cube(fname1)
cube2=iris.load(fname2)

#print(cube)
#print(cube2)

sed_cube = cube[0] # Sediment area
vol_cube = cube2[0] # Volume gridbox

time,depth,lat,lon =np.shape(vol_cube.data)
depth, lat, lon = np.shape(sed_cube.data)

# mask the zeros in the vol cube to avoid dividng by 0.0 
vol_cube.data = np.ma.masked_equal(vol_cube.data, 0.0)


#print(sed_cube.shape)
#print(vol_cube.shape)

sed_data= sed_cube.data
vol_data=vol_cube.data




# return functions to get area_ratio_cube and volume_ratio_cube
# plot returned cubes


area_ratio_cube = area_ratio(sed_cube)
area_ratio_cube.rename('Area_Ratio_cube')
print(area_ratio_cube.data)
qplt.contourf(area_ratio_cube[0,:,:])
plt.show()

qplt.contourf(area_ratio_cube[10,:,:])
plt.show()

qplt.contourf(area_ratio_cube[16,:,:])
plt.show()
    


volume_ratio_cube = volume_ratio(vol_cube)
volume_ratio_cube.rename('Volume_Ratio_cube')
print(volume_ratio_cube.data)
qplt.contourf(volume_ratio_cube[0,0,:,:])
plt.show()


qplt.contourf(volume_ratio_cube[0,10,:,:])
plt.show()

qplt.contourf(volume_ratio_cube[0,16,:,:])
plt.show()



#*******************************************************
# Calculate the source of boundary exchange
# S(BE)/ source_boundary_cube
# This is the density flux of Nd from the boundary source 
# Units g(Nd)/m3/yr

# ******************************************************


""" Units for the s(be)/ source_boundary_cube here are
                g(Nd)/m3/yr
                
              
"""

source_boundary_cube = vol_cube.copy()
source_boundary_cube = area_ratio_cube[:,:,:] * volume_ratio_cube[0,:,:,:] * BOUNDARY_SOURCE
source_boundary_cube.rename('Boundary_Source_Nd_g_m3_yr)')

#print(source_boundary_cube.shape)
#print(source_boundary_cube.data)
#print(source_boundary_cube)


qplt.pcolormesh(source_boundary_cube[0,:,:])
plt.show()

qplt.pcolormesh(source_boundary_cube[6,:,:])
plt.show()

qplt.pcolormesh(source_boundary_cube[16,:,:])
plt.show()



# ***************************************************

""" Units for the s(be)/ source boundary exchange here are
                
    kg(Nd)/m3/s (scaled 1E18 as consistent with dust and river flux)
                
    Converted grams to kg (1E-3)
    Converted year to seconds (3.21502E-8)
    FAMOUS has 360 model days in a model year
    Scaled by 1E18
    
    These are the units consistent with the mod            
                
"""
# Making 4 dimension cube by filling 4d array with SourceBoundayCube.data
array=np.zeros((time,depth,lat,lon))
array[0,:,:,:]=source_boundary_cube.data 

source_boundary = vol_cube.copy(data=array.data)



source_boundary = (source_boundary[:,:,:,:] * 1E-3 * 3.21502E-8 * 1E18)
source_boundary.rename('Boundary_Source_Nd_kg_m3_s_1E18')
source_boundary.coord("depth_1").bounds = None

#print(source_boundary.shape)
#print(source_boundary.data)
#print(source_boundary)





print('Sum of [Nd] from boundary source (kg/s) (*1E18)',np.nansum(source_boundary.data))


# Save sediment area per gridbox as a netCDF file
iris.save(source_boundary, "/Users/suzie/Documents/PhD_2ndYear/Self_Isolation_Work/Git/BoundarySource/input_files/boundarySource_Nd.nc", fill_value=np.nan)






"""
PLotting the boundary conditions for [Nd] for the boundary source

"""

# Create a plot with 6 subplots for Nd
cmap = cmocean.cm.tempo

fig = plt.figure(figsize=(10, 12))
plt.suptitle("Bounday conditions for [Nd] from boundary exchange kg/$m^{3}$/s (*1E18)", fontsize=16)



# These are subplot grid parameters encoded as a single integer. 
#For example, "111" means "1x1 grid, first subplot" and "234" means "2x3 grid, 4th subplot".
ax1 = fig.add_subplot(321)
ax1.set_title('Surface')
iplt.pcolormesh(source_boundary[0,0,:,:], vmin=0.0, vmax=0.0015, cmap=cmap)
plt.colorbar(ax=ax1, orientation="horizontal", extend='max')
ax1 = plt.gca()
ax1.coastlines()

ax2 = fig.add_subplot(322)
ax2.set_title('67 m depth ')
iplt.pcolormesh(source_boundary[0,5,:,:], vmin=0.0, vmax=1.0, cmap=cmap)
plt.colorbar(ax=ax2, orientation="horizontal", extend='max')
ax2 = plt.gca()
ax2.coastlines()

ax3 = fig.add_subplot(323)
ax3.set_title('447 m depth ')
iplt.pcolormesh(source_boundary[0,10,:,:], vmin=0.0, vmax=1.0, cmap=cmap)
plt.colorbar(ax=ax3, orientation="horizontal", extend='max')
ax3 = plt.gca()
ax3.coastlines()

ax4 = fig.add_subplot(324)
ax4.set_title('4577 m depth')
iplt.pcolormesh(source_boundary[0,18,:,:], vmin=0.0, vmax=1.0, cmap=cmap)
plt.colorbar(ax=ax4, orientation="horizontal", extend='max')
ax4 = plt.gca()
ax4.coastlines()

ax5 = fig.add_subplot(325)
ax5.set_title('5192.65 m depth')
iplt.pcolormesh(source_boundary[0,19,:,:], vmin=0.0, vmax=1.0, cmap=cmap)
plt.colorbar(ax=ax5, orientation="horizontal", extend='max')
ax5 = plt.gca()
ax5.coastlines()

ax6 = fig.add_subplot(326)
ax6.set_title(' cross section')
iplt.pcolormesh(source_boundary[0,:,:,0], vmin=0.0, vmax=0.006, cmap=cmap)
plt.colorbar(ax=ax6, orientation="horizontal", extend='max')



plt.tight_layout()
plt.savefig('/Users/suzie/Documents/PhD_2ndYear/Self_Isolation_Work/Git/BoundarySource/plots/TotalNdBoundarySourceBoundaryConditions.png', bbox_inches='tight',format='png', dpi=300)
plt.show()


