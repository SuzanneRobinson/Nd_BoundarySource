#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Suzanne Robinson: Tue Mar 17 12:02:14 2020

Calculates the area of sediment in contact with the ocean in FAMOUS


"""

# Load packages
import numpy as np
import matplotlib.pyplot as plt
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import iris.analysis.cartography



# Identifying ocean sediment and then calculate area
def land_to_ew(cube):
    """ Check if it is an ocean point and has a land point to the east or west
    
    """

    nt, nz, ny, nx = np.shape(cube.data)
    alldata = cube.data
   
    # Guess depth boundaries
    cube.coord('depth_1').guess_bounds()
    z_bounds = cube.coord('depth_1').bounds[:,:]
    # Guess the difference and fill array with difference
    z_diff = z_bounds[:,1]-z_bounds[:,0]

    print(nx,ny,nz)

    land_ew_data = np.zeros((nt, nz, ny, nx))
   
   
    for i in range(0,nx):
        im1 = i-1 # i minus 1
        ip1 = i+1 # i plus 1
        if im1 < 0:
            im1 = im1 + nx
        if ip1 >= nx:
            ip1 = ip1 - nx


        for j in range(0,ny):
            # Define arctic edge as in FAMOUS lat[0] is set to land even though ocean
            arctic_edge= ny-2
            for k in range(0,nz):
                dz = z_diff[k]
                if (alldata[0,k,j,i] < 1.0E10):
                       # ocean point
                    if ((alldata.mask[0,k,j,im1]) or
                        (alldata.mask[0,k,j,ip1])):
                        latsize = 2.5/360. * 2.0 * np.pi * 6371000. #calculating latsize 2(Pi)r
                       
                         
                        land_ew_data[0,k,j,i] = 1.0 * latsize  * dz
                        # Remove the arctic edge calculated as this is ocean not sediment 
                        land_ew_data[0,k,arctic_edge,i]=0.0
                        
                       
                        
                       
    land_ew_cube = cube.copy(data=land_ew_data)
    return land_ew_cube




def land_to_ns(cube):
    """ Check if it is an ocean point and has a land point to the north or south
    
    """

    nt, nz, ny, nx = np.shape(cube.data)
    alldata = cube.data
    latitude = cube.coord('latitude').points

    # Guess depth boundaries
    z_bounds = cube.coord('depth_1').bounds[:,:]
    # Guess the difference and fill array with difference
    z_diff = z_bounds[:,1]-z_bounds[:,0]
    
    print(nx,ny,nz)


    land_ns_data = np.zeros((nt, nz, ny, nx))
    
    # see if land point to north
    for j in range(1,ny):
        coslat = np.cos(latitude[j] * 2.0 * np.pi /360.) # Calculate cos(lat)
        for i in range(0,nx):
            for k in range(0,nz):
                dz = z_diff[k]
                if (alldata[0,k,j,i] < 1.0E10):
                   # ocean point
                    if alldata.mask[0,k,j-1,i]:
                        lonsize = (3.75/360. * 2.0 * np.pi * 6371000. * 
                                       coslat) #calculating lonsize
                        
                         
                        land_ns_data[0,k,j,i] = (1.0 * lonsize * dz)

    # see if land point to south
    for j in range(0,ny-1):
        coslat = np.cos(latitude[j] * 2.0 * np.pi /360.)
        # Define  arctic edge as in FAMOUS lat[0] is set to land even though ocean
        arctic_edge= ny-2
        for i in range(0,nx):
            for k in range(0,nz):
                dz = z_diff[k]
                if (alldata[0,k,j,i] < 1.0E10):
                   # ocean point
                    if alldata.mask[0,k,j+1,i]:
                        lonsize = (3.75/360. * 2.0 * np.pi * 6371000. * 
                                       coslat)
                        
                         
                        land_ns_data[0,k,j,i] = (land_ns_data[0,k,j,i] + 
                                                 (lonsize * dz))
                        # Remove the arctic edge calculated as this is ocean not sediment
                        land_ns_data[0,k,arctic_edge,i]=0.0
                        
                        
    land_ns_cube = cube.copy(data=land_ns_data)
    return land_ns_cube 



def land_to_oceanBottom(cube):
    """ Check if it is an ocean point and has a land point underneath it 
    
    """

    nt, nz, ny, nx = np.shape(cube.data)
    alldata = cube.data
    latitude = cube.coord('latitude').points
    
    print(nx,ny,nz)

    land_bot_data = np.zeros((nt, nz, ny, nx))
    
    # see if land point to the bottom
    for j in range(0,ny):
        coslat = np.cos(latitude[j] * 2.0 * np.pi /360.) # Calculate cos(lat)
        for i in range(0,nx):
            for k in range(0,nz-1):  
                if (alldata[0,k,j,i] < 1.0E10):
                    if (alldata.mask[0, k+1,j,i]):
                        latsize = (2.5/360. * 2.0 * np.pi * 6371000.)
                        lonsize = (3.75/360. * 2.0 * np.pi * 6371000. * coslat)
                   
                    
                        land_bot_data[0,k,j,i] = (1.0 * lonsize  * latsize)
                        
    for j in range(0,ny):
        coslat = np.cos(latitude[j] * 2.0 * np.pi /360.)
        for i in range(0,nx):
            for k in range(0,nz):  
                bot= nz-1
                if (alldata[0,bot,j,i] < 1.0E10):
                    latsize = (2.5/360. * 2.0 * np.pi * 6371000.)
                    lonsize = (3.75/360. * 2.0 * np.pi * 6371000. * coslat)
                    # Calculate sediment area of all ocean gridboxes at bottom
                    land_bot_data[0,bot,j,i] = (land_bot_data[0,bot,j,i] + (1.0*lonsize * latsize))
                    
                    
                    
                        
    land_bot_cube = cube.copy(data=land_bot_data)
    return land_bot_cube 





def sediment_depth(cube):
    """ Define the depth (m) at which boundary exchange will occur
    to create a mask depenedent on the depth specified
    
    Set to 6000 m if want to use whole ocean
    
    """
    nt, nz, ny, nx = np.shape(cube.data)

    
    depth_mask_data = np.zeros((nt, nz, ny, nx))
    depth= cube.coord('depth_1').points
    
    depthin= depth[:]
    
    depth_bounds = np.logical_and(depthin >=0.0, depthin <=3000.)
    
    depth_mask_data[0,depth_bounds,:,:]=1.0
    
    depth_mask_cube = cube.copy(data=depth_mask_data)
    
    return depth_mask_cube


#expt='xoiad'
fname = '/Users/suzie/Documents/PhD_2ndYear/Self_Isolation_Work/Git/BoundarySource/input_files/xoiado#pg000000003c1+.nc'

cube=iris.load_cube(fname,'sea_water_salinity')
#print(cube)


qplt.contourf(cube[0,0,:,:])
plt.show()

land_ew_cube = land_to_ew(cube)
#print(land_ew_cube.data)
qplt.contourf(land_ew_cube[0,0,:,:])
plt.show()

land_ns_cube = land_to_ns(cube)
#print(land_ns_cube.data)
qplt.contourf(land_ns_cube[0,0,:,:], levels =np.arange(0,400000,10000))
plt.show()

land_bot_cube = land_to_oceanBottom(cube)
#print(land_bot_cube.data)
qplt.contourf(land_bot_cube[0,:,:,0])
plt.show()


depth_mask_cube = sediment_depth(cube)
depth_mask_cube.rename('Boundary_Source_Mask')
#print(depth_mask_cube.data)
qplt.pcolormesh(depth_mask_cube[0,:,:,0])
plt.show()




# Create cube with total ocean sediment area in each grid box

# Add the area of each grid box
time, depth, lat, lon =cube.shape
sed_area_cubeData = np.zeros((time, depth, lat, lon))
sed_area_cubeData[0,:,:,:]= (land_ew_cube.data + land_ns_cube.data + land_bot_cube.data) * depth_mask_cube.data

sed_area_cube=cube.copy(data=sed_area_cubeData)
sed_area_cube.rename('SedimentArea_m2')

#print(sed_area_cube.shape)
#print(sed_area_cube)
#print(sed_area_cube.data)



# Calculate the total sediment area 
print ('Global ocean sediment area for boundary source calculation=', np.sum(sed_area_cube.data),'m-2')

# Save sediment area per gridbox as a netCDF file
iris.save(sed_area_cube, "/Users/suzie/Documents/PhD_2ndYear/Self_Isolation_Work/Git/BoundarySource/input_files/sedimentAreaPerGridbox.nc")








# PLOTTING!!

# Create a plot with 4 subplots for area ocean sediment 

fig = plt.figure(figsize=(10, 12))
plt.suptitle("Sediment area per gridbox FAMOUS $m^{2}$", fontsize=16)

cmap=plt.get_cmap('viridis')

# These are subplot grid parameters encoded as a single integer. 
#For example, "111" means "1x1 grid, first subplot" and "234" means "2x3 grid, 4th subplot".


ax1 = fig.add_subplot(321)
ax1.set_title('Surface')
iplt.pcolormesh(sed_area_cube[0,0,:,:], cmap=cmap)
plt.colorbar(ax=ax1, orientation="horizontal", extend='max')
ax1 = plt.gca()
#ax1.coastlines()

ax2 = fig.add_subplot(322)
ax2.set_title('447 m depth ')
iplt.pcolormesh(sed_area_cube[0,10,:,:], cmap=cmap)
plt.colorbar(ax=ax2, orientation="horizontal", extend='max')
ax2 = plt.gca()
#ax2.coastlines()

ax3 = fig.add_subplot(323)
ax3.set_title('4577 m  depth ')
iplt.pcolormesh(sed_area_cube[0,19,:,:],  cmap=cmap)
plt.colorbar(ax=ax3, orientation="horizontal", extend='max')
ax3 = plt.gca()
#ax3.coastlines()

ax4 = fig.add_subplot(324)
ax4.set_title('cross section')
iplt.pcolormesh(sed_area_cube[0,:,:,0], cmap=cmap , vmin=0.0, vmax=0.001)
plt.colorbar(ax=ax4, orientation="horizontal", extend='max')
ax4 = plt.gca()

ax5 = fig.add_subplot(325)
ax5.set_title('Boundary source mask at the surface')
iplt.pcolormesh(depth_mask_cube[0,0,:,:],  cmap=cmap)
plt.colorbar(ax=ax5, orientation="horizontal", extend='max')
ax5 = plt.gca()

ax6 = fig.add_subplot(326)
qplt.pcolormesh(depth_mask_cube[0,:,:,0])
ax6 = plt.gca()


plt.tight_layout()
plt.savefig('/Users/suzie/Documents/PhD_2ndYear/Self_Isolation_Work/Git/BoundarySource//plots/FAMOUS_AreaSediment.png', bbox_inches='tight',format='png', dpi=300)
plt.show()


