#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 14:48:17 2020

Create boundary conditions [144Nd] and [143Nd] for the boundary source
Units: kg(Nd)/m3/s(*1E18)
Ready to create ancil files to input into FAMOUS/HadCM3


Loads two files
1. eNd value for the boundary source sediment (eNd_continentalExtrapDepth.nc)
2. [Nd] total Nd from boundary source (BoundarySource_Nd.nc)

Note:
Check the depth of the sediment specified
Chck eNd map, this will vary with depth, updated following seafloor eNd map ect


Converting [Nd] and eNd into [144Nd] and [143Nd] based on the equations:
    
eNd = ((IR/CHUR)-1)*10^4 (where CHUR = 0.512638)
1. isotopic_ratio = ((eNd/10^4)+1) * CHUR 
2. [Nd144]= [Nd]/(IR+1)
3. [Nd143] = [Nd]/((1/IR)+1))


@author: Suzanne Robinson
"""

# Load packages
import numpy as np
import matplotlib.pyplot as plt
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import iris.analysis.cartography
import cmocean


# LOAD FILES 

# eNd
#fname1 = '/Users/suzie/Documents/PhD_2ndYear/Self_Isolation_Work/Git/BoundarySource/input_files/eNd/eNd_continentalExtrapAllOcean.nc'
fname1 = '/Users/suzie/Documents/PhD_2ndYear/Self_Isolation_Work/Git/BoundarySource/input_files/eNd/continentalMargin_eNd.nc'

# [Nd]
fname2 = '/Users/suzie/Documents/PhD_2ndYear/Self_Isolation_Work/Git/BoundarySource/input_files/boundarySource_Nd.nc'

cube=iris.load_cube(fname1)
cube2=iris.load(fname2)

print(cube)
print(cube2)

eNd_cube = cube[0] # Sediment area
Nd_cube = cube2[0] # Volume gridbox



""" 1.  Calculating the IR using eNd
    IR = ((eNd/10^4)+1)*CHUR (where IR = the isotopic ratio)
    
"""

isotopic_ratio = (((eNd_cube/10000)+1)*0.512638)
isotopic_ratio.rename('isotopic_ratio_Nd')

# Plot the isotopic ratio (should be around the value of CHUR) 
qplt.pcolormesh(isotopic_ratio[0,:,:])
plt.show()



""" 2. Calculating [144Nd]
     [Nd144]= [Nd]/(IR+1)

"""

Nd144_cube = (Nd_cube/(isotopic_ratio+1.0))
Nd144_cube.rename('Nd144_kg_m3_s_1E18')
Nd144_cube.coord("depth_1").bounds = None


# Plot [144Nd]
qplt.pcolormesh(Nd144_cube[0,0,:,:])
plt.show()

qplt.pcolormesh(Nd144_cube[0,10,:,:])
plt.show()

qplt.pcolormesh(Nd144_cube[0,18,:,:])
plt.show()


""" 3. Calculating [143Nd]
    [Nd143] = [Nd]/((1/IR)+1))
        
"""

# [Nd]/((IR^-1)+1)
Nd143_cube = Nd_cube/((isotopic_ratio**-1.0)+1.0)
Nd143_cube.rename('Nd143_kg_m3_s_1E18')
Nd144_cube.coord("depth_1").bounds = None



# Plot [143Nd]
qplt.pcolormesh(Nd143_cube[0,0,:,:])
plt.show()

qplt.pcolormesh(Nd143_cube[0,2,:,:])
plt.show()

qplt.pcolormesh(Nd143_cube[0,18,:,:])
plt.show()

qplt.pcolormesh(Nd143_cube[0,:,:,0], vmin=0.00, vmax=0.001)
plt.show()


print(Nd144_cube.shape)
print(Nd143_cube.shape)



""" Save [Nd144] and [Nd143] cubes as seperate .nc files
units are kg(Nd)/m3/s (*1E18)

"""

# Where there are gaps in the ocean use fill value as 0.0 rather than -999.999
# then when creating ancil file do not add a land sea mask
# Use -999.999 when complete land mask and add land sea mask when creating ancil files 

# Save [144Nd] as a netCDF file
iris.save(Nd144_cube, "/Users/suzie/Documents/PhD_2ndYear/Self_Isolation_Work/Git/BoundarySource/boundary_conditions/Nd144_boundarySource_3000m.nc", fill_value=0.0, netcdf_format='NETCDF3_CLASSIC')
# Save [143Nd] as a netCDF file as a netCDF file
iris.save(Nd143_cube, "/Users/suzie/Documents/PhD_2ndYear/Self_Isolation_Work/Git/BoundarySource/boundary_conditions/Nd143_boundarySource_3000m.nc", fill_value=0.0, netcdf_format='NETCDF3_CLASSIC')



    





"""
Plotting the boundary conditions for [144Nd] and [143Nd] for the boundary source

"""

 # Create a plot with 6 subplots for Nd144
cmap = cmocean.cm.tempo
fig = plt.figure(figsize=(10, 12))
plt.suptitle("Bounday conditions for $^{144}Nd$ from boundary exchange kg/$m^{3}$/s (*1E18)", fontsize=16)



# These are subplot grid parameters encoded as a single integer. 
#For example, "111" means "1x1 grid, first subplot" and "234" means "2x3 grid, 4th subplot".
ax1 = fig.add_subplot(321)
ax1.set_title('Surface')
iplt.pcolormesh(Nd144_cube[0,0,:,:], vmin=0.0, vmax=0.001, cmap=cmap)
plt.colorbar(ax=ax1, orientation="horizontal", extend='max')
ax1 = plt.gca()
ax1.coastlines()

ax2 = fig.add_subplot(322)
ax2.set_title('67 m depth ')
iplt.pcolormesh(Nd144_cube[0,5,:,:], vmin=0.0, vmax=0.6, cmap=cmap)
plt.colorbar(ax=ax2, orientation="horizontal", extend='max')
ax2 = plt.gca()
ax2.coastlines()

ax3 = fig.add_subplot(323)
ax3.set_title('447 m depth ')
iplt.pcolormesh(Nd144_cube[0,10,:,:], vmin=0.0, vmax=0.6, cmap=cmap)
plt.colorbar(ax=ax3, orientation="horizontal", extend='max')
ax3 = plt.gca()
ax3.coastlines()

ax4 = fig.add_subplot(324)
ax4.set_title('4577 m depth')
iplt.pcolormesh(Nd144_cube[0,18,:,:], vmin=0.0, vmax=0.6, cmap=cmap)
plt.colorbar(ax=ax4, orientation="horizontal", extend='max')
ax4 = plt.gca()
ax4.coastlines()

ax5 = fig.add_subplot(325)
ax5.set_title('5192.65 m depth')
iplt.pcolormesh(Nd144_cube[0,19,:,:], vmin=0.0, vmax=0.6, cmap=cmap)
plt.colorbar(ax=ax5, orientation="horizontal", extend='max')
ax5 = plt.gca()
ax5.coastlines()

ax6 = fig.add_subplot(326)
ax6.set_title(' cross section')
iplt.pcolormesh(Nd144_cube[0,:,:,0], vmin=0.0, vmax=0.003, cmap=cmap)
plt.colorbar(ax=ax6, orientation="horizontal", extend='max')


plt.tight_layout()
plt.savefig('/Users/suzie/Documents/PhD_2ndYear/Self_Isolation_Work/Git/BoundarySource/plots/Nd144BoundarySourceBoundaryConditions.png', bbox_inches='tight',format='png', dpi=300)
plt.show()





 # Create a plot with 6 subplots for Nd143
fig = plt.figure(figsize=(10, 12))
plt.suptitle("Bounday conditions for $^{143}Nd$ from boundary exchange kg/$m^{3}$/s (*1E18)", fontsize=16)



# These are subplot grid parameters encoded as a single integer. 
#For example, "111" means "1x1 grid, first subplot" and "234" means "2x3 grid, 4th subplot".
ax1 = fig.add_subplot(321)
ax1.set_title('Surface')
iplt.pcolormesh(Nd143_cube[0,0,:,:], vmin=0.0, vmax=0.001, cmap=cmap)
plt.colorbar(ax=ax1, orientation="horizontal", extend='max')
ax1 = plt.gca()
ax1.coastlines()

ax2 = fig.add_subplot(322)
ax2.set_title('67 m depth ')
iplt.pcolormesh(Nd143_cube[0,5,:,:], vmin=0.0, vmax=0.6, cmap=cmap)
plt.colorbar(ax=ax2, orientation="horizontal", extend='max')
ax2 = plt.gca()
ax2.coastlines()

ax3 = fig.add_subplot(323)
ax3.set_title('447 m depth ')
iplt.pcolormesh(Nd143_cube[0,10,:,:], vmin=0.0, vmax=0.6, cmap=cmap)
plt.colorbar(ax=ax3, orientation="horizontal", extend='max')
ax3 = plt.gca()
ax3.coastlines()

ax4 = fig.add_subplot(324)
ax4.set_title('4577 m depth')
iplt.pcolormesh(Nd143_cube[0,18,:,:], vmin=0.0, vmax=0.6, cmap=cmap)
plt.colorbar(ax=ax4, orientation="horizontal", extend='max')
ax4 = plt.gca()
ax4.coastlines()

ax5 = fig.add_subplot(325)
ax5.set_title('5192.65 m depth')
iplt.pcolormesh(Nd143_cube[0,19,:,:], vmin=0.0, vmax=0.6, cmap=cmap)
plt.colorbar(ax=ax5, orientation="horizontal", extend='max')
ax5 = plt.gca()
ax5.coastlines()

ax6 = fig.add_subplot(326)
ax6.set_title(' cross section')
iplt.pcolormesh(Nd143_cube[0,:,:,0], vmin=0.0, vmax=0.003, cmap=cmap)
plt.colorbar(ax=ax6, orientation="horizontal", extend='max')



plt.tight_layout()
plt.savefig('/Users/suzie/Documents/PhD_2ndYear/Self_Isolation_Work/Git/BoundarySource/plots/Nd143BoundarySourceBoundaryConditions.png', bbox_inches='tight',format='png', dpi=300)
plt.show()


