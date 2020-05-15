Contains the eNd signature of the sediment for the boundary source
Detail all files put into this directory here so clear which one to pick to create boundary conditions


1. end_continentalExtrapAllOcean.nc
this file has the continental eNd map extrapolated across the whole ocean using NN horizontally 
Extrapolated so eNd is constant with depth
This is to create boundary conditions for a test boundary source run with depth set to the whole ocean. 
Once new seafloor eNd map has been created this map file will be updated. 
The depth to which boundary source occurs is the whole ocean

2. continentalMargin_eNd.nc
this file has the continental eNd map extrapolated to the continental margins
With continental margins defined as where sediment thickness is 
>= 2000 m 
The ocean depth to which boundary source occurs is restricted to <= 3000 m