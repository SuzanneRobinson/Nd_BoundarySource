This directory contains the scripts and input files to produce 
the boundary conditions for the boundary source input

Creates 2 .nc files in boundary_conditions 
Boundary source of [144Nd] and [143Nd] in kg(Nd)/m3/s *1E18


Scripts
1. area_ocean_sediment.py
calculates the area of sediment in contact with the ocean for the boundary source
need to define the depth of boundary source in m 
outputs total sediment area (m2) (the sediment area where boundary source occurs)

2. calculating_[Nd].py
calculates the total [Nd] for the boundary source
kg(Nd)/m3/s *1E18
Requires:
- total_sediment_area (the sediment area where boundary source occurs)
- f(bs) (tuning parameter, currently following Rempfer et al., 2011: 5.5*10^9 g(Nd)/yr)

3. calculating_[144Nd]_[143Nd].py
creates the boundary conditions for each Nd isotope- outputted to boundary_conditions
Requires:
- eNd input file for the eNd signature of the sediment


input files
1. oceanVolumePerGridbox.nc 
volume of each ocean grid box (m3)

2. xoiado#pg0000000003c1.nc
sea_water_salinity used to create land-sea mask and
calculate sediment area

3. sedimentAreaPerGridbox.nc
calculated from area_ocean_sediment.py

4. boundarySource_Nd.nc
calculated from calculating_[Nd]_boundary_source.py


input_files/eNd
These are the files representing the eNd signature of the ocean sediment 
1. eNd_continentalExtrapAllOcean.nc
eNd from continental map extrapolated using NN across whole ocean
and extrapolated at depth so constant with depth

2. continentalMargin_eNd.nc
eNd from continental map extrapolated to the continental margins
margins defined where sediment thickness >= 2000 m 