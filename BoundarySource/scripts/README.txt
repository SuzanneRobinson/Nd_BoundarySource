Script order

1. area_ocean_sediment 
Need to define depth of boundary source here
Obtain total sediment area (m2) to be used in next script

2. calculating_[Nd].py
Input 
- total sediment area (m2)
- f(bs) - tuning parameter g(Nd)/yr

3. calculating_[144Nd]_[143Nd].py
Input
- eNd of the sediment (eg. extrapolated margin eNd map, seafloor eNd map)
Output
- /boundary_conditions/Nd144_boundarySource.nc
- /boundary_conditions/Nd143_boundarySource.nc