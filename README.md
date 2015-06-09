# r_movement_homerange
collection of animal movement and utilization distribution scripts in R

site_fidelity.r

This script creates 2 types of utilization distribution estimations from point locations, minimum convex polygon (also know as convex hull) and kernel density estimation. It uses mainly functions from the adehabitatHR package, with some customization to allow for rescaling of locations prior to kernel density estimation, replicating the option in the Kernel Density Estimation tool in the Home Range Tools extension developed for ArcView and ArcGIS 9 (http://flash.lakeheadu.ca/~arodgers/hre/).

The script requires a specific input file containing animal ID, track ID, date, time, and X/Y coordinates. It outputs a summary of the run and shapefiles for each utilization distribution.



mcp_and_kde.r

This script creates 2 types of utilization distribution estimations from point locations, minimum convex polygon (also know as convex hull) and kernel density estimation. It uses mainly functions from the adehabitatHR package, with some customization to allow for rescaling of locations prior to kernel density estimation, replicating the option in the Kernel Density Estimation tool in the Home Range Tools extension developed for ArcView and ArcGIS 9 (http://flash.lakeheadu.ca/~arodgers/hre/).

The script requires a specific input file containing animal ID, track ID, date, time, and X/Y coordinates. It outputs a summary of the run and shapefiles for each utilization distribution.

