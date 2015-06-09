# r_movement_homerange

#### collection of animal movement and utilization distribution scripts in R

**site_fidelity.r**

This script compares animal movement trajectories to a set of randomly generated trajectories (generated by randomizing the true trajectories' steps' angles and distances. It compares the true trajectories' R-squared and linearity value to the distribution of R-squared and linearity values of the random replicates, following a procedure implemented in the Site Fidelity tool of Animal Movement Extension for ArcView (http://alaska.usgs.gov/science/biology/spatial/gistools/animalmovementdoca.pdf).

It takes inputs of animal locations and a shapefile with one polygon which represents the area available to the animal (this constrains the random walks). Plots are output displaying the trajectory and random walks, and a output table with r-squared, linearity, and p-values for each track.



**mcp_and_kde.r**

This script creates 2 types of utilization distribution estimations from point locations, minimum convex polygon (also know as convex hull) and kernel density estimation. It uses mainly functions from the adehabitatHR package, with some customization to allow for rescaling of locations prior to kernel density estimation, replicating the option in the Kernel Density Estimation tool in the Home Range Tools extension developed for ArcView and ArcGIS 9 (http://flash.lakeheadu.ca/~arodgers/hre/).

The script requires a specific input file containing animal ID, track ID, date, time, and X/Y coordinates. It outputs a summary of the run and shapefiles for each utilization distribution.

