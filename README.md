# Memory_Movement

## Codes for How to fit memory-models 
### Sandhill_Crane 
- 0.calculate.UDs.Rmd
  - use ctmm package to fit the tracking data to a continuous-time movement model 
  - calculate a Kriged occurrence distribution (OD) estimate per each year model 
  - convert each year ODs as a raster layer and save it for further analysis 

- 1.crane_trk.Rmd:
  - extract land cover layer (NLCD2016) and each year ODs raster layers based on tracking data per each year 
  - convert the tracking data to amt object for fitting ssf 

- 2.crane_memory.Rmd: 
  - call each year tracking data (amt object) and fit SSFs per each yr data including the effects of landcovers + ODs for hatched year + ODs for the previous years 

### Mechanistic_approach
