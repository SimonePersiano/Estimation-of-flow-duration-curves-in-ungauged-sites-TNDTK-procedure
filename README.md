https://zenodo.org/badge/latestdoi/366670107

This repository contains an example of application of the following procedure:
1. extraction of the POR-FDC (period-of-record flow-duration-curve) from the daily streamflow series observed at the given gauged sites;
2. computation of the total negative deviation (TND, as defined in Pugliese et al., 2014, 2016);
3. application of total negative deviation top-kriging (TNDTK; see e.g. Pugliese et al., 2014, 2016) for computing POR-FDCs at the given (ungauged) target sites.

The required inputs are:
* daily streamflow series observed for different river cross-sections;
* shapefile of the catchment boundaries associated with the gauged river cross-sections;
* shapefile of the catchment boundaries associated with the target ungauged river cross-sections.

The folder contains: 
* FDC_TND_TNDTK.R: main R script performing the operations described above;
* tnd.R: function for computing the TND;
* resample_FDC.R: funtion for resampling FDC (for the durations of interest);
* GaugedCatchments: .zip folder with the shapefile of the gauged catchment boundaries;
* GaugedStations: .zip folder with the shapefile of the gauged stations;
* UngaugedCatchments: .zip folder with the shapefile of the ungauged (target) catchment boundaries;
* StreamflowData: .zip folder with the .csv files of the daily streamflow series observed at the given gauged sites.
(Please do remember to unzip the appropriate folders before running the code)

The example application considers 27 gauged catchments and 3 target ungauged catchments located in Tyrol (Austria) and South Tyrol (Italy), within an area considered for the activities of the European INTERREG Italy-Austria project named ALFFA (Holistic multiscale AnaLysis of the factors and their effect on the Fish Fauna in inner‐Alpine space).
