
by Jaap Nienhuis, Utrecht University, 2019
by Jaap Nienhuis, Utrecht University, 2021, version 2.0
(Version 2.0 includes the newest land/water change data from GSW, local wave estimates from local wind fetch, submarine and subaerial elevation, river names, and more.)


This readme explains how to generate a global delta dataset including relevant river/wave/tidal parameters that can be used to predict delta morphology.
It also contains a lot of other useful functions to retrieve delta morphology, recent delta area change, and other delta attributes.

Input datasets:

- HydroSheds 15 arcsec drainage direction (DIR), flow accumulation (ACC), and basin outline (BAS) files
source: https://www.hydrosheds.org/

- DIVA typology_coastline
source: AT Vafeidis, G Boot, J Cox, R Maatens, L McFadden, RJ Nicholls, T Spencer, RSJ Tol, (2006) The DIVA database documentation, DINAS-COAST Consortium

- DURR dataset
source: Dürr, H.H., Laruelle, G.G., van Kempen, C.M. et al. Estuaries and Coasts (2011) 34: 441. https://doi.org/10.1007/s12237-011-9381-y

- NOAA vectorized shoreline
source: https://www.ngdc.noaa.gov/mgg/shorelines/

- WBMSed global discharge, pristine, and disturbed sediment fluxes
source: https://sdml.ua.edu/datasets-2/

- Global directional wave statistics (WaveWatch), and global tides (TOPEX)
source: https://jhnienhuis.users.earthengine.app/view/changing-shores

- SRTM, 1 arcsec (30 meter) resolution global topography
source: https://lpdaac.usgs.gov/products/srtmgl1v003/

- River Names, from FAO Aquamaps
source: http://www.fao.org/nr/water/aquamaps/

(note, I don't store these here because of versioning and file size limitations. Please get in touch if you can't find them, I will send them to you)

MAIN FUNCTIONS:
(1) find_river_mouth.m
    uses hydrosheds, DIVA, Durr, and SRTM to find all alluvial river mouths globally, furtheron referred to as deltas. Starts the GlobalDeltaData.mat file

(2) get_QRiver.m
    uses WBMSED to get a pristine and disturbed sediment and water flux to each delta. Optionally you can use get_QRiver_timeseries to get daily QRiver and Discharge output

(3) get_channel_slope.m
    uses SRTM and hydrosheds to extract river elevation profiles for all deltas up to 30 meters elevation
    
(4) get_bathy_profile.m
    uses etopo data to get steepest descent profiles of the underwater basin depths, from the river mouth to -100m
    
(5) get_Qwave.m
    adds wave data to each delta from WaveWatch. For deltas that are (partially) sheltered from wave approach angles, it estimates a fetch based on shoreline orientation.
    It uses the bretschneider fetch formula and WaveWatch wind data to estimate wave heights in sheltered locations. Uses get_global_fetch.m. 
    Optionally you can use get_QWave_timeseries to get daily wave statistics, or get_QWave_future to get estimates of future wave heights (up to 2100).

(6) get_Qtide.m
    adds tide data to each delta, based on TOPEX data
    
(7) get_hydrobasins_id.m
    adds identifiers from the new WWF HydroATLAS, HydroBasins, and HydroRIVERS datasets

(8) add_names_to_deltas.m
    Uses FAO data to find river names for deltas, where available.

SUPPLEMENTAL FUNCTIONS

land_area_change/get_aquamonitor_data
    defines polygons for each river delta, and retrieves aquamonitor and earthsurfacewater explorer data to get delta coastal area land gain and loss within those regions. 
    These data are noisy, so use with caution and with appropriate estimates of data uncertainty. The GEE code can be found at:
    https://code.earthengine.google.com/21dd5f216c625b8696b4d9af6ee55215
    We manually define polygons for the 100 largest deltas (see GlobalDeltaMax100.kml), and use proxies for delta area size for the remaining deltas.
    
export_data/create_kml, create_netcdf, create_shapefile, create_shapefile_deltaland
    various functions to export relevant data to kml, netcdf, and shapefile formats
    
misc/galloway_predictor
    function to plot output in the galloway triangle.
