**************
GlobalDeltaChange
**************

.. image:: https://github.com/DeltaRCM/pyDeltaRCM/actions/workflows/build.yml/badge.svg
    :target: https://github.com/DeltaRCM/pyDeltaRCM/actions

.. image:: https://badge.fury.io/gh/jhnienhuis%2FGlobalDeltaChange.svg
    :target: https://badge.fury.io/py/pyDeltaRCM

.. image:: https://codecov.io/gh/DeltaRCM/pyDeltaRCM/branch/develop/graph/badge.svg
  :target: https://codecov.io/gh/DeltaRCM/pyDeltaRCM

.. image:: https://app.codacy.com/project/badge/Grade/1c137d0227914741a9ba09f0b00a49a7
    :target: https://www.codacy.com/gh/DeltaRCM/pyDeltaRCM?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DeltaRCM/pyDeltaRCM&amp;utm_campaign=Badge_Grade

*GlobalDeltaChange* is a (1) theoretical framework to predict delta morphology and delta change, and (2) a set of codes to make this predictions on a global scale for ~11,000 deltas. Results and methods are described in `Nienhuis et al., 2020 <https://www.nature.com/articles/s41586-019-1905-9>`_

.. figure:: https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41586-019-1905-9/MediaObjects/41586_2019_1905_Fig1_HTML.png?as=webp
    
    Global delta morphology, as predicted by three sediment fluxes (Qwave, Qtide, and Qriver), within a ternary space and along Earths' coast.

Documentation
#############

Versioning
**********

by Jaap Nienhuis, Utrecht University, 2019, version 1.0
by Jaap Nienhuis, Utrecht University, 2021, version 2.0
(Version 2.0 includes the newest land/water change data from GSW, local wave estimates from local wind fetch, submarine and subaerial elevation, river names, and more.)

Use the data
#############

The data can be viewed interactively in `a GEE App <https://jhnienhuis.users.earthengine.app/view/globaldelta>`_.
Raw data is available here on github, formatted as `MATLAB .mat <https://github.com/jhnienhuis/GlobalDeltaChange/blob/master/GlobalDeltaData.mat>`_, `Shapefiles <https://github.com/jhnienhuis/GlobalDeltaChange/blob/master/export_data/GlobalDeltaMouth_shp.zip>`_, `NetCDF .nc <https://github.com/jhnienhuis/GlobalDeltaChange/blob/master/export_data/GlobalDeltaData.nc>`_, and `.kml <https://github.com/jhnienhuis/GlobalDeltaChange/blob/master/export_data/GlobalDeltaData.kml>`_ files. 

Reproduce the data
#############

To reproduce the GlobalDeltaData.mat file, run the following functions in this order: 

Main functions
**********
(1) find_river_mouth.m
    uses hydrosheds, DIVA, Durr, and SRTM to find all alluvial river mouths globally, furtheron referred to as deltas. Initiates the GlobalDeltaData.mat file

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

Supplemental functions
**********

land_area_change/get_aquamonitor_data
    defines polygons for each river delta, and retrieves aquamonitor and earthsurfacewater explorer data to get delta coastal area land gain and loss within those regions. 
    These data are noisy, so use with caution and with appropriate estimates of data uncertainty. The GEE code can be found at:
    https://code.earthengine.google.com/21dd5f216c625b8696b4d9af6ee55215
    We manually define polygons for the 100 largest deltas (see GlobalDeltaMax100.kml), and use proxies for delta area size for the remaining deltas.
    
export_data/create_kml, create_netcdf, create_shapefile, create_shapefile_deltaland
    various functions to export relevant data to kml, netcdf, and shapefile formats
    
misc/galloway_predictor
    function to plot output in the galloway triangle.
    
Input datasets
#############

Reproducing the data can be done with the following input datasets:

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

