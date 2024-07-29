**************
GlobalDeltaChange
**************
.. image:: https://badge.fury.io/gh/jhnienhuis%2FGlobalDeltaChange.svg
    :target: https://badge.fury.io/gh/jhnienhuis%2FGlobalDeltaChange

.. image:: https://app.codacy.com/project/badge/Grade/0ae4939efdcd43b9b70e3ac605619f50
    :target: https://www.codacy.com/gh/jhnienhuis/GlobalDeltaChange/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jhnienhuis/GlobalDeltaChange&amp;utm_campaign=Badge_Grade
    
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
by Jaap Nienhuis, Utrecht University, 2022, version 3.0
(Version 3.0 includes better delta slopes)
by Jaap Nienhuis, Utrecht University, 2023, version 4.0
(Version 4.0 includes better fluvial sediment flux data, including bedload fluxes from Cohen et al 2022 and modern suspended load values from Dethier et al)

Use the data
#############

The data can be viewed interactively in `a arcGIS App <https://jhnienhuis.github.io/globaldeltas>`_.
Raw data is available here on github, formatted as `MATLAB .mat <https://github.com/jhnienhuis/GlobalDeltaChange/blob/master/GlobalDeltaData.mat>`_, `Shapefiles <https://github.com/jhnienhuis/GlobalDeltaChange/blob/master/export_data/GlobalDeltaMouth_shp.zip>`_, `NetCDF .nc <https://github.com/jhnienhuis/GlobalDeltaChange/blob/master/export_data/GlobalDeltaData.nc>`_, and `.kml <https://github.com/jhnienhuis/GlobalDeltaChange/blob/master/export_data/GlobalDeltaData.kml>`_ files. 

BasinArea: Drainage Basin Area (km2)
BasinID: HydroSheds ID of the drainage basin
BasinID2: HydroSheds ID of the drainage basin, with the last number added that specifies the continent
BasinID_ATLAS: ID of the drainage basin of the HydroATLAS product
channel_len: 30 point vector of each delta specifying the length to the next elevation (e.g., the 20th value represents to distance from the mouth to the 20 meter elevation contour)
channel_len_lat: 30 point vector of each delta specifying the latitude of the delta channel (in decimal degrees)
channel_len_lon: 30 point vector of each delta specifying the longitude of the delta channel (in decimal degrees)
channel_slope: overall delta surface slope (in m/m)
Continent: value from 1:8 specifying which continent the delta lies on
delta_name: the delta name, not complete
depth_mouth: channel depth at the river mouth (m)
depth_upstream: channel depth at the delta apex (m)
Discharge_dist: Modern river water discharge (m3/s)
Discharge_prist: Pristine (non human modified) river water discharge (m3/s)
Discharge_tide: Average tide-driven discharge at the river mouth (m3/s)
Hs: Average offshore significant wave height (m)
MouthLat: Latitude of the river mouth (decimal degrees)
MouthLon: Longitude of the river mouth (decimal degrees)
QRiver_bedload: Bedload sediment flux at the delta apex (kg/s)
QRiver_dist: Modern suspended load sediment flux at the delta apex (kg/s)
QRiver_prist: Pristine (non human modified) suspended load sediment flux at the delta apex (kg/s)
QTide: Average tide-driven sediment flux at the river mouth (kg/s)
QWave: Potential wave-driven sediment flux away from the river mouth (kg/s)
Region: Number indicating a certain coastal region
Region_str: Name of the Region
RiverID_ATLAS: ID of the river from the RiverATLAS product
shelf_depth: depth of the shelfbreak (m)
shelf_len: 31 point vector of each delta specifying the length to the next shelf contour line (e.g., the 20th value represents to distance from the mouth to the 20th value in shelf_lines)
shelf_len_lat: 31 point vector of each delta specifying the latitude of the delta steepest descent into the basin (in decimal degrees)
shelf_len_lon: 31 point vector of each delta specifying the longitude of the delta steepest descent into the basin (in decimal degrees)
shelf_lines: contour lines of the shelf
shelf_slope: average slope of the continental shelf (m/m)
shelf_width: distance to the shelf break (km)
TidalAmp: Average tidal amplitude (m)
Tp: Average wave period (s)
wave_lat: latitude of the wave data
wave_lon: longitude of the wave data
width_mouth: channel width at the river mouth (m)
width_upstream: channel width at the delta apex (m)


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
    Uses FAO data to find river names for deltas, where available. Needs updating.

Supplemental functions
**********

land_area_change/get_aquamonitor_data
    defines polygons for each river delta, and retrieves aquamonitor and earthsurfacewater explorer data to get delta coastal area land gain and loss within those regions. 
    These data are noisy, so use with caution and with appropriate estimates of data uncertainty. The GEE code can be found at:
    https://code.earthengine.google.com/21dd5f216c625b8696b4d9af6ee55215
    We manually define polygons for the 100 largest deltas (see GlobalDeltaMax100.kml), and use proxies for delta area size for the remaining deltas.
    
export_data/create_kml, create_netcdf, create_shapefile, create_shapefile_deltaland
    various functions to export relevant data to kml, netcdf, xlsx, and shapefile formats
    
misc/galloway_predictor
    function to plot output in the galloway triangle.

validation/global_delta_validation
    function to compare predictions against observations and put the resulting accuracy in the readme.rst file on github
    
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

Global Delta Accuracy
#############

The accuracy of the global delta dataset is assessed through comparison against field measurements and other datasets, scipts are validation data are in the subfolder "validation".

We compare the total number of predicted deltas (~11,000) against field observations of deltas that meet our definition (see the publication). We also compare the predicted morphology and give accuracy for individual predictions and for the global total. Lastly, we compare the delta land area change against a set of other datasets and observations.

For deltas on Madagascar, and additional deltas drawn at random from the dataset, we obtain the following confusion matrix:

+-----------+------------+------------+-----------+---------+
|           |              Observed                         |
+===========+============+============+===========+=========+
|           |            | Wave       | River     | Tide    |
+-----------+------------+------------+-----------+---------+
|           | Wave       | 241        |  012      | 033     |
+-----------+------------+------------+-----------+---------+
| Predicted | River      | 024        |  024      | 018     |
+-----------+------------+------------+-----------+---------+
|           | Tide       | 002        |  001      | 017     |
+-----------+------------+------------+-----------+---------+

For individual predictions, we retrieve the following accuracies

================    =======================
Morphology          Prediction accuracy (%)
----------------    -----------------------
Wave dominated               88%
River dominated              63%
Tide dominated               23%
================    =======================

Scaling up to the globe, we retrieve the following estimates for the global number of deltas and their morphologies

================    ==============  =======================
Morphology          Global number   Uncertainty (+/- 1std)
----------------    --------------  -----------------------
All deltas            10848             0371
Wave dominated        08234             0990 
River dominated       01840             0689
Tide dominated        00774             0598
================    ==============  =======================

The accuracy of our Aquamonitor-derived land area change estimats for global deltas is assessed by comparison against other models, and individual delta assessments.

================    ==============  =======================
Selection            Percentage of      Expressed in 
                     delta change       Area (km2/yr)  
----------------    --------------  -----------------------
Detection error         001%                001.00
Mapping error           153%                152.64
Intermodel error        092%                092.16
----------------    --------------  -----------------------
One delta (mean)        246%                245.80
All deltas (SE)         103%                103.16 
================    ==============  =======================









