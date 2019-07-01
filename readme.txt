(c) Jaap Nienhuis, Utrecht University, 2019

This readme explains how to generate a global delta dataset including relevant river/wave/tidal parameters that can be used to predict delta morphology.

Required datasets:

- HydroSheds 15 arcsec drainage direction (DIR), flow accumulation (ACC), and basin outline (BAS) files
source: https://www.hydrosheds.org/

- DIVA typology_coastline.shp
source: AT Vafeidis, G Boot, J Cox, R Maatens, L McFadden, RJ Nicholls, T Spencer, RSJ Tol, (2006) The DIVA database documentation, DINAS-COAST Consortium

- DURR dataset
source: Dürr, H.H., Laruelle, G.G., van Kempen, C.M. et al. Estuaries and Coasts (2011) 34: 441. https://doi.org/10.1007/s12237-011-9381-y

- WBMSed global discharge, pristine, and disturbed sediment fluxes
source: https://sdml.ua.edu/datasets-2/

- Global directional wave statistics (WaveWatch), and global tides (TOPEX)
source: https://jhnienhuis.users.earthengine.app/view/changing-shores

- SRTM, 1 arcsec (30 meter) resolution global topography
source: https://lpdaac.usgs.gov/products/srtmgl1v003/


(1) run find_river_mouth.m

(2) run get_QRiver.m

(3) run get_channel_slope.m

(4) run combine_continent_files.m

(5) run get_Qwave

(6) run get_Qtide

(7a) run get_aquamonitor_data

(7b) run https://code.earthengine.google.com/451465dc579747ae0b39b5f2f8ad1b12

(7c) run get_aquamonitor_data

(8) run create_netcdf
