(c) Jaap Nienhuis, Utrecht University, 2019

to generate the GlobalDeltaData.mat file:
required datasets:
- HydroSheds 15s acc_15s.bil files for every continent
- DIVA typology_coastline.shp (found on this website)
- DURR dataset
- WBMSed
- TPXO
- WaveDirGlo
- SRTM


(1) run find_river_mouth.m

(2) run get_QRiver.m

(3) run combine_continent_files.m

(4) run get_channel_slope.m

(5) 