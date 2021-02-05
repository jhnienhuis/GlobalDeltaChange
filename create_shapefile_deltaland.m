function create_shapefile_deltaland
%DeltaEarthEngineShapefile
load([dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep 'GlobalDeltaData.mat'],'Discharge_prist','QRiver_prist','MouthLon','MouthLat','channel_len_lat','channel_len_lon','BasinID2','channel_len');
res=1;
remfun = @(lon) (rem(res*360-1+lon,res*360)+1);

%create buffer per delta
DeltaRadius = max(2,sqrt(1.07.*Discharge_prist.^0.7.*QRiver_prist.^0.45/pi))./111; %.*(cos(deg2rad(MouthLat)));

%load shoreline
shore = shaperead('D:\OneDrive - Universiteit Utrecht\WorldCoastline\osm_coastline\coastline_z8.shp');
idx = zeros(length(shore),1);
island = true(length(shore),1);
for ii=1:length(shore)
   idx(ii) = length(shore(ii).X);
   bx = ceil(shore(ii).X*100);
   by = ceil(shore(ii).Y*100);
   island(ii) = (bx(1)==bx(end-1) && by(1)==by(end-1));
   if max(abs(diff(bx)))>10000;
       shore(ii).X = rad2deg(unwrap(deg2rad(shore(ii).X)));
   end
   
end
1
shore2 = shore(~island | idx>200);
shore2 = struct2cell(shore2);

%put lat and lon in imaginary numbers to make it a 1 line search
shore_lon = remfun(cell2mat(shore2(3,:)));
shore_lat = cell2mat(shore2(4,:));
shore_both_all = shore_lon + 1i*shore_lat;
[shore_both,~,shoreline_idx] = unique(round(0.5*shore_both_all));
shore_both = shore_both*2;

%put delta lat and lon also in im
delta_both = (remfun(MouthLon)) + 1i*(MouthLat);
2
%find minimum distance
idx = zeros(size(delta_both));
min_dist = zeros(size(delta_both));
for ii=1:length(delta_both),
   [min_dist(ii),idx(ii)] = min(abs(shore_both-delta_both(ii)));
end

%find OSM shoreline per delta within a search radius depending on the buffer
delta_shoreline = cell(size(MouthLon));
delta_shoreline_lat = cell(size(MouthLon));
delta_shoreline_lon = cell(size(MouthLon));
for ii=1:length(delta_both),
    x = find(shoreline_idx==idx(ii));
    [~,idxx] = min(abs(shore_both_all(x)-delta_both(ii)));
    delta_shoreline{ii} = x(abs(shore_both_all(x)-shore_both_all(x(idxx)))<DeltaRadius(ii));
    if length(delta_shoreline{ii})==1,
        delta_shoreline{ii}(2) = delta_shoreline{ii}(1)+1; 
    end
    delta_shoreline_lat{ii} = shore_lat(delta_shoreline{ii});
    delta_shoreline_lon{ii} = shore_lon(delta_shoreline{ii});
end
3

delta_land_lat = cell(size(MouthLon));
delta_land_lon = cell(size(MouthLon));

%add inland river profile to get entire river delta
for ii=1:length(MouthLon),
    delta_land_lat{ii} = [delta_shoreline_lat{ii} MouthLat(ii) channel_len_lat(ii,:) ]; 
    delta_land_lon{ii} = [delta_shoreline_lon{ii} MouthLon(ii) channel_len_lon(ii,:) ]; 
end

DeltaRadius = sqrt(1.07.*Discharge_prist.^0.7.*QRiver_prist.^0.45/pi/100);
DeltaBuffer = 3000+(DeltaRadius*2000);

delta_polygons = mapshape(delta_land_lon,delta_land_lat,'BasinID2',double(BasinID2),'Buffer',(DeltaBuffer));

shapewrite(delta_polygons,'GlobalDeltaLand.shp')
zip('GlobalDeltaLand.zip',{'GlobalDeltaLand.shp','GlobalDeltaLand.shx','GlobalDeltaLand.dbf'})
delete('GlobalDeltaLand.shp','GlobalDeltaLand.shx','GlobalDeltaLand.dbf')
end