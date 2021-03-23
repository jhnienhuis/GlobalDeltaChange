%% DeltaEarthEngineShapefile
load([dropbox filesep 'WorldDeltas' filesep 'scripts' filesep 'GlobalDeltaData.mat']);
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
%make shape file
%improve BasinID(!) 
BasinID = (10*BasinID+Continent);

DeltaRadius = sqrt(1.07.*Discharge_prist.^0.7.*QRiver_prist.^0.45/pi/100);
DeltaBuffer = 1000+(DeltaRadius*1000);

gee_shapefile = mapshape(delta_shoreline_lon,delta_shoreline_lat,'BasinID',(BasinID),'Buffer',(DeltaBuffer));

shapewrite(gee_shapefile,'GlobalDeltaShorelineData.shp')

%% Create shapefile for manual data for 100 largest deltas
load([dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep 'GlobalDeltaData.mat'],'MouthLon','MouthLat','BasinID2','QRiver_prist');
%first write kml file
[~,idx] = sort(QRiver_prist,'descend');
kname = 'GlobalDeltaMax100';

%kmlwritepoint([kname '.kml'],MouthLat(idx(1:100)),MouthLon(idx(1:100)),'name',string(BasinID2(idx(1:100))))

%retrieve from google earth
p = kml2struct([kname '_poly.kml']);

%match with earlier struct
MouthLon(MouthLon>180) = MouthLon(MouthLon>180)-360;
[~,idxx] = min((MouthLat(idx(1:100))-arrayfun(@(x) nanmean(x.Lat),p)).^2+(MouthLon(idx(1:100))-arrayfun(@(x) nanmean(x.Lon),p)).^2,[],2);
%scatter(MouthLon(idx(1:100)),MouthLat(idx(1:100))), hold on,
%scatter(arrayfun(@(x) nanmean(x.Lon),p),arrayfun(@(x) nanmean(x.Lat),p))
for ii=1:100,
    delta_shoreline_lat{ii} = p(idxx(ii)).Lat;
    delta_shoreline_lon{ii} = p(idxx(ii)).Lon;
end

gee_shapefile100 = mapshape(delta_shoreline_lon,delta_shoreline_lat,'BasinID2',double(BasinID2(idx(1:100))));    
gee_shapefile100.Geometry = 'polygon';

shapewrite(gee_shapefile100,[kname '.shp'])
zip([kname '.zip'],{[kname '.shp'],[kname '.shx'],[kname '.dbf']})
delete([kname '.shp'],[kname '.shx'],[kname '.dbf'])

%% put results back into mat file
clr
f = [dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep];
ee = load([f 'GlobalDeltaData.mat'],'BasinID2','delta_name');

%file exported from earth engine
fileID = fopen([f 'land_area_change' filesep 'GlobalDeltaChange.csv'],'r');
data = textscan(fileID, '%q%f%f%f%f%f%f%f%f%q%[^\n\r]', 'Delimiter', ',', 'HeaderLines' ,1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

pekel2_dry = string(data{10});
pekel2_dry = cell2mat(arrayfun(@str2num,pekel2_dry,'UniformOutput',false));
pekel2_wet = string(data{11});
pekel2_wet = cell2mat(arrayfun(@(x) (str2num(str2num(x))),pekel2_wet,'UniformOutput',false));


data = cell2mat(data(2:9)); %basin ID, deltaArea, change, deposition, erosion
data_aqua = data(:,[4 6 8])/28; %aquamonitor change deposition erosion per year
data_pekel = data(:,[3 5 7])/31; data_pekel(:,3) = data_pekel(:,3).*-1;

[~,idx] = ismember(ee.BasinID2,data(:,1));
idx(idx==0) = 2;

%put in structure
ee.net_pekel = data_pekel(idx,1);
ee.dep_pekel = data_pekel(idx,2);
ee.ero_pekel = data_pekel(idx,3);

ee.net_aqua = data_aqua(idx,1);
ee.dep_aqua = data_aqua(idx,2);
ee.ero_aqua = data_aqua(idx,3);


%convert pekel2 to rates
pekel2_dry = pekel2_dry(idx,:);
pekel2_wet = pekel2_wet(idx,:);

la_max = max(pekel2_dry+pekel2_wet,[],2);
fr_dry = pekel2_dry./(pekel2_dry+pekel2_wet);
null = max(0,la_max-pekel2_dry-pekel2_wet);
dry_corr = (pekel2_dry+fr_dry.*null);
dry_corr = dry_corr-nanmean(dry_corr,2);

t = datenum(1984:2019,1,1);
[ee.net_pekel2_t,~,~] = datevec(t);
ee.net_pekel2_y = zeros(size(dry_corr,1),36);
ee.net_pekel2 = zeros(size(dry_corr,1),1);
ee.net_pekel2_unc = zeros(size(dry_corr,1),2);

fitType = fittype('poly1');

for ii=1:size(dry_corr,1),
    idxnan = ~isnan(dry_corr(ii,:));
    %[rate(ii,[1 2]),S] = polyfit(tt(~idx),dry_corr(ii,~idx),1);
    %[~,rate_unc(ii)] = polyval(rate(ii,:),0,S);
    if sum(idxnan)<2, 
        ee.net_pekel2(ii)=0; 
        ee.net_pekel2_unc(ii,:)=0; 
        ee.net_pekel2_y(ii,1:36) = 0; 
        continue, 
    end
    
        
    fi = fit(t(idxnan)'./365,dry_corr(ii,idxnan)',fitType);
    ee.net_pekel2(ii) = fi.p1;
    unc = confint(fi);
    ee.net_pekel2_unc(ii,:) = unc(:,1);

    ee.net_pekel2_y(ii,:) = dry_corr(ii,:);
    %accumarray(ee.net_pekel2_t'-1984,dry_corr(ii,:)',[],@nanmean);
    
end

%add to .mat file
save([f 'land_area_change' filesep 'GlobalDeltaData_AreaChange.mat'],'-struct','ee');


%% put 100 largest deltas also in mat file

clr
f = [dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep];
ee = load([f 'land_area_change' filesep 'GlobalDeltaData_AreaChange.mat']);

%file exported from earth engine
fileID = fopen([f 'land_area_change' filesep 'GlobalDeltaMax100.csv'],'r');
data = textscan(fileID, '%q%f%f%f%f%f%f%f%q%q%[^\n\r]', 'Delimiter', ',', 'HeaderLines' ,1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

pekel2_dry = string(data{9});
pekel2_dry = cell2mat(arrayfun(@str2num,pekel2_dry,'UniformOutput',false));
pekel2_wet = string(data{10});
pekel2_wet = cell2mat(arrayfun(@str2num,pekel2_wet,'UniformOutput',false));


data = cell2mat(data(2:8)); %basin ID, change, deposition, erosion
data_aqua = data(:,[3 5 7])/28; %aquamonitor change deposition erosion per year
data_pekel = data(:,[2 4 6])/31; data_pekel(:,3) = data_pekel(:,3).*-1;

[~,idx] = ismember(data(:,1),ee.BasinID2);

%put in structure
ee.net_pekel(idx) = data_pekel(:,1);
ee.dep_pekel(idx) = data_pekel(:,2);
ee.ero_pekel(idx) = data_pekel(:,3);

ee.net_aqua(idx) = data_aqua(:,1);
ee.dep_aqua(idx) = data_aqua(:,2);
ee.ero_aqua(idx) = data_aqua(:,3);


%convert pekel2 to rates
%pekel2_dry = pekel2_dry(idx,:);
%pekel2_wet = pekel2_wet(idx,:);

la_max = max(pekel2_dry+pekel2_wet,[],2);
fr_dry = pekel2_dry./(pekel2_dry+pekel2_wet);
null = max(0,la_max-pekel2_dry-pekel2_wet);
dry_corr = (pekel2_dry+fr_dry.*null);
dry_corr = dry_corr-nanmean(dry_corr,2);

t = datenum(1984:2019,1,1);
[ee.net_pekel2_t,~,~] = datevec(t);
ee.net_pekel2_y(idx,:) = zeros(size(dry_corr,1),36);
ee.net_pekel2(idx,:) = zeros(size(dry_corr,1),1);
ee.net_pekel2_unc(idx,:) = zeros(size(dry_corr,1),2);

fitType = fittype('poly1');

for ii=1:size(dry_corr,1),
    idxnan = ~isnan(dry_corr(ii,:));
    %[rate(ii,[1 2]),S] = polyfit(tt(~idx),dry_corr(ii,~idx),1);
    %[~,rate_unc(ii)] = polyval(rate(ii,:),0,S);
    if sum(idxnan)<2, 
        ee.net_pekel2(idx(ii))=0; 
        ee.net_pekel2_unc(idx(ii),:)=0; 
        ee.net_pekel2_y(idx(ii),1:35) = 0; 
        continue, 
    end

    
    fi = fit(t(idxnan)'./365,dry_corr(ii,idxnan)',fitType);
    ee.net_pekel2(idx(ii)) = fi.p1;
    unc = confint(fi);
    ee.net_pekel2_unc(idx(ii),:) = unc(:,1);

    ee.net_pekel2_y(idx(ii),:) = dry_corr(ii,:);
    %accumarray(ee.net_pekel2_t'-1984,dry_corr(ii,:)',[],@nanmean);
    
end

%add to .mat file
save([f 'land_area_change' filesep 'GlobalDeltaData_AreaChange.mat'],'-struct','ee');