%DeltaEarthEngineShapefile
load([dropbox filesep 'WorldDeltas' filesep 'scripts' filesep 'GlobalDeltaData.mat']);
res=1;
remfun = @(lon) (rem(res*360-1+lon,res*360)+1);

%create buffer per delta
DeltaRadius = max(2,sqrt(1.07.*Discharge_prist.^0.7.*QRiver_prist.^0.45/pi))./111.*(cos(deg2rad(MouthLat)));
DeltaBuffer = 1000+(DeltaRadius*10000);

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
rate = mapshape(delta_shoreline_lon,delta_shoreline_lat,'BasinID',(BasinID),'Buffer',(DeltaBuffer));

shapewrite(rate,'GlobalDeltaShorelineData.shp')

%% put results back into mat file
load([dropbox filesep 'WorldDeltas' filesep 'GlobalDeltaData.mat']);

%file exported from earth engine
fileID = fopen('D:\Dropbox\WorldDeltas\EarthEngine\GlobalDeltaChange.csv','r');
data = textscan(fileID, '%q%f%f%f%f%f%f%f%f%q%[^\n\r]', 'Delimiter', ',', 'HeaderLines' ,1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
data = cell2mat(data(2:9)); %basin ID, deltaArea, change, deposition, erosion
data_aqua = data(:,[4 6 8])/28; %aquamonitor change deposition erosion per year
data_pekel = data(:,[3 5 7])/31; data_pekel(:,3) = data_pekel(:,3).*-1;

[~,idx] = ismember(BasinID,data(:,1));
idx(idx==0) = 2;

%put in structure
ee.net_pekel = data_pekel(idx,1);
ee.dep_pekel = data_pekel(idx,2);
ee.ero_pekel = data_pekel(idx,3);

ee.net_aqua = data_aqua(idx,1);
ee.dep_aqua = data_aqua(idx,2);
ee.ero_aqua = data_aqua(idx,3);

%add to .mat file
save('GlobalDeltaData.mat','ee','-append');
%% Delta Earth Engine shape file from NOAA maps
res=1;
remfun = @(lon) (rem(res*360-1+lon,res*360)+1);

%load shoreline
shore = shaperead('D:\GlobalDatasets\WorldCoastline\GSHHS_shp\f\GSHHS_f_L1.shp','Selector',{@(v1) (v1>100),'area'});
shore2 = struct2cell(shore);

%put lat and lon in imaginary numbers to make it a 1 line search
shore_cellidx = cellfun(@length,shore2(3,:));
shore_lon = remfun(cell2mat(shore2(3,:)));
shore_lat = cell2mat(shore2(4,:));
shore_both_all = shore_lon + 1i*shore_lat;
[shore_both,~,shoreline_idx] = unique(round(shore_both_all));
delta_both = (remfun(MouthLon)) + 1i*(MouthLat);
plot(shore_both(1:100),'o'), hold on, plot(delta_both(1),'or','MarkerFaceColor','r','MarkerEdgeColor','none')

%find closest shoreline for all 
idx = zeros(size(delta_both));
min_dist = zeros(size(delta_both));
for ii=1:length(delta_both),
   [min_dist(ii),idx(ii)] = min(abs(shore_both-delta_both(ii)));
end

delta_shoreline = cell(size(MouthLon));
delta_shoreline_lat = cell(size(MouthLon));
delta_shoreline_lon = cell(size(MouthLon));
for ii=1:length(delta_both),
    x = find(shoreline_idx==idx(ii));
    [~,idxx] = min(abs(shore_both_all(x)-delta_both(ii)));
    
    %which continents in this?
    co = find(x(idxx)<cumsum(shore_cellidx),1,'first');
    
    delta_shoreline{ii} = shore_both_all(x(idxx));
    %<DeltaRadius(ii));
    delta_shoreline_lat{ii} = shore_lat(delta_shoreline{ii});
    delta_shoreline_lon{ii} = shore_lon(delta_shoreline{ii});
end


%plot(shore_both_all(delta_shoreline),'o'), hold on, plot(delta_both,'or','MarkerFaceColor','r','MarkerEdgeColor','none')
%make shape file
rate = mapshape(delta_shoreline_lon,delta_shoreline_lat,'DeltaArea',DeltaArea,'BasinID',BasinID);

mapshow(rate)

shapewrite(rate,'GlobalDeltaShorelineData.shp')


%% get monthly data
out = load([dropbox filesep 'WorldDeltas' filesep 'scripts' filesep 'GlobalDeltaData.mat']);

%file exported from earth engine
filename = 'D:\Dropbox\WorldDeltas\EarthEngine\GlobalDeltaChangeMonth5.csv';
delimiter = ',';
formatSpec = '%q%q%q%q%q%q%[^\n\r]';
fileID = fopen(filename,'r');
T = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
%mississippiID = 4267691
BasinID = double(T{2}(2:end));
dry = arrayfun(@str2num,T{4}(2:end),'UniformOutput',0);
wet = arrayfun(@str2num,T{5}(2:end),'UniformOutput',0);


%re-order table
BasinID2 = int64(out.BasinID*10+out.Continent);
[~,idx] = (ismember(int64(BasinID2),BasinID));
dry = vertcat(dry{idx});
wet = vertcat(wet{idx});
BasinID = BasinID(idx);


%how to do this??
%we have km2 water, km2 land, km2 nothing
la_max = max(dry+wet,[],2);
fr_dry = dry./(dry+wet);
null = max(0,la_max-dry-wet);
dry_corr = (dry+fr_dry.*null);
dry_corr = dry_corr-nanmean(dry_corr,2);

%do yearly thing also
t = datenum(1985,3:420,1);
tt=(0:length(t)-1)./12;

save GlobalDeltaData_monthly BasinID BasinID2 dry_corr t tt 

%% 
load GlobalDeltaData_monthly

[yy,mm,~] = (datevec(t));
rate_y = zeros(size(dry_corr,1),35);
rate = zeros(size(dry_corr,1),1);
rate_unc = zeros(size(dry_corr,1),2);

fitType = fittype('poly1');

for ii=1:size(dry_corr,1),
    if mod(ii,100)==1, ii, end
    idx = ~isnan(dry_corr(ii,:));
    %[rate(ii,[1 2]),S] = polyfit(tt(~idx),dry_corr(ii,~idx),1);
    %[~,rate_unc(ii)] = polyval(rate(ii,:),0,S);
    if sum(idx)<2, rate(ii)=0; rate_unc(ii,:)=0; rate_y(ii,1:35) = 0; continue, end
    f = fit(tt(idx)',dry_corr(ii,idx)',fitType);
    rate(ii) = f.p1;
    unc = confint(f);
    rate_unc(ii,:) = unc(:,1);

    rate_y(ii,1:35) = accumarray(yy'-1984,dry_corr(ii,:)',[],@nanmean);
    
end
save GlobalDeltaData_monthly -append rate rate_unc rate_y

%% analysis
load GlobalDeltaData_monthly
out = load([dropbox filesep 'WorldDeltas' filesep 'scripts' filesep 'GlobalDeltaData.mat']);

%scatter(out.MouthLon,out.MouthLat,max(1,10*abs(p(:,1))),sign(p(:,1)))

%get change rate for couple of well-known deltas
[~,idx] = ismember(int64(out.delta_name_id*10+out.delta_name_continent),BasinID2);
[~,idx] = (ismember(int64(BasinID2(idx)),BasinID));
table(out.delta_name,idx,int64(BasinID(idx)),rate(idx),rate_unc(idx,1),rate_unc(idx,2))

out.Region(out.Region>20)=11;
table(out.Region_str',accumarray(out.Region,rate(:,1),[],@sum),accumarray(out.Region,rate(:,1),[],@(x) (mean(abs(x)))))

table(accumarray(out.Continent,rate(:,1),[],@(x) (mean(abs(x)))))

%compare against aquamonitor and other pekel change metric
scatter(out.ee.net_aqua,rate)

sqrt(mean((rate-out.ee.net_aqua).^2))

T = table(BasinID2,rate,rate_unc,rate_y);
writetable(T,'GlobalDeltaData_monthly.csv')





