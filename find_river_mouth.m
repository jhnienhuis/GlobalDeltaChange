function find_river_mouth
%get started with obtaining all global deltas

%%new river mouth algorithm, all cells
continents = {'na','af','ca','sa','eu','as','au'};

res=240;
remfun = @(lon) (rem(res*360-1+lon,res*360)+1);

%load SRTM data to check if basin is within 40 m.
load('D:\GlobalDatasets\GlobalDEM\CoastalZone40m.mat');

%load coastal topology, durr et al
durr = shaperead('D:\GlobalDatasets\DIVA\typology_coastline.shp','UseGeoCoords',true);
durr_lon = remfun(res.*[durr(:).Lon]);
durr_lat = res.*(90+[durr(:).Lat]);
durr_size = arrayfun(@(x) length(x.Lon),durr);
durr_type = cell2mat(arrayfun(@(x,y) (x.*ones(1,y)),[durr(:).FIN_TYP],durr_size','UniformOutput',false));
durr_nan = isnan(durr_lon); durr_lon(durr_nan) = []; durr_lat(durr_nan) = []; durr_type(durr_nan) = [];
durr_both = single(durr_lon + 1i*durr_lat);

%load DIVA dataset
diva = shaperead('D:\GlobalDatasets\DIVA\cls_p18_2.shp','UseGeoCoords',true);
diva_lon = remfun(res.*[diva(:).Lon]);
diva_lat = res.*(90+[diva(:).Lat]);
diva_size = arrayfun(@(x) length(x.Lon),diva);
diva_type = cell2mat(arrayfun(@(x,y) (x.*ones(1,y)),[diva(:).CPC],diva_size','UniformOutput',false));
diva_cont = cell2mat(arrayfun(@(x,y) (x.*ones(1,y)),[diva(:).GVARID],diva_size','UniformOutput',false));
diva_nan = isnan(diva_lon); diva_lon(diva_nan) = []; diva_lat(diva_nan) = []; diva_type(diva_nan) = []; diva_cont(diva_nan) = [];
diva_both = single(diva_lon + 1i*diva_lat);
clear durr diva

%now do HydroSheds retrieve all deltas
for jj=1:length(continents)
%load accumulated drainage area (# cells)
d = ['D:\GlobalDatasets\HydroSheds\' continents{jj} '_acc_15s_bil\' continents{jj} '_acc_15s'];
fileID = fopen([d '.hdr'],'r');
hdr = textscan(fileID, '%s%s%[^\n\r]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true); fclose(fileID);
hdr = str2double(hdr{2});
ulX = round(remfun(res*hdr(11)));
ulY = round(res*(90+hdr(12)));
a = multibandread([d '.bil'],[hdr(3) hdr(4) hdr(5)],'int32',0,'bil','ieee-le');
%latitude and longitude vectors of grid
lat_vec = (ulY/res-90):(-1/res):((ulY-hdr(3))/res-90);
lon_vec = (ulX/res):(1/res):((ulX+hdr(4))/res);

%get points in HydroSheds gridded drainage area dataset and convert to actual km2!
areapercell = 6371.^2.*2*pi/360/res*(sin(deg2rad(lat_vec(1:end-1)))-sin(deg2rad(lat_vec(2:end))))';
disp('Loading drainage accumulation...')
a = int32(bsxfun(@times,a,areapercell));

%remove river deltas with drainage area under 50 km2
a(a<50) = 0;

x = size(a,1); y = size(a,2);

%find all basins that empty in ocean (bil==0, necessary to get specific
%basin lat/lon of river mouth and to eliminate deep endorheic lakes
d = ['D:\GlobalDatasets\HydroSheds\' continents{jj} '_dir_15s_bil\' continents{jj} '_dir_15s'];
fileID = fopen([d '.hdr'],'r');
hdr = textscan(fileID, '%s%s%[^\n\r]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true); fclose(fileID);
hdr = str2double(hdr{2});
ulX = res*rem(hdr(11)+360,360);
ulY = res*(90+hdr(12));
rm = int32(multibandread([d '.bil'],[hdr(3) hdr(4) hdr(5)],'int8',0,'bil','ieee-le')==0);

%find all river mouths
rm =  rm .* a;
disp('Calculating regional maxima...')
rm = rm .* int32(imregionalmax(rm)); %ocean emptying points that are regional maxima in terms of accumulated drainage area

[rmlat,rmlon] = find(rm);

%convert to actual lat/lon
MouthLat_c = lat_vec(rmlat);
MouthLon_c = rem(lon_vec(rmlon),360);

%turn into linear vector
mouth_both = int32((res*MouthLon_c) + 1i*res*(90+MouthLat_c));

%get basin area for all river mouths
BasinArea_c = a(sub2ind(size(a),rmlat,rmlon))'; %km2

disp(['River Mouths found: ' num2str(length(mouth_both))])
clearvars a rm

%now do the matching with the basin outline vector dataset
disp('Loading Basins...')
basins = shaperead(['D:\GlobalDatasets\HydroSheds\' continents{jj} '_bas_15s_beta\' continents{jj} '_bas_15s_beta.shp'],'Selector',{@(v1) (v1>40),'AREA_SQKM'});

BAS_X = floor(remfun(res.*[basins(:).X]));
BAS_Y = floor(res.*(90+[basins(:).Y]));
BAS_Both = single(BAS_X + 1i*BAS_Y);
BAS_ID_V = zeros(size(BAS_Y),'int32');
BAS_AREA_V = zeros(size(BAS_Y),'int32');
BAS_AREA = [basins(:).AREA_SQKM];
BAS_ID = [basins(:).BASIN_ID];

d=0;
for ii=1:length(basins),
    sz = length(basins(ii).X);
    BAS_ID_V(d+(1:sz)) = basins(ii).BASIN_ID;
    BAS_AREA_V(d+(1:sz)) = basins(ii).AREA_SQKM;
    d=d+sz;
end

disp('Matching basins...')
%first do match, then coastal check!
BasinID_c = zeros(size(MouthLat_c));
BasinAreaBAS_c = zeros(size(MouthLat_c));

[a,idx] = ismember(mouth_both,int32(BAS_Both));
BasinAreaBAS_c(a) = BAS_AREA_V(idx(idx>0));

BasinID_c(a) = BAS_ID_V(idx(idx>0));

BasinID_c(abs(log10(single(BasinAreaBAS_c))-log10(single(BasinArea_c)))>0.3) = 0;
BAS_Both(ismember(BAS_ID_V,BasinID_c)) = nan;
%do combination of minimum distance and equal drainage area

for ii=1:length(BasinID_c)
    if BasinID_c(ii)>0,continue, end
    
    AreaDiff = abs(log10(single(BAS_AREA_V))-log10(single(BasinArea_c(ii))))<0.4;

    [a,idx] = min(abs(single(mouth_both(ii))-BAS_Both(AreaDiff)));
    idx = find(cumsum(AreaDiff)==idx);
    if isempty(a) || a>10, BasinAreaBAS_c(ii) = nan; continue, end 
    BasinID_c(ii) = BAS_ID_V(idx);
    BasinAreaBAS_c(ii) = BAS_AREA_V(idx);
    BAS_Both(BAS_ID_V==BAS_ID_V(idx)) = nan;
end

%accumulate non coastal zone cells for all basins coordinates. 
remnan = ~isnan(BAS_X);
czbasins = accumarray(single(BAS_ID_V(remnan))',cz(sub2ind(size(cz),BAS_Y(remnan),BAS_X(remnan))));

%find if any basins have no non-coastal cells
idx = (~ismember(BasinID_c,find(czbasins==0)) | BasinArea_c > 1000) & BasinID_c>0;

disp(['Lowland drainage eliminated: ' num2str(sum(~idx))])

%get rid of fjords, only for small drainage areas
durr_idx = false(size(BasinID_c));
for ii=1:length(BasinID_c)
    %find 5 coastal types closest to river mouth, only discard if all are
    %fjord (4) and is relatively close (10 degrees) to durr dataset
    durr_dis = abs(single(mouth_both(ii))-durr_both);
    [durr_dis,durr_sort] = sort(durr_dis);
    durr_idx(ii) = ~(all(durr_type(durr_sort(1:5))==4) && durr_dis(1)<2400);
    
end
durr_idx = durr_idx | BasinArea_c > 500;
disp(['Coastal Type Eliminated: ' num2str(sum(~durr_idx))])

idx = idx & durr_idx;
%{
basins2 = shaperead(['D:\GlobalDatasets\HydroSheds\' continents{jj} '_bas_15s_beta\' continents{jj} '_bas_15s_beta.shp'],'Selector',{@(v1) (v1<40),'AREA_SQKM'});
BAS_X2 = floor(remfun(res.*[basins2(:).X]));
BAS_Y2 = floor(res.*(90+[basins2(:).Y]));

BASBOTH = setdiff(BAS_X2+i*BAS_Y2,BAS_X+i*BAS_Y,'stable');

data = load('D:\GlobalDatasets\GlobalDEM\CoastalZone_int8.mat');
xl = [262 275]*240;
yl = [107 112]*240;


imagesc(xl(1):xl(2),yl(1):yl(2),data.cz(yl(1):yl(2),xl(1):xl(2))-1),hold on
contour(xl(1):xl(2),yl(1):yl(2),cz(yl(1):yl(2),xl(1):xl(2)),1,'k')
plot(real(BASBOTH),imag(BASBOTH),'r')
plot(BAS_X,BAS_Y,'w')
scatter(remfun(res*MouthLon(~idx)),res*(90+MouthLat(~idx)),30,'w','filled','MarkerEdgeColor','k'), axis xy
scatter(remfun(res*MouthLon(idx)),res*(90+MouthLat(idx)),30,'g','filled','MarkerEdgeColor','k'), axis xy
set(gca,'CLim',[0 127],'DataAspectRatio',[1 1 1],'XLim',xl,'YLim',yl,...
'XTick',xl(1):240:xl(2),'YTick',yl(1):240:yl(2),'XTickLabel',(262:275),'YTickLabel',(17:22))
colormap(jet)
colorbar('east','Color','w')

xl = [264 267]*240;
yl = [108 109]*240;
set(gca,'CLim',[0 127],'DataAspectRatio',[1 1 1],'XLim',xl,'YLim',[107.5*240 109*240],...
'XTick',xl(1):240:xl(2),'YTick',yl(1):240:yl(2),'XTickLabel',(264:267),'YTickLabel',(18:19))
%}

BasinID_a{jj} = BasinID_c(idx);
BasinArea_a{jj} = BasinArea_c(idx);
BasinAreaBAS_a{jj} = BasinAreaBAS_c(idx);
MouthLat_a{jj} = MouthLat_c(idx);
MouthLon_a{jj} = MouthLon_c(idx);

%check result
%figure, scatter(log10(single(BasinArea)),log10(single(BasinAreaBAS)))


end



%% find river mouths >59 lat based on WBMSed

%do this with the WBMSed dataset, based on the 6min data that is available
%north of the wall
direc = 'D:\GlobalDatasets\WBMSed\Prist\';
res=240;
remfun = @(lon) (rem(res*360-1+lon,res*360)+1);

%load discharge (seems better at punctuating at the coast
fname = 'Discharge\Global_Discharge_FixedB+Prist_06min_aTS2010.nc';
Discharge = ncread([direc fname],'discharge');
Discharge(isnan(Discharge)) = -999;

lat = ncread([direc fname],'latitude');
lon = ncread([direc fname],'longitude');

bqart_a = fliplr(imread('D:\GlobalDatasets\WBMSed\bqart_sed-flow_dir\bqart_a.tif')');
lin_idx = find(imregionalmax(bqart_a));
[lon_idx,lat_idx] = ind2sub(size(bqart_a),lin_idx);
rm_lat = lat(lat_idx);
rm_lon = lon(lon_idx);
dr_a = bqart_a(lin_idx);
rm_dis = Discharge(lin_idx);
rm_both = round(5*rm_lon)+1i*round(5*rm_lat);

Dist = load(['D:\GlobalDatasets\WBMSed' filesep 'Dist_Total.mat']);
Prist = load(['D:\GlobalDatasets\WBMSed' filesep 'Prist_Total.mat']);

Dist.dis_total = Dist.dis_total(lin_idx);
Dist.sed_total = Dist.sed_total(lin_idx);

Prist.dis_total = Prist.dis_total(lin_idx);
Prist.sed_total = Prist.sed_total(lin_idx);

%only select above 59 deg north and bigger than 1000 km2
idx_arctic = rm_lat>59 & dr_a>1000;

load('D:\GlobalDatasets\WorldCoastline\shore_vec.mat');
shore_lat = cell2mat(shore_vec(4,:));
shore_lon = cell2mat(shore_vec(3,:));
shore_both_12m = round(5*shore_lon) + 1i*round(5*shore_lat);

river_mouth_idx = zeros(size(rm_both));
river_mouth_lat = zeros(size(rm_both));
river_mouth_lon = zeros(size(rm_both));


for jj=1:length(rm_lon)
    
    if mod(jj,1000)==1, jj, end
    
    if idx_arctic(jj)
    
    [shore_pixels] = ismember(int16(shore_both_12m),int16(rm_both(jj)));
        
    if any(shore_pixels)
        
        %[basin_pixels] = ismember(int16(basin_both),int16(shore_both_12m(shore_pixels)));
        lis = find(shore_pixels);
        [~,idx] = min((rm_lon(jj)-shore_lon(shore_pixels)).^2+...
            (rm_lat(jj)-shore_lat(shore_pixels)).^2);
        river_mouth_idx(jj) = lis(idx);
        
        river_mouth_lat(jj) = shore_lat(river_mouth_idx(jj));
        river_mouth_lon(jj) = shore_lon(river_mouth_idx(jj));
    elseif dr_a(jj) > 100000,
        river_mouth_idx(jj) = 1;
        river_mouth_lat(jj) = rm_lat(jj);
        river_mouth_lon(jj) = rm_lon(jj);
    end
    end
end

wbm_idx_arctic = lin_idx(river_mouth_idx~=0);
MouthLat_c = river_mouth_lat(river_mouth_idx~=0);
MouthLon_c = river_mouth_lon(river_mouth_idx~=0);

QRiver_prist_arc = Prist.sed_total(river_mouth_idx~=0);
Discharge_prist_arc = Prist.dis_total(river_mouth_idx~=0);
QRiver_dist_arc = Dist.sed_total(river_mouth_idx~=0);
Discharge_dist_arc = Dist.dis_total(river_mouth_idx~=0);
BasinArea_c = dr_a(river_mouth_idx~=0);

%get rid of fjords, only for small drainage areas
%load coastal topology, durr et al
durr = shaperead('D:\GlobalDatasets\DIVA\typology_coastline.shp','UseGeoCoords',true);
durr_lon = remfun(res.*[durr(:).Lon]);
durr_lat = res.*(90+[durr(:).Lat]);
durr_size = arrayfun(@(x) length(x.Lon),durr);
durr_type = cell2mat(arrayfun(@(x,y) (x.*ones(1,y)),[durr(:).FIN_TYP],durr_size','UniformOutput',false));
durr_nan = isnan(durr_lon); durr_lon(durr_nan) = []; durr_lat(durr_nan) = []; durr_type(durr_nan) = [];
durr_both = single(durr_lon + 1i*durr_lat);
durr_idx = false(size(BasinArea_c));

mouth_both = remfun(res*MouthLon_c) + 1i*res*(90+MouthLat_c);
for ii=1:length(BasinArea_c)
    %find 5 coastal types closest to river mouth, only discard if all are
    %fjord (4) and is relatively close (10 degrees) to durr dataset
    durr_dis = abs(single(mouth_both(ii))-durr_both);
    [durr_dis,durr_sort] = sort(durr_dis);
    durr_idx(ii) = ~(all(durr_type(durr_sort(1:5))==4) && durr_dis(1)<2400);
    
end
idx = durr_idx | BasinArea_c > 10000;
disp(['Coastal Type Eliminated: ' num2str(sum(~durr_idx))])

BasinArea_a{8} = BasinArea_c(idx);
MouthLat_a{8} = MouthLat_c(idx);
MouthLon_a{8} = MouthLon_c(idx);
QRiver_prist_arc = QRiver_prist_arc(idx);
Discharge_prist_arc = Discharge_prist_arc(idx);
QRiver_dist_arc = QRiver_dist_arc(idx);
Discharge_dist_arc = Discharge_dist_arc(idx);

%% combine continent files

%merge all files into one mat file
continents = {'na','af','ca','sa','eu','as','au','ao'};
remfun = @(lon) (rem(360-1+lon,360)+1);

k=0;
for jj = 1:8,
    %load([dropbox filesep 'WorldDeltas\scripts\Rivermouth' continents{jj} '.mat'],'BasinArea');
    sz = length(BasinArea{jj});
    Continent(k+(1:sz)) = jj;
    k=k+sz;
end
Continent = Continent';
BasinArea = zeros(k,1);
MouthLat = zeros(k,1);
MouthLon = zeros(k,1);
QRiver_dist = zeros(k,1);
QRiver_prist = zeros(k,1);
Discharge_dist = zeros(k,1);
Discharge_prist = zeros(k,1);
BasinID = zeros(k,1);

for jj = 1:8,
    %out = load([dropbox filesep 'WorldDeltas\scripts\Rivermouth' continents{jj} '.mat']);
    
    MouthLon(Continent==jj) = MouthLon_a{jj};
    MouthLat(Continent==jj) = MouthLat_a{jj};
        
    if jj==8,
        BasinID(Continent==jj) = 1:length(MouthLat_a{jj});
        
        QRiver_dist(Continent==jj) = QRiver_dist_arc;
        QRiver_prist(Continent==jj) = QRiver_prist_arc;
        Discharge_dist(Continent==jj) = Discharge_dist_arc;
        Discharge_prist(Continent==jj) = Discharge_prist_arc;
        
    else,
        BasinID(Continent==jj) = BasinID_a{jj};    
    end
        
    BasinArea(Continent==jj) = BasinArea_a{jj};
    
end



% add continent flag
diva = shaperead('D:\OneDrive - Universiteit Utrecht\DIVA\cls_p18_2.shp','UseGeoCoords',true);
remfun = @(lon) (rem(360-1+lon,360)+1);
diva_lon = remfun([diva(:).Lon]);
diva_lat = ([diva(:).Lat]);
diva_size = arrayfun(@(x) length(x.Lon),diva);
diva_type = cell2mat(arrayfun(@(x,y) (x.*ones(1,y)),[diva(:).CPC],diva_size','UniformOutput',false));
diva_cont = cell2mat(arrayfun(@(x,y) (x.*ones(1,y)),[diva(:).GVARID],diva_size','UniformOutput',false));
[~,ia] = unique(int32(diva_lat) + 1000*int32(diva_lon));
diva_lon = diva_lon(ia); diva_lat = diva_lat(ia); diva_type = diva_type(ia); diva_cont = diva_cont(ia);
%scatter(diva_lon,diva_lat,40,diva_cont)

diva_both = diva_lon + 1i*diva_lat;
delta_both = MouthLon + 1i*MouthLat;
[~,x] = min(abs(repmat(delta_both,1,length(diva_both))-diva_both),[],2);

%{
Region = int32(diva_cont(x))';
%scatter(MouthLon,MouthLat,30,cont)
hold on,
l = unique(diva_cont);
for ii=1:length(l),
    
    scatter(MouthLon(Region==l(ii)),MouthLat(Region==l(ii)),30,Region(Region==l(ii)),'filled')
end

Region(Region==0) = 3;
Region(Region==23) = 18;
Region(Region>14) = Region(Region>14)-1;
Region_str = {'East Africa', 'South Asia','West Africa','Baltic','Eastern Central America','Western Central America','Carribean','Russia','East Asia','Middle-East','Northern Africa','Eastern North America','Western North America','Northern Mediterranean','Western Europe','Oceania','Pacific Islands','Eastern South America','Western South America','Southeast Asia'};
%}

% save files
BasinID2 = int64(BasinID.*10+Continent);

save([dropbox filesep 'github\GlobalDeltaChange\GlobalDeltaData.mat'],'BasinArea','MouthLat','MouthLon','BasinID','BasinID2','Continent','QRiver_prist','QRiver_dist','Discharge_prist','Discharge_dist');