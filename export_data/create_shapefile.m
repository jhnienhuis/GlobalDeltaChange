function create_shapefile

% load all deltas
load([dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep 'GlobalDeltaData.mat'])
ee = load([dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep 'land_area_change' filesep 'GlobalDeltaData_AreaChange.mat']);

MouthLon = rem(MouthLon+360,360);
mouthboth = MouthLon+1i*MouthLat;

%get list of basin id
opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["HYBAS_ID", "MAIN_BAS"];
opts.VariableTypes = ["int64", "int64"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
mainbasid = table2array(readtable("D:\OneDrive - Universiteit Utrecht\HydroSheds\BasinATLAS_Data_v10_shp\BasinATLAS_v10_shp\main_bas_id.csv", opts));
clear opts

% match rivers to deltas (again...)
con = {'RiverATLAS_v10_eu','RiverATLAS_v10_na','RiverATLAS_v10_ar','RiverATLAS_v10_si','RiverATLAS_v10_au','RiverATLAS_v10_as','RiverATLAS_v10_af','RiverATLAS_v10_sa_north','RiverATLAS_v10_sa_south'};

for ii=1:length(con)
    ii
%rivers = shaperead(['D:\OneDrive - Universiteit Utrecht\HydroSheds\RiverATLAS_Data_v10_shp\RiverATLAS_v10_shp\' con{ii}],'Attributes',{'HYBAS_L12','UPLAND_SKM'},'Selector',{@(v1,v2,v3) (v1==0) && (v2==0) && (v3>50),'NEXT_DOWN','ENDORHEIC','UPLAND_SKM'},'UseGeoCoords',true);
rivers = shaperead(['D:\OneDrive - Universiteit Utrecht\HydroSheds\RiverATLAS_Data_v10_shp\RiverATLAS_v10_shp\' con{ii} '_filt'],'UseGeoCoords',true);
%last one is on coastline
riverLat = arrayfun(@(x) (x.Lat(end-1)),rivers);
riverLon = arrayfun(@(x) (x.Lon(end-1)),rivers);
riverLon = rem(riverLon+360,360);

riverID{ii} = [rivers.HYRIV_ID];
riverBasID{ii} = [rivers.HYBAS_L12];
riverArea{ii} = [rivers.UPLAND_SKM];
riverboth{ii} = rot90(riverLon+1i*riverLat);
end
riverID = int64([riverID{:}]);
riverBasID = int64([riverBasID{:}]);
riverArea = [riverArea{:}];
riverboth = [riverboth{:}];


% do fancy minimum
blub = (abs(riverboth-mouthboth)+...
    (abs((riverArea-BasinArea)./BasinArea.*log10(BasinArea)))+...
    (abs(riverboth-mouthboth)>8).*10);
idx = zeros(size(blub,1),1);
for ii=1:size(blub,1)
    [~,idx(ii)] = min(blub(ii,:));
    blub(:,idx(ii)) = inf;
end

%find 310319
%{
plot(riverBoth,'*r'), hold on
plot(mouthboth(Continent==8),'o'), hold on
plot(riverBoth(idx),'d','MarkerSize',15);
%}
% match
Basins = riverBasID(idx);
RiverBasinArea = riverArea(idx);
riverboth = riverboth(idx);
%link river mouth to basin, seems correct?
%plot([real(riverboth)' real(mouthboth)]',[imag(riverboth)' imag(mouthboth)]','-o')
[~,ida] = ismember(Basins,mainbasid(:,1));

BasinID_ATLAS = mainbasid(ida,2);

basins = shaperead('D:\OneDrive - Universiteit Utrecht\HydroSheds\BasinATLAS_Data_v10_shp\BasinATLAS_v10_shp\BasinATLAS_v10_lev12_mainbas4','UseGeoCoords',true);
basins_add = shaperead('D:\OneDrive - Universiteit Utrecht\HydroSheds\BasinATLAS_Data_v10_shp\BasinATLAS_v10_shp\BasinATLAS_v10_lev12_mainbas4add','UseGeoCoords',true);

[~,idx] = setxor(int64([basins_add.MAIN_BAS]),int64([basins.MAIN_BAS]));
basins = [basins; basins_add(idx)];
main_bas = int64([basins.MAIN_BAS]);

[~,ida] = ismember(BasinID_ATLAS,main_bas);


%sum(ismember(BasinID_ATLAS2,main_bas))
%plot(mouthboth(Continent==8),'or'), hold on
%plot([basins(ismember(main_bas,BasinID_ATLAS2)).Lon],[basins(ismember(main_bas,BasinID_ATLAS2)).Lat])
%hold on, plot(mouthboth(~ismember(BasinID_ATLAS2,main_bas)),'or')
basins_deltas = basins(ida);

basins_lat = {basins_deltas.Lat};
basins_lon = {basins_deltas.Lon};

for ii=1:length(basins_lat),
    if rem(ii,100)==1, ii, end
    idx = [1 find(isnan(basins_lat{ii}))];
    [~,idx2] = max(diff(idx));
    
    %if length(idx)>1, [~,idx2] = max(diff(idx)); idx3 = [idx(idx2)+1 idx(idx2+1)]; else, idx3 = [1 idx]; end
    %idx2 = 1:5:length(basins_arctic_lat{ii});
    %idx = union(idx2,idx);
    tol = ceil(log10(max(idx)))*0.002;
    [basins_lat{ii}, basins_lon{ii}] = reducem(basins_lat{ii}(idx(idx2):([-1 0]+idx(idx2+1)))', basins_lon{ii}(idx(idx2):([-1 0]+idx(idx2+1)))',tol);
    
    
    
    %basins_lat{ii} = basins_lat{ii}(idx(idx2):10:([-1 0]+idx(idx2+1)));
    %basins_lon{ii} = basins_lon{ii}(idx(idx2):10:([-1 0]+idx(idx2+1)));
end
%save GlobalDeltaData -append BasinID_ATLAS 
%}
% write to shapefiles
MouthLon(MouthLon>180) = MouthLon(MouthLon>180)-360;
channel_len_lon(channel_len_lon>180) = channel_len_lon(channel_len_lon>180)-360;
channel_len_lat(channel_len_lat==0) = NaN; channel_len_lon(channel_len_lon==0) = NaN;
idx = isnan(channel_len_lat(:,1));
channel_len_lat(idx,1) = MouthLat(idx)+0.001; channel_len_lon(idx,1) = MouthLon(idx)+0.001;

%shapefile limited to 10 character attribute names..
p = geoshape(basins_lat,basins_lon,'BasinID',BasinID,'BasinID2',double(BasinID2),'BasinID_ATL',double(BasinID_ATLAS),'BasinArea',BasinArea,'Dis_pris',Discharge_prist,'Dis_dist',Discharge_dist,'delta_name',delta_name,'Dis_tide',Discharge_tide,...
    'MouthLat',MouthLat,'MouthLon',MouthLon,'QRiver_dis',QRiver_dist,'QRiver_pri',QRiver_prist,'QTide',QTide,'QWave',QWave,'Hs',Hs,'TidalAmp',TidalAmp,...
    'DeltaAreaGa',ee.dep_aqua,'DeltaAreaLo',ee.ero_aqua);
p.Geometry = 'polygon';
fname = 'GlobalDeltaBasins';
shapewrite(p,fname)
zip([fname '_shp'],{[fname '.dbf'],[fname '.shx'],[fname '.shp']})
delete([fname '.dbf'],[fname '.shx'],[fname '.shp'])

b = geopoint(MouthLat,MouthLon,'BasinID',double(BasinID),'BasinID2',double(BasinID2),'BasinID_ATL',double(BasinID_ATLAS),'BasinArea',BasinArea,'Dis_pris',Discharge_prist,'Dis_dist',Discharge_dist,'delta_name',delta_name,'Dis_tide',Discharge_tide,...
    'MouthLat',MouthLat,'MouthLon',MouthLon,'QRiver_dis',QRiver_dist,'QRiver_pri',QRiver_prist,'QTide',QTide,'QWave',QWave,'Hs',Hs,'TidalAmp',TidalAmp,...
    'DeltaAreaGa',ee.dep_aqua,'DeltaAreaLo',ee.ero_aqua);
b.Geometry = 'point';
fname = 'GlobalDeltaMouth';
shapewrite(b,fname)
zip([fname '_shp'],{[fname '.dbf'],[fname '.shx'],[fname '.shp']})
delete([fname '.dbf'],[fname '.shx'],[fname '.shp'])

c = geoshape(num2cell([MouthLat channel_len_lat],2),num2cell([MouthLon channel_len_lon],2),'BasinID',BasinID,'BasinID2',double(BasinID2),'BasinID_ATL',double(BasinID_ATLAS),'BasinArea',BasinArea,'Dis_pris',Discharge_prist,'Dis_dist',Discharge_dist,'delta_name',delta_name,'Dis_tide',Discharge_tide,...
    'MouthLat',MouthLat,'MouthLon',MouthLon,'QRiver_dis',QRiver_dist,'QRiver_pri',QRiver_prist,'QTide',QTide,'QWave',QWave,'Hs',Hs,'TidalAmp',TidalAmp,...
    'DeltaAreaGa',ee.dep_aqua,'DeltaAreaLo',ee.ero_aqua);
c.Geometry = 'line';
fname = 'GlobalDeltaRivers';
shapewrite(c,fname)
zip([fname '_shp'],{[fname '.dbf'],[fname '.shx'],[fname '.shp']})
delete([fname '.dbf'],[fname '.shx'],[fname '.shp'])

%% create shapefile for earth engine app
p = geoshape(basins_lat,basins_lon,'BasinID2',double(BasinID2),'BasinArea',BasinArea,'Dis_pris',Discharge_prist,'Dis_dist',Discharge_dist,'delta_name',delta_name,'Dis_tide',Discharge_tide,...
    'MouthLat',MouthLat,'MouthLon',MouthLon,'QRiver_dis',QRiver_dist,'QRiver_pri',QRiver_prist,'QTide',QTide,'QWave',QWave);

p.rate = ee.net_pekel2;
for ii=1:36,
    p.(['rate_y_' num2str(ii)]) = ee.net_pekel2_y(:,ii);
end

p.Geometry = 'polygon';
fname = 'GlobalDeltaBasins_app';
shapewrite(p,fname)
zip([fname '_shp'],{[fname '.dbf'],[fname '.shx'],[fname '.shp']})
delete([fname '.dbf'],[fname '.shx'],[fname '.shp'])