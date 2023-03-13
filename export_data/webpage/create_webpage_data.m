function create_webpage_data
%
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData_Web','DeltaSL_SSP126','DeltaSL_SSP245','DeltaSL_SSP585','DeltaSL')
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaProfile','bed_h'); 
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaArea.mat','src','delta_area','apex_latlon','sho1_latlon','mouth_latlon','sho2_latlon');
load('D:\Drive\github\GlobalDeltaChange\GlobalDeltaData.mat','QRiver_dist','MouthLon','MouthLat','delta_name','Discharge_prist','Discharge_tide','Discharge_dist','QTide','QWave','QRiver_prist','BasinID2','BasinID_ATLAS');
load('D:\Drive\github\GlobalDeltaChange\land_area_change\GlobalDeltaData_AreaChange.mat','net_aqua');

%remove caspian sea and other non deltas and deltas > 10km2
idx = (~isnan(bed_h) & bed_h<0) & ~(MouthLon>46 & MouthLon<56 & MouthLat>34 & MouthLat<48) & delta_area>1e7 &src==1;

%addpath('D:\Drive\github\GlobalDeltaSeaLevel')

%% add drainage basins outlines
basins = shaperead('D:\OneDrive - Universiteit Utrecht\HydroSheds\BasinATLAS_Data_v10_shp\BasinATLAS_v10_shp\BasinATLAS_v10_lev12_mainbas4','UseGeoCoords',true);
basins_add = shaperead('D:\OneDrive - Universiteit Utrecht\HydroSheds\BasinATLAS_Data_v10_shp\BasinATLAS_v10_shp\BasinATLAS_v10_lev12_mainbas4add','UseGeoCoords',true);

[~,kk] = setxor(int64([basins_add.MAIN_BAS]),int64([basins.MAIN_BAS]));
basins = [basins; basins_add(kk)];
main_bas = int64([basins.MAIN_BAS]);

[~,ida] = ismember(BasinID_ATLAS(idx),main_bas);

basins_deltas = basins(ida);

basins_lat = {basins_deltas.Lat};
basins_lon = {basins_deltas.Lon};

for ii=1:length(basins_lat),
    if rem(ii,100)==1, ii, end
    kk = [1 find(isnan(basins_lat{ii}))];
    [~,idx2] = max(diff(kk));
    
    tol = ceil(log10(max(kk)))*0.002;
    [basins_lat{ii}, basins_lon{ii}] = reducem(basins_lat{ii}(kk(idx2):([-1 0]+kk(idx2+1)))', basins_lon{ii}(kk(idx2):([-1 0]+kk(idx2+1)))',tol);
    
end

p = geoshape(basins_lat,basins_lon,'BasinID2',double(BasinID2(idx)));
p.Geometry = 'polygon';

%p = p(res.idx);

fname = 'GlobalDeltaBasinsWeb';
shapewrite(p,fname)
zip([fname '_shp'],{[fname '.dbf'],[fname '.shx'],[fname '.shp']})
delete([fname '.dbf'],[fname '.shx'],[fname '.shp'])

%% create delta land shapefile
delta_area_lat = cell(size(BasinID2));
delta_area_lon = cell(size(BasinID2));
% write to shapefiles
for ii=1:length(BasinID2),
    delta_area_lat{ii} = [apex_latlon(ii,1) sho1_latlon(ii,1) mouth_latlon(ii,1) sho2_latlon(ii,1)];
    delta_area_lon{ii} = [apex_latlon(ii,2) sho1_latlon(ii,2) mouth_latlon(ii,2) sho2_latlon(ii,2)];
end

delta_name(strlength(delta_name)==0) = "N/A";
%shapefile limited to 10 character attribute names..
p = geoshape(delta_area_lat,delta_area_lon,'BasinID2',double(BasinID2));
p.Geometry = 'polygon';

p = p(idx);

fname = 'D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaAreaWeb';
shapewrite(p,fname)
zip([fname '_shp'],{[fname '.dbf'],[fname '.shx'],[fname '.shp']})
delete([fname '.dbf'],[fname '.shx'],[fname '.shp'])

%% delta slr table

idx = (~isnan(bed_h) & bed_h<0) & ~(MouthLon>46 & MouthLon<56 & MouthLat>34 & MouthLat<48) & delta_area>1e7 &src==1;

%src = [126+zeros(size(BasinID2));245+zeros(size(BasinID2));585+zeros(size(BasinID2))];
%slr = [DeltaSLR DeltaSLR_SSP126;DeltaSLR DeltaSLR_SSP245;DeltaSLR DeltaSLR_SSP585].*1e3; %convert to mm/yr

yrs = repmat((1900:50:2300)',size(BasinID2));
[SSP126, SSP245, SSP585] = deal([DeltaSL DeltaSL_SSP126]',[DeltaSL DeltaSL_SSP245]',[DeltaSL DeltaSL_SSP585]');
[SSP126 ,SSP245, SSP585] = deal(reshape(SSP126,[9*length(BasinID2) 1]),reshape(SSP245,[9*length(BasinID2) 1]),reshape(SSP585,[9*length(BasinID2) 1]));

id = repmat(BasinID2',9,1); id = id(:);
idx2 = repmat(idx',9,1); idx2 = idx2(:);

T = [table(id), table(yrs), table(SSP126,SSP245,SSP585)];
T.Properties.VariableNames = (["BasinID2";"yrs";"SSP126";"SSP245";"SSP585"]);
writetable(T(idx2,:),'GlobalDelta_DeltaSLRWeb.csv')



%% delta land area change table
ff = 365*24*3600/1600;
fr = 0.8;
idx = (~isnan(bed_h) & bed_h<0) & ~(MouthLon>46 & MouthLon<56 & MouthLat>34 & MouthLat<48) & delta_area>1e7 &src==1;

%validation
%land_pred = (ff.*QRiver_dist.*fr - (delta_area.*(DeltaSLR+DeltaSub)))./-bed_h;
%nansum(land_pred(idx))./1e6
%sum(net_aqua(idx))
%corrcoef(net_aqua(idx),land_pred(idx)./1e6).^2

%future


t = [1985 2015 2050 2100 2150 2200 2250 2300];
slr = diff([zeros(size(BasinID2)) DeltaSL_SSP126],1,2)./50; %slr in m/yr
dA = (ff.*QRiver_dist.*fr - (delta_area.*slr))./-bed_h;
Ayr_126 = delta_area+[zeros(size(MouthLat)), cumsum([net_aqua.*1e6 dA].*diff(t),2)];
%plot(t,(nansum(Ayr_585(idx,:),1))./1e6,'-o'), hold on

slr = diff([zeros(size(BasinID2)) DeltaSL_SSP245],1,2)./50; %slr in m/yr
dA = (ff.*QRiver_dist.*fr - (delta_area.*slr))./-bed_h;
Ayr_245 = delta_area+[zeros(size(MouthLat)), cumsum([net_aqua.*1e6 dA].*diff(t),2)];

slr = diff([zeros(size(BasinID2)) DeltaSL_SSP585],1,2)./50; %slr in m/yr
dA = (ff.*QRiver_dist.*fr - (delta_area.*slr))./-bed_h;
Ayr_585 = delta_area+[zeros(size(MouthLat)), cumsum([net_aqua.*1e6 dA].*diff(t),2)];

%src = [126+zeros(size(BasinID2));245+zeros(size(BasinID2));585+zeros(size(BasinID2))];

%Ayr = max(0,[Ayr_126;Ayr_245;Ayr_585])./1e6; %in km2

[Ayr_126 ,Ayr_245, Ayr_585] = deal(reshape(Ayr_126',[8*length(BasinID2) 1])./1e6,reshape(Ayr_245',[8*length(BasinID2) 1])./1e6,reshape(Ayr_585',[8*length(BasinID2) 1])./1e6);

id = repmat(BasinID2',8,1); id = id(:);
idx2 = repmat(idx',8,1); idx2 = idx2(:);
yrs = repmat(t',size(BasinID2));

T = [table(id), table(yrs), table(Ayr_126 ,Ayr_245, Ayr_585)];

T.Properties.VariableNames = (["BasinID2";"yrs";"SSP126";"SSP245";"SSP585"]); %num2str([2015 DeltaSLR_SSPt]')]);
writetable(T(idx2,:),'GlobalDelta_DeltaAreaWeb.csv')

%% delta water and sediment fluxes


T = table(BasinID2,Discharge_prist, Discharge_dist, Discharge_tide, QRiver_prist, QRiver_dist, QTide, QWave,delta_name,delta_area);
writetable(T(idx,:),'GlobalDelta_FluxesWeb.csv')








%{
%adapt to show cumulative
slr.DeltaSLR_series = cumsum(-slr.DeltaSLR_series,2,'reverse');
slr.DeltaSLR_series = slr.DeltaSLR_series-slr.DeltaSLR_series(:,108); %put 2007 as reference year
slr.DeltaSLR_RCP26_series = cumsum(slr.DeltaSLR_RCP26_series,2);
slr.DeltaSLR_RCP45_series = cumsum(slr.DeltaSLR_RCP45_series,2);
slr.DeltaSLR_RCP85_series = cumsum(slr.DeltaSLR_RCP85_series,2);

T = [table(BasinID2) array2table(slr.DeltaSLR_series)];
T.Properties.VariableNames = (["BasinID2";num2str(slr.DeltaSLR_time')]);
writetable(T(res.idx,:),'GlobalDelta_DeltaSLR.csv')

T = [table(BasinID2) array2table(slr.DeltaSLR_RCP26_series)];
T.Properties.VariableNames = (["BasinID2";num2str(slr.DeltaSLR_RCP_time')]);
writetable(T(res.idx,:),'GlobalDelta_DeltaSLR_RCP26.csv')

T = [table(BasinID2) array2table(slr.DeltaSLR_RCP45_series)];
T.Properties.VariableNames = (["BasinID2";num2str(slr.DeltaSLR_RCP_time')]);
writetable(T(res.idx,:),'GlobalDelta_DeltaSLR_RCP45.csv')

T = [table(BasinID2) array2table(slr.DeltaSLR_RCP85_series)];
T.Properties.VariableNames = (["BasinID2";num2str(slr.DeltaSLR_RCP_time')]);
writetable(T(res.idx,:),'GlobalDelta_DeltaSLR_RCP85.csv')
%}




%show 
%{
idx=5390;
plot(slr.DeltaSLR_time,slr.DeltaSLR_series(idx,:))
hold on
plot(slr.DeltaSLR_RCP_time,slr.DeltaSLR_RCP26_series(idx,:))
plot(slr.DeltaSLR_RCP_time,slr.DeltaSLR_RCP85_series(idx,:))
%}




