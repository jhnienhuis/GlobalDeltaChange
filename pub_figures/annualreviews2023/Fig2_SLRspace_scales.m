%https://www.seanoe.org/data/00637/74862/ finally!
figure
h = ncread('D:\Drive\github\GlobalDeltaChange\pub_figures\annualreviews2023\prandi_2021.nc','trend');
%gia = ncread('D:\Drive\2022 DeltaSLR_AnnualRev\Data\gia_mean.nc','rslrate');

slrate = h; %+(1000*gia(2:2:end,2:2:end)');
slrate(slrate>999) = nan;
slrate(slrate<-20) = nan;
lim = [-10 10];
%map = flipud(cbrewer('div', 'RdYlBu', 64));

%map = [flipud(cbrewer('seq','Blues',9)); cbrewer('seq','YlOrBr',45)];

%map = [flipud(cbrewer('seq','YlGnBu',9)); cbrewer('seq','YlOrRd',20)];
map = [flipud(cbrewer('seq','Blues',270)); ones(60,3); cbrewer('seq','YlOrRd',270)];

%map(6,:) = [1 1 1];
%map(64,:) = [1 1 1];
blub = gray2ind(mat2gray(slrate,lim),length(map));

m = axesm('MapProjection','eckert4','MapLonLimit', [-180 180], 'MapLatLimit',[-72.5 80.5]);
g = geoshow(blub,map,[0.5 89 0]);
land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(m, land, 'FaceColor', [0.9 0.9 0.9])
colormap(m,map);
colorbar
set(gca, 'FontSize', 8,'FontName','Myriad Pro')
vecrast(gcf, 'Fig2a_SLRspace_scales', 300, 'bottom', 'pdf');
%figure
%imagesc(slrate,[-2 10]), axis xy, colormap(map)

%%
%add subsidence map

[A,R] = readgeoraster('D:\Drive\2017\2017 SubsidenceMap GSA\SubsidenceMapGeoTiff.tif', 'OutputType','double');
shorelines = shaperead('D:\OneDrive - Universiteit Utrecht\WorldCoastline\GSHHS_shp\h\GSHHS_h_L1.shp','UseGeoCoords', true,'BoundingBox',[lonlim(1) latlim(1);lonlim(2) latlim(2)]);
s = kml2struct('D:\Drive\Students\Joey ODell\LeveedAreas.kml');

close all
lim = [-10 10];
figure
latlim = R.LatitudeLimits;
lonlim = R.LongitudeLimits;
m = axesm('MapProjection','eckert4','MapLonLimit', lonlim, 'MapLatLimit',latlim);
hold on
geoshow(A,R,'DisplayType','surface');
%map = flipud(cbrewer('div', 'RdYlBu', 64));
%map = [flipud(cbrewer('seq','Blues',9)); ones(9,3); cbrewer('seq','YlOrRd',27)];
map = [flipud(cbrewer('seq','Blues',270)); ones(60,3); cbrewer('seq','YlOrRd',270)];
colormap(gca,map);
set(gca,'CLim',lim)
geoshow(m, s, 'FaceColor', [0.9 0.9 0.9])
geoshow([shorelines.Lat],[shorelines.Lon])
colorbar
vecrast(gcf, 'Fig2b_SLRspace_scales', 300, 'bottom', 'pdf');






