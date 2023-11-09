%sea level risk map
%% arrows per delta 2
load('D:\Drive\github\GlobalDeltaChange\GlobalDeltaData.mat','MouthLon','MouthLat','QRiver_prist','delta_name','BasinID2');
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaArea.mat','delta_area')

t = readtable("Nicholls_subs.csv");

[~,idx] = ismember(t{:,1},BasinID2);


slr = ncread('77853.nc','trend')';
slr_lat = ncread('77853.nc','latitude');
slr_lon = ncread('77853.nc','longitude');


MouthLon = mod(MouthLon,360);

[slr_lat_grid, slr_lon_grid] = meshgrid(slr_lat,slr_lon);
slr_lon_grid = slr_lon_grid(:);
slr_lat_grid = slr_lat_grid(:);
slr_lon_grid(slr>999 | slr<-60) = [];
slr_lat_grid(slr>999 | slr<-60) = [];
slr(slr>999 | slr<-60) = [];

k = dsearchn([slr_lon_grid,slr_lat_grid],[MouthLon(idx) MouthLat(idx)]); 


b = [-inf, -10:10 inf];
map = flipud(cbrewer('div', 'RdYlBu', length(b)+3));
map(13,:) = [1,1,1];

ax = worldmap('World');
setm(ax,'mapprojection','eqdcylin');

land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.9 0.9 0.9]);
hold on,

SUB = 1000*(t{:,3});

xx  = discretize(2.*slr(k),b);
xy  = discretize(-2.*t{:,3},b);
%SLR = max(-10,min(20,SLR));

for ii=1:length(idx),
   scatterm(MouthLat(idx(ii)),MouthLon(idx(ii)),'Color','k')
   plotm([MouthLat(idx(ii)) MouthLat(idx(ii))+(2*slr(k(ii)))]',[MouthLon(idx(ii)) MouthLon(idx(ii))]','Color',map(xx(ii),:),'LineWidth',4),
   plotm([MouthLat(idx(ii)) MouthLat(idx(ii))-(2*t{ii,3})]',[MouthLon(idx(ii)) MouthLon(idx(ii))]','Color',map(xy(ii),:),'LineWidth',4),

end
%legend:
plotm([-70; -70+-20],[170;170],'Color',map(2,:),'LineWidth',4),
plotm([-70; -70+-10],[170;170],'Color',map(7,:),'LineWidth',4),
plotm([-70; -70+20],[170;170],'Color',map(end,:),'LineWidth',4),
plotm([-70; -70+10],[170;170],'Color',map(find(b>5,1),:),'LineWidth',4),

set(gcf, 'Units', 'Centimeters','Renderer','painters');
set(gca, 'FontSize', 8,'FontName','Helvetica')
saveas(gcf,'Fig1_SLR_riskmap.pdf')