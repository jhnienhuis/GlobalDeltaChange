clr
load('D:\Drive\github\GlobalDeltaChange\GlobalDeltaData.mat','BasinID2','delta_name','MouthLon','MouthLat','QTide','QWave','QRiver_dist');

load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData.mat','DeltaSub','DeltaSLR','DeltaSLR_RCP45_2100');
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelResponse.mat','idx');
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaProfile.mat','r','s')
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaArea.mat','w','delta_area')

% galloway
QRiver_future = max(0,(QRiver_dist - (1600.*w.*(s-r).*(DeltaSLR_RCP45_2100)./2/365/24/3600)));

s_now = QRiver_dist+QWave+QTide;
[Rriver,~,Rtide] = deal(QRiver_dist./s_now,QWave./s_now,QTide./s_now);
y0 = Rriver*sqrt(3)/2;
x0 = Rtide + y0*sqrt(3)/3;

hold on,

s_future = QRiver_future+QWave+QTide;
[Rriver,~,Rtide] = deal(QRiver_future./s_future,QWave./s_future,QTide./s_future);
y1 = Rriver*sqrt(3)/2;
x1 = Rtide + y1*sqrt(3)/3;

x_change = x1-x0;
y_change = y1-y0;

quiver(gca,x0(idx(1:3:end)),y0(idx(1:3:end)),x_change(idx(1:3:end)),y_change(idx(1:3:end)),0,'MaxHeadSize',0.04,'Color',[0.7 0.7 0.7],'LineWidth',0.1)


list = ["MacKenzie";"Mississippi";"Mekong";"Parana";"Nile";"Ganges-Brahmaputra";"Dvina";"Niger";"Yellow";"Waipaoa"];
[~,idx] = ismember(list,delta_name);
text(x0(idx),y0(idx),list)
quiver(gca,x0(idx),y0(idx),x_change(idx),y_change(idx),0,'MaxHeadSize',0.1,'Color','k')



plot([0,0.5,1,0],[0,sqrt(3)/2,0,0]), set(gca,'DataAspectRatio',[1 1 1])

set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 10, 10]);
set(gca, 'FontSize', 8,'FontName','Myriad Pro')
saveas(gcf,'Fig5a_MorphologicChange.svg')



%% histogram change compared to delta area
out = load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelResponse.mat');
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaArea.mat','delta_area')
load('D:\Drive\github\GlobalDeltaChange\GlobalDeltaData.mat','MouthLon','MouthLat','BasinID2','delta_name');
ee = load('D:\Drive\github\GlobalDeltaChange\land_area_change\GlobalDeltaData_AreaChange.mat');
out.delta_change_obs = mean([ee.net_aqua ee.net_pekel],2).*1e6;

hold on,
histogram(out.delta_change_obs(out.idx)./delta_area(out.idx),[-0.2:0.002:0.2],'facecolor','g','facealpha',.5)
histogram(out.delta_change_RCP85_2100(out.idx)./delta_area(out.idx),[-0.2:0.002:0.2],'facecolor','r','facealpha',.5)

list = ["MacKenzie";"Mississippi";"Mekong";"Parana";"Nile";"Ganges-Brahmaputra";"Dvina";"Niger";"Yellow";"Waipaoa"];
[~,idx_small] = ismember(list,delta_name);

%scatter(out.delta_change_RCP85_2100(idx_small)./delta_area(idx_small),zeros(size(idx_small)))
scatter(delta_area(out.idx),out.delta_change_RCP85_tot(out.idx)./delta_area(out.idx))
set(gca,'YLim',[-1 1],'XScale','log')

%% Paola Map
%{
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData.mat','DeltaSub','DeltaSLR_RCP45_2100');
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelResponse.mat','idx');
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaArea.mat','delta_area');

load('D:\Drive\github\GlobalDeltaChange\GlobalDeltaData.mat','MouthLon','MouthLat','BasinID2','QRiver_dist');

%paola area (km2):
Aw = max(1,QRiver_dist.*1600.*(365*24*3600)./(DeltaSLR_RCP45_2100+DeltaSub))./1e6;

m = {'obs','RCP26_2100','RCP45_2100','RCP85_2100'}; %1985_2015
%a = tight_subplot(2,2);
land = shaperead('landareas.shp', 'UseGeoCoords', true);
[Lat, Lon] = reducem([land(:).Lat]', [land(:).Lon]',0.2);
land = geoshape(Lat,Lon,'Geometry','Polygon');

MouthLon(MouthLon>180) = MouthLon(MouthLon>180) - 360;

for ii=1,
    
%axes(a(ii))

%rivers = shaperead('worldrivers', 'UseGeoCoords', true);
title(m{ii},'interpreter','none')
%setm(gca,'Origin',[0 180])
axesm('MapProjection','pcarree','MapLonLimit', [-180 180],'MapLatLimit',[-60 90])
axis tight
hold on
geoshow(land, 'FaceColor', [1 0.92 0.8]);
scatterm(MouthLat(idx),MouthLon(idx),5*abs(Aw(idx))./1e9,Aw(idx)-delta_area(idx),'filled')
%colormap((cbrewer('div', 'RdYlGn', 64)))
%colormap([0.5 0 0;0 0.8 0]);
%set(a(ii),'CLim',[-2 2])

if ii==3,
   scatterm([-20 -10 0],[-170 -170 -170],15*abs([1 10 50]),sign([1 10 50]),'filled','MarkerEdgeColor','k')
end
   
end


set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 18.3, 10]);
set(gca, 'FontSize', 8,'FontName','Helvetica')
%saveas(gcf,'Fig5b_MorphologicChange.svg')
%}