load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaArea.mat','w','delta_area')
load('D:\Drive\github\GlobalDeltaChange\GlobalDeltaData.mat','QRiver_dist','MouthLon','MouthLat','QTide','QWave','delta_name');
load('D:\Drive\github\GlobalDeltaChange\land_area_change\GlobalDeltaData_AreaChange.mat','net_aqua');

%remove caspian sea and other non deltas and deltas > 10km2
idx= ~(MouthLon>46 & MouthLon<56 & MouthLat>34 & MouthLat<48) & delta_area>1e7;


% global map
ax = worldmap('World');
land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])

scatterm(ax,MouthLat(idx),MouthLon(idx),'r','filled')
%sortrows(table(list,MouthLat(idx),wrapTo180(MouthLon(idx))))%

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       saveas(gcf,'Fig3.svg')

%% galloway
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelResponse.mat')
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaArea.mat')




[~,mor] = max([QRiver_dist(idx),QTide(idx),QWave(idx)],[],2);

s = QRiver_dist(idx)+QWave(idx)+QTide(idx);
[Rriver,Rwave,Rtide] = deal(QRiver_dist(idx)./s,QWave(idx)./s,QTide(idx)./s);

y0 = Rriver*sqrt(3)/2;
x0 = Rtide + y0*sqrt(3)/3;

scatter(x0,y0,sqrt(t.Qsettling)*2,t.Qsettling,'o');
hold on,
text(x0,y0,t.xuejiao_list)

plot([0,0.5,1,0],[0,sqrt(3)/2,0,0])



