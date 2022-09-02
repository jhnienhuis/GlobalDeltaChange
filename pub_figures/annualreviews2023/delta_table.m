load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','BasinID2','delta_name','MouthLon','MouthLat','QTide','QWave','QRiver');

list = ["MacKenzie";"Mississippi";"Mekong";"Parana";"Nile";"Ganges-Brahmaputra";"Dvina";"Niger";"Yellow";"Waipaoa"];


[~,idx] = ismember(list,delta_name)

%% global map
ax = worldmap('World');
land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])

scatterm(ax,MouthLat(idx),MouthLon(idx),'r','filled')
sortrows(table(list,MouthLat(idx),wrapTo180(MouthLon(idx))))

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



