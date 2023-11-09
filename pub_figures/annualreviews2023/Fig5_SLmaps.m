%% global map

f = [cd '\anu_rsl\rsl.7.xyz'];
gunzip([f '.gz'])
fid = fopen(f); 
SL7 = textscan(fid,'%f %f %f'); 
fclose(fid);
pause(0.1)
delete(f)



F7 = scatteredInterpolant(SL7{1},SL7{2},SL7{3});
[xx,yy] = meshgrid(0:360,-90:90);
RSL7 = F7(xx,yy)./7;
sl_dif_8 = gray2ind(mat2gray(-RSL7,[-10 10]));

f = [cd '\anu_rsl\rsl.11.xyz'];
gunzip([f '.gz'])
fid = fopen(f); 
SL8 = textscan(fid,'%f %f %f'); 
fclose(fid);
pause(0.1)
delete(f)

F8 = scatteredInterpolant(SL7{1},SL7{2},SL8{3});
[xx,yy] = meshgrid(0:360,-90:90);
sl_dif_8 = gray2ind(mat2gray((F7(xx,yy)-F8(xx,yy))./5,[-10 10]));

tiledlayout('flow')

nexttile
ax = worldmap('World');
map = flipud(cbrewer('div', 'RdYlBu', 64));
g = geoshow(sl_dif_8,map,[1 90 0],'DisplayType','image');
land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.9 0.9 0.9])

nexttile
ax = worldmap('World');
map = flipud(cbrewer('div', 'RdYlBu', 64));
g = geoshow(sl_dif_8,map,[1 90 0],'DisplayType','image');
land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.9 0.9 0.9])


nexttile
ax = worldmap('World');
map = flipud(cbrewer('div', 'RdYlBu', 64));
load('gridded-SL-IPCC-AR6.mat','SL_245LC_2200_50')

sl_dif_future = (SL_245LC_2200_50)./200;
sl_dif_future = gray2ind(mat2gray(sl_dif_future,[-10 10]));

g = geoshow(rot90(sl_dif_future),map,[1 90 0],'DisplayType','image');

land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.9 0.9 0.9])


vecrast(gcf, 'Fig3_SLmaps', 300, 'bottom', 'pdf');

%% arrows per delta
f = ['D:\Drive\2022 DeltaSLR_AnnualRev\Data\anu_rsl\rsl.7.xyz'];
gunzip([f '.gz'])
fid = fopen(f); 
SL7 = textscan(fid,'%f %f %f'); 
fclose(fid);
pause(0.1)
delete(f)
F7 = scatteredInterpolant(SL7{1},SL7{2},SL7{3});


%}


%% arrows per delta 2
load('D:\Drive\github\GlobalDeltaChange\GlobalDeltaData.mat','MouthLon','MouthLat','QRiver_prist','Continent');
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaArea.mat','delta_area')
idx = delta_area>1e7;

load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData.mat')

b = [-inf, -10:1:10 inf];
map = flipud(cbrewer('div', 'RdYlBu', length(b)+3));
map([25,20,15],:) = [];
% 7ka to now
tiledlayout('flow')
nexttile

ax = worldmap('World');
setm(ax,'mapprojection','eqdcylin');

land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.9 0.9 0.9])
hold on,


SLR = HoloceneSL(:,4)./7;

xx  = discretize(SLR,b);

SLR = max(-10,min(20,SLR));

for ii=1:length(b),
    if sum(idx & xx==ii)==0,
        continue,
    end
   %     scatterm(MouthLat(idx & xx==ii),MouthLon(idx & xx==ii),5,map(ii,:),'filled'), 
    plotm([MouthLat(idx & xx==ii) MouthLat(idx & xx==ii)+SLR(idx & xx==ii)]',[MouthLon(idx & xx==ii) MouthLon(idx & xx==ii)]','Color',map(ii,:),'LineWidth',2),
end

% past century

nexttile

ax = worldmap('World');
setm(ax,'mapprojection','eqdcylin');

land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.9 0.9 0.9])
hold on,



SLR = 1000*(DeltaSLR+DeltaSub);

xx  = discretize(SLR,b);

SLR = max(-10,min(20,SLR));

for ii=1:length(b),
    if sum(idx & xx==ii)==0,
        continue,
    end
    
  %  scatterm(MouthLat(idx & xx==ii),MouthLon(idx & xx==ii),5,map(ii,:),'filled'), 
    
    plotm([MouthLat(idx & xx==ii) MouthLat(idx & xx==ii)+SLR(idx & xx==ii)]',[MouthLon(idx & xx==ii) MouthLon(idx & xx==ii)]','Color',map(ii,:),'LineWidth',2),
end
%legend:
plotm([-70; -70+-10],[170;170],'Color',map(2,:),'LineWidth',2),
plotm([-70; -70+20],[170;170],'Color',map(end,:),'LineWidth',2),
plotm([-70; -70+10],[170;170],'Color',map(find(b>10,1),:),'LineWidth',2),

%% projections
nexttile

ax = worldmap('World');
setm(ax,'mapprojection','eqdcylin');

land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.9 0.9 0.9])
hold on,

SLR = mean(1000*(DeltaSLR_SSP245(:,1:2)),2);

xx  = discretize(SLR,b);

SLR = max(-10,min(20,SLR));

for ii=1:length(b),
    if sum(idx & xx==ii)==0,
        continue,
    end
       % scatterm(MouthLat(idx & xx==ii),MouthLon(idx & xx==ii),5,map(ii,:),'filled'), 
    plotm([MouthLat(idx & xx==ii) MouthLat(idx & xx==ii)+SLR(idx & xx==ii)]',[MouthLon(idx & xx==ii) MouthLon(idx & xx==ii)]','Color',map(ii,:),'LineWidth',2),
end

nexttile
ax = worldmap('World');
setm(ax,'mapprojection','eqdcylin');

land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.9 0.9 0.9])
hold on,

SLR = mean(1000*(DeltaSLR_SSP245(:,3:4)),2);

xx  = discretize(SLR,b);

SLR = max(-10,min(20,SLR));

for ii=1:length(b),
    if sum(idx & xx==ii)==0,
        continue,
    end
       % scatterm(MouthLat(idx & xx==ii),MouthLon(idx & xx==ii),5,map(ii,:),'filled'), 
    plotm([MouthLat(idx & xx==ii) MouthLat(idx & xx==ii)+SLR(idx & xx==ii)]',[MouthLon(idx & xx==ii) MouthLon(idx & xx==ii)]','Color',map(ii,:),'LineWidth',2),
end

set(gcf, 'Units', 'Centimeters','Renderer','painters');
set(gca, 'FontSize', 8,'FontName','Helvetica')
saveas(gcf,'Fig5_SLmaps.pdf')

