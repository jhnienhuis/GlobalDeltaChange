%% ED Fig. 1
load('CoastalZone_int8.mat')
imagesc(cz((107*240):(240*114),(260*240):(268*240))), axis xy, axis equal, axis tight
demcmap([-1 124])

%% QWave
a = tight_subplot(2,1,[0.02,0.02],[0.09 0.01],[0.09 0.01]);
axes(a(1))
h = imagesc(grb_lon,grb_lat,log10(QWave),[0 4]); axis xy equal tight
colormap(viridis)
set(h, 'AlphaData', QWave~=0)
ylabel('Latitude'),
set(a(1),'XTickLabel',{''})
colorbar('YTick',0:4,'YTickLabel',10.^(0:4))

axes(a(2))
h = imagesc(grb_lon,grb_lat,log10(QDiff),[-4 0]); axis xy equal tight
colormap(viridis)
set(h, 'AlphaData', QWave~=0)
xlabel('Longitude'), ylabel('Latitude')
colorbar
colorbar('YTick',-4:0,'YTickLabel',10.^(-4:0))

print -dpdf -r2000 Fig_ED3_QWave.pdf


%% Qriver
load('D:\GlobalDatasets\WBMSed\Prist_Yield.mat')

%area per cell of WBMSed in m2
res_wbm = 10;
areapercell = repmat(6371.^2.*2*pi/360/res_wbm*(sin(deg2rad(-89.9:0.1:90))-sin(deg2rad(-90:0.1:89.9))),360*res_wbm,1);
areapercell = areapercell*1e6;

a = tight_subplot(3,1,[0.02,0.02],[0.09 0.01],[0.09 0.01]);
axes(a(1))
h = imagesc(0:0.1:360,-89.9:0.1:90,real(log10((WBMdis./areapercell)')),[-10 -7]); axis xy equal tight
colormap(viridis); %
set(h, 'AlphaData', WBMdis'~=0)
ylabel('Latitude'),
set(a(1),'XTickLabel',{''},'YLim',[-60 90])
h = colorbar;
ylabel(h,'Discharge yield (m s-1)')

axes(a(2))
h = imagesc(0:0.1:360,-89.9:0.1:90,real(log10(WBMsed./areapercell)'),[-10 -7]); axis xy equal tight
colormap(viridis); %flipud(cbrewer('div', 'RdYlBu', 64)))
set(h, 'AlphaData', WBMsed'~=0)
ylabel('Latitude'),
set(a(2),'XTickLabel',{''},'YLim',[-60 90])
h = colorbar;
ylabel(h,'Sediment yield (kg s-1 m-2)')
print -dpdf -r2000 Fig_ED2.pdf
figure
%axesm('MapProjection','pcarree','FLonLimit', [0 360])
load([dropbox filesep 'WorldDeltas' filesep 'GlobalDeltaData.mat'])

%h = worldmap('World');
%delete(findobj(allchild(h),'Type','text'));
land = shaperead('landareas.shp', 'UseGeoCoords', true);
%rivers = shaperead('worldrivers', 'UseGeoCoords', true);
%mapshow(land, 'FaceColor', [1 0.92 0.8]);
%setm(gca,'Origin',[0 180])
%axesm('MapProjection','pcarree','MapLonLimit', [0 360])

%geoshow(rivers,'Color','blue');

hold on, axis equal, axis tight
blub = struct2cell(land); plot(mod([blub{3,:}],360),[blub{4,:}])

y = QRiver_dist./QRiver_prist; y(y>2) = 2;
%MouthLon(MouthLon>180) = MouthLon(MouthLon>180)-360;
scatter(MouthLon(y>0.5 & y<1.5),MouthLat(y>0.5 & y<1.5),20,y(y>0.5 & y<1.5),'filled')
scatter(MouthLon(y>1.5),MouthLat(y>1.5),20,y(y>1.5),'filled')
scatter(MouthLon(y<0.5),MouthLat(y<0.5),20,y(y<0.5),'filled')

caxis([0 2])
colormap(flipud(cbrewer('div', 'RdYlBu', 64)))
colorbar('SouthOutside')
ylabel('Latitude'),
set(gca,'XTickLabel',{''},'YLim',[-60 90]), box on

print -dpdf -r2000 Fig_ED3.pdf

%% Qtide
direc = 'D:\GlobalDatasets\Tides\';
lat_tide = ncread([direc 'hf.m2_tpxo8_atlas_30c_v1.nc'],'lat_z');
lon_tide = ncread([direc 'hf.m2_tpxo8_atlas_30c_v1.nc'],'lon_z');

m2 = int32(abs(single(ncread([direc 'hf.m2_tpxo8_atlas_30c_v1.nc'],'hRe'))+1i*single(ncread([direc 'hf.m2_tpxo8_atlas_30c_v1.nc'],'hIm'))));
s2 = int32(abs(single(ncread([direc 'hf.s2_tpxo8_atlas_30c_v1.nc'],'hRe'))+1i*single(ncread([direc 'hf.s2_tpxo8_atlas_30c_v1.nc'],'hIm'))));
k1 = int32(abs(single(ncread([direc 'hf.k1_tpxo8_atlas_30c_v1.nc'],'hRe'))+1i*single(ncread([direc 'hf.k1_tpxo8_atlas_30c_v1.nc'],'hIm'))));
o1 = int32(abs(single(ncread([direc 'hf.o1_tpxo8_atlas_30c_v1.nc'],'hRe'))+1i*single(ncread([direc 'hf.o1_tpxo8_atlas_30c_v1.nc'],'hIm'))));

%hamax = reshape(max([m2(:),s2(:),k1(:),o1(:)],[],2),size(s2));
hamax = (m2+s2+o1+k1)/1.3;

hamax = double(hamax)./1000;
h = imagesc(lon_tide,lat_tide,hamax,[0 2]); axis xy equal tight
colormap(viridis)
set(h, 'AlphaData', hamax~=0)
ylabel('Latitude'),xlabel('Longitude')
h = colorbar;
ylabel(h,'average tidal amplitude (m)')

print -dpdf -r2000 Fig_ED3_Qtide.pdf
