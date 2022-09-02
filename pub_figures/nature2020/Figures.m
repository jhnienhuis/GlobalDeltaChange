%% fig. 1b
clr
load([dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep 'GlobalDeltaData'],'QRiver_dist','QRiver_prist','QWave','QTide');

QTide(QTide<=0 | isnan(QTide)) = 1;
QWave(QWave<=0 | isnan(QWave)) = 1;
QRiver_prist(QRiver_prist<=0 | isnan(QRiver_prist)) = 1;

[QRiver_prist_log,QWave_prist_log,QTide_prist_log] = DeltaLogMaker(QRiver_prist,QWave,QTide);

colormap(flipud(cbrewer('div', 'RdYlBu', 64)))
ternplot(QTide_prist_log,QRiver_prist_log,QWave_prist_log,'scatter',3+(QRiver_prist).^0.7,log10(QRiver_prist),'filled')
caxis([0,4])
set(gca,'Layer','top')
set([findall(gcf,'String','  25'); findall(gcf,'String','25')],'String','0.1')
set([findall(gcf,'String','  50'); findall(gcf,'String','50')],'String','0.5')
set([findall(gcf,'String','  75'); findall(gcf,'String','75')],'String','0.9')
h = colorbar;
set(h,'YTick',0:5,'YTickLabel',cellstr(num2str((10.^(0:5))','%1.0f'))); ylabel(h,'QRiver (kgs^{-1})')

%% fig. 1c
clr
load([dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep 'GlobalDeltaData'],'QRiver_dist','QRiver_prist','QWave','QTide','MouthLat','MouthLon');

[~,mor] = max([QWave,QRiver_prist,QTide],[],2);

axesm('Robinson','MapLonLimit', [0 360])
setm(gca,'Origin',[0 0])
land = shaperead('landareas.shp', 'UseGeoCoords', true);
rivers = shaperead('worldrivers', 'UseGeoCoords', true);
geoshow(land, 'FaceColor', [1 0.92 0.8]);
framem('on')

geoshow(rivers,'Color','blue');

scatterm(MouthLat(mor==1),MouthLon(mor==1),'filled','b')
scatterm(MouthLat(mor==2),MouthLon(mor==2),'filled','g')
scatterm(MouthLat(mor==3),MouthLon(mor==3),'filled','r')

%% fig. 2a
clr
load([dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep 'GlobalDeltaData'],'QRiver_dist','QRiver_prist','QWave','QTide','MouthLat','MouthLon','BasinID2');

[QRiver_dist_log,QWave_dist_log,QTide_dist_log] = DeltaLogMaker(QRiver_dist,QWave,QTide);
[QRiver_prist_log,QWave_prist_log,QTide_prist_log] = DeltaLogMaker(QRiver_prist,QWave,QTide);

[~,x0,y0] =ternplot(QTide_prist_log,QRiver_prist_log,QWave_prist_log,'scatter','SizeData',5,'Marker','o','MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor',[0.4 0.4 0.4]);
[~,x1,y1] =ternplot(QTide_dist_log,QRiver_dist_log,QWave_dist_log,'scatter','SizeData',5,'Marker','o','MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor',[0.4 0.4 0.4]);
close all
x_change = x1-x0;
y_change = y1-y0;
change = sqrt(x_change.^2+y_change.^2)>0.01;

figure
ii = [0 logspace(0,3,62) inf];
hold on
blub = flipud(cbrewer('div', 'RdYlBu', 64));



for k= [2:64],
    iplot = QRiver_prist<ii(k) & QRiver_prist>ii(k-1);
    quiver(x0(iplot & change),y0(iplot & change),x_change(iplot & change),y_change(iplot & change),0,'Color',blub(k,:),'LineWidth',max(1.5,k/30),'MaxHeadSize',0.1)
    scatter(x0(iplot & ~change),y0(iplot & ~change),k/10,blub(k,:),'filled')
end
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'Layer','top')
set([findall(gcf,'String','  25'); findall(gcf,'String','25')],'String','0.1')
set([findall(gcf,'String','  50'); findall(gcf,'String','50')],'String','0.5')
set([findall(gcf,'String','  75'); findall(gcf,'String','75')],'String','0.9')
colormap(flipud(cbrewer('div', 'RdYlBu', 64)));
h = colorbar;
set(h,'YTick',linspace(0,1,4),'YTickLabel',cellstr(num2str((10.^(0:3))','%1.0f'))); ylabel(h,'QRiver (kgs^{-1})')
axis off
%% fig. 2b
close all
[~,~,x] = xlsread([dropbox filesep 'WorldDeltas' filesep 'FamousDeltaData.xlsx'],'A1:D69');
[~,idx] = ismember([x{3:end,1}].*10+[x{3:end,4}],BasinID2);
text(x0(idx),y0(idx),x(3:end,2))
hold on
quiver(x0(idx),y0(idx),x_change(idx),y_change(idx),0,'MaxHeadSize',0.1)
scatter(x1(idx),y1(idx),'r','filled')
set(gca,'Layer','top')

set([findall(gcf,'String','  25'); findall(gcf,'String','25')],'String','0.1')
set([findall(gcf,'String','  50'); findall(gcf,'String','50')],'String','0.5')
set([findall(gcf,'String','  75'); findall(gcf,'String','75')],'String','0.9')

%% fig. 3a
clr
load([dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep 'GlobalDeltaData'],'QRiver_dist','QRiver_prist','QWave','QTide','MouthLat','MouthLon','BasinID2');
ee = load([dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep 'land_area_change\GlobalDeltaData_AreaChange.mat']);


land = ee.net_aqua;
change = QRiver_dist./QRiver_prist;

edges = linspace(0,2,21);
[xbin] = discretize(change,edges);
ybin = accumarray(xbin(~isnan(xbin) & QRiver_prist>10),land(~isnan(xbin) & QRiver_prist>10),[],@(x) (trimmean(x,15)));

scatter(edges(2:end)-0.05,ybin,40,'k','filled'), grid on, hold on
scatter(edges(2)-0.05,-0.02,40,'k','filled'); text(edges(2)-0.02,-0.018,num2str(round(ybin(1),1)));
scatter(change,land,1,[0.4 0.4 0.4],'filled')
set(gca,'YLim',[-0.02 0.08]), box on, set(gca,'XLim',[0 2])
plot([0 2],polyval(polyfit(change(change>edges(1)&change<edges(end)),land(change>edges(1)&change<edges(end)),1),[0 2]),'r');

xlabel('Q_{dist}/Q_{prist}')
ylabel('Land Area Change (km2/yr)')
%% fig. 3b
clr
load([dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep 'GlobalDeltaData'],'QRiver_dist','QRiver_prist','QWave','QTide','MouthLat','MouthLon','BasinID2');
ee = load([dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep 'land_area_change\GlobalDeltaData_AreaChange.mat']);

s = QWave+QTide+QRiver_dist;

fW = QWave./s;
fT = QTide./s;
fR = 1-(fW+fT);
[~,mor] = max([QWave,QRiver_prist,QTide],[],2);

fd = @(x) (sum(ee.net_aqua(x)));

delta_gain(1) = fd(fW>=0.9);
delta_gain(2) = fd(fR>=0.9);
delta_gain(3) = fd(fT>=0.9);
delta_gain(4) = fd(fW<=0.9 & mor==1);
delta_gain(5) = fd(fR<=0.9 & mor==2);
delta_gain(6) = fd(fT<=0.9 & mor==3);

fd = @(x) (mean(ee.net_aqua(x)));

delta_gain_av(1) = fd(fW>=0.9);
delta_gain_av(2) = fd(fR>=0.9);
delta_gain_av(3) = fd(fT>=0.9);
delta_gain_av(4) = fd(fW<=0.9 & mor==1);
delta_gain_av(5) = fd(fR<=0.9 & mor==2);
delta_gain_av(6) = fd(fT<=0.9 & mor==3);

table({'WW','RR','TT','W','R','T'}',delta_gain',delta_gain_av')
