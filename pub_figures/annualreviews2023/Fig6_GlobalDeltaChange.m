clr
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData','DeltaSLR_SSPt','DeltaSLR_SSP126','DeltaSLR_SSP245','DeltaSLR_SSP585','DeltaSLR')
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaProfile','s','r','bed_h')
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaArea.mat','delta_width','delta_area','src')
load('D:\Drive\github\GlobalDeltaChange\GlobalDeltaData.mat','QRiver_dist','MouthLon','MouthLat','QTide','QWave','delta_name');
load('D:\Drive\github\GlobalDeltaChange\land_area_change\GlobalDeltaData_AreaChange.mat','net_aqua');


addpath('D:\Drive\github\GlobalDeltaSeaLevel')

m = {'Aobs','DeltaSLR_SSP126','DeltaSLR_SSP245','DeltaSLR_SSP585'}; %1985_2015

ff = 365*24*3600/1600;
fr = 0.9;

%remove caspian sea and other non deltas and deltas > 10km2
idx = (~isnan(bed_h) & bed_h<0) & ~(MouthLon>46 & MouthLon<56 & MouthLat>34 & MouthLat<48) & delta_area>1e7 &src==1;


for jj=2:length(m),
    slr = eval(m{jj});
    dA = (ff.*QRiver_dist.*fr - (delta_width.*(s-r).*0.5.*(slr)))./-bed_h;
    
    %dA = func_delta_areachange(ff.*QRiver_dist,slr,s,r,delta_width,bed_h,fr);
    dA(~idx) = nan;
    
    %Afut(jj,1:size(slr,2),:) = dA';
    out.(m{jj}) = dA';
end
idx = squeeze(~isnan(out.DeltaSLR_SSP126(1,:)));
t = [1985 2015 2050 2100 2200 2300];
out.Aobs = net_aqua'.*1e6;




%% Map
tiledlayout('flow')


land = shaperead('landareas.shp', 'UseGeoCoords', true);
[Lat, Lon] = reducem([land(:).Lat]', [land(:).Lon]',0.2);
land = geoshape(Lat,Lon,'Geometry','Polygon');

MouthLon(MouthLon>180) = MouthLon(MouthLon>180) - 360;

for ii=1:4,
    nexttile
%axes(a(ii))

%rivers = shaperead('worldrivers', 'UseGeoCoords', true);
title(m{ii},'interpreter','none')
%setm(gca,'Origin',[0 180])
axesm('MapProjection','pcarree','MapLonLimit', [-180 180],'MapLatLimit',[-60 90])
axis tight
hold on
geoshow(land, 'FaceColor', [1 0.92 0.8]);
x = out.(m{ii});
if ii>1,
   x = x(2,:); 
end
scatterm(MouthLat(idx &x>0),MouthLon(idx &x>0),5*abs(x(idx &x>0))./1e6,[0 0.8 0],'filled')
scatterm(MouthLat(idx &x<0),MouthLon(idx &x<0),5*abs(x(idx &x<0))./1e6,[0.5 0 0],'filled')

colormap([0.5 0 0;0 0.8 0]);

if ii==3,
   scatterm([-20 -10 0],[-170 -170 -170],15*abs([1 10 50]),sign([1 10 50]),'filled','MarkerEdgeColor','k')
end
   
end


set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 18.3, 10]);
set(gca, 'FontSize', 8,'FontName','Helvetica')

%% Timeseries
nexttile
for ii=2:4,
   hold on,
   Ayr = [zeros(1,length(MouthLat)); cumsum([out.Aobs; out.(m{ii})].*diff(t)',1)];
   plot(t,(sum(delta_area(idx))+sum(Ayr(:,idx),2))./1e6,'-o')
end
xlabel('Time (yr)')
ylabel('Global delta area (km2)')
set(gcf,'Renderer','Painters')
set(gca,'YGrid','on')


%% morphologic change

nexttile

QRiver_future = max(0,(QRiver_dist - (1600.*delta_area.*sum(DeltaSLR_SSP245.*diff([2000, DeltaSLR_SSPt])./300,2)/365/24/3600)));
QRiver_now = max(0,(QRiver_dist - (1600.*delta_area.*DeltaSLR./365/24/3600)));

s_now = QRiver_now+QWave+QTide;
[Rriver,~,Rtide] = deal(QRiver_now./s_now,QWave./s_now,QTide./s_now);
y0 = Rriver*sqrt(3)/2;
x0 = Rtide + y0*sqrt(3)/3;

hold on,

s_future = QRiver_future+QWave+QTide;
[Rriver,~,Rtide] = deal(QRiver_future./s_future,QWave./s_future,QTide./s_future);
y1 = Rriver*sqrt(3)/2;
x1 = Rtide + y1*sqrt(3)/3;

x_change = x1-x0;
y_change = y1-y0;

%quiver(gca,x0(idx(1:3:end)),y0(idx(1:3:end)),x_change(idx(1:3:end)),y_change(idx(1:3:end)),0,'MaxHeadSize',0.04,'Color',[0.7 0.7 0.7],'LineWidth',0.1)

ii = [0 logspace(7,10,62) inf];
blub = flipud(cbrewer('div', 'RdYlBu', 64));

colormap(flipud(cbrewer('div', 'RdYlBu', 64)))
hold on

for k= [2:64],
    iplot = (delta_area<ii(k) & delta_area>ii(k-1))';
    quiver(gca,x0(iplot & idx),y0(iplot & idx),x_change(iplot & idx),y_change(iplot & idx),0,'Color',blub(k,:),'LineWidth',max(1.5,k/30),'MaxHeadSize',0.1)
    %scatter(gca,x0(iplot & ~change),y0(iplot & ~change),k/10,blub(k,:),'filled')
end

list = ["Danube","Niger"];
[~,idx2] = ismember(list,delta_name);
text(x0(idx2),y0(idx2),list)
quiver(gca,x0(idx2),y0(idx2),x_change(idx2),y_change(idx2),0,'Color','k','LineWidth',2,'MaxHeadSize',0.1)

plot([0,0.5,1,0],[0,sqrt(3)/2,0,0]), set(gca,'DataAspectRatio',[1 1 1])

saveas(gcf,'Fig6_GlobalDeltaChange.svg')