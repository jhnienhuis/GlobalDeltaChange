%function get_QWave_future

%% cowclip
f = 'D:\OneDrive - Universiteit Utrecht\Waves\COWCLIP\';
ti = {'historical','rcp45','rcp85'};
yr = {'197901-200412','208101-209912','208101-209912'};

%ncdisp([f 'Dm_glob_CSIRO_ACCESS1-0_historical_r1i1p1_mon_197901-200412.nc'])
lat_cow = ncread([f 'Hs_glob_CSIRO_ACCESS1-0_historical_r1i1p1_mon_197901-200412.nc'],'latitude');
lon_cow = ncread([f 'Hs_glob_CSIRO_ACCESS1-0_historical_r1i1p1_mon_197901-200412.nc'],'longitude');
[lat_cow, lon_cow] = meshgrid(lat_cow,lon_cow);
lon_cow(lon_cow==0) = 360;

Hs_cow = zeros([3 size(lon_cow)]);
Tp_cow = zeros([3 size(lon_cow)]);
for ii=1:3,
    Hs_cow_h = ncread([f 'Hs_glob_CSIRO_ACCESS1-0_' ti{ii} '_r1i1p1_mon_' yr{ii} '.nc'],'Hs_avg');
    Hs_cow(ii,:,:) = sqrt(mean(Hs_cow_h.^2,3));
    Tp_cow_h = ncread([f 'Tm_glob_CSIRO_ACCESS1-0_' ti{ii} '_r1i1p1_mon_' yr{ii} '.nc'],'Tm_avg');
    Tp_cow(ii,:,:) = mean(Tp_cow_h,3);
end

Hs_cow_h(:,lat_cow(1,:)>60,1) = nan; %replace w/ arctic

x = find(isnan(Hs_cow_h(:,:,1)));
lat_cow(x) = [];
lon_cow(x) = [];
Hs_cow(:,x) = [];
Tp_cow(:,x) = [];

%% casas-prat future arctic
f = 'D:\OneDrive - Universiteit Utrecht\Waves\Casas_Prat_futurearctic\Arctic_monthly_stat.nc';
lat_arc = ncread(f,'lat'); lon_arc = ncread(f,'lon'); Hs_arc = ncread(f,'hs'); Tp_arc = ncread(f,'tp');
Hs_arc = squeeze(sqrt(mean(Hs_arc(:,6,:,:).^2,3)))';
Tp_arc = squeeze(mean(Tp_arc(:,6,:,:),3))';

%% merge
latx = [lat_cow lat_arc'];
lonx = [lon_cow lon_arc'];
idxs = [ones(1,length(lat_cow)) 2*ones(1,length(lat_arc))];
Hs_fut = [Hs_cow  [Hs_arc(1,:); nan(size(Hs_arc(1,:))) ;Hs_arc(2,:)]]; 
Tp_fut = [Tp_cow  [Tp_arc(1,:); nan(size(Tp_arc(1,:))) ;Tp_arc(2,:)]]; 
%% matching w/ globaldelta and doing statistics
load('GlobalDeltaData.mat','wave_lat','wave_lon','QWave','Hs','Tp','Continent');

Hs_h = zeros(size(Hs)); Hs_rcp45 = zeros(size(Hs)); Hs_rcp85 = zeros(size(Hs));
Tp_h = zeros(size(Hs)); Tp_rcp45 = zeros(size(Hs)); Tp_rcp85 = zeros(size(Hs));
d = zeros(size(Hs)); Hs_src = zeros(size(Hs));

for ii=1:length(wave_lat),
    if Hs(ii)==0, continue, %far away etc,
    end
    [d(ii),idx] = min((latx-wave_lat(ii)).^2+(lonx-wave_lon(ii)).^2);
    
    Hs_src(ii) = idxs(idx);
    
    Hs_h(ii) = Hs_fut(1,idx);
    Hs_rcp45(ii) = Hs_fut(2,idx);
    Hs_rcp85(ii) = Hs_fut(3,idx);
    
    Tp_h(ii) = Tp_fut(1,idx);
    Tp_rcp45(ii) = Tp_fut(2,idx);
    Tp_rcp85(ii) = Tp_fut(3,idx);
    
end

%correct for differences between new datasets and wavewatch
Hs_rcp45 = Hs_rcp45.*Hs./Hs_h;
Hs_rcp85 = Hs_rcp85.*Hs./Hs_h;

Tp_rcp45 = Tp_rcp45.*Tp./Tp_h;
Tp_rcp85 = Tp_rcp85.*Tp./Tp_h;

g = 9.81;
k = 5.3e-6*1050*(g^1.5)*(0.5^1.2)*(sqrt(g*0.78)/(2*pi))^0.2;
sedtrans = 2650 * (1-0.4) * k .*  0.46*2;
wave_multiplier = 3; %winnowing of fine sediments by waves
QWave_h = wave_multiplier.*(Tp_h.^0.2.*Hs_h.^2.*sedtrans).^1.2;
QWave_rcp85 = (wave_multiplier.*(Tp_rcp85.^0.2.*Hs_rcp85.^2.*sedtrans).^1.2).*QWave./QWave_h;
QWave_rcp45 = (wave_multiplier.*(Tp_rcp45.^0.2.*Hs_rcp45.^2.*sedtrans).^1.2).*QWave./QWave_h;


%% arctic map


land = shaperead('landareas.shp', 'UseGeoCoords', true);
rivers = shaperead('worldrivers', 'UseGeoCoords', true);
axesm('eqdazim','MapLatLimit',[40 90])
geoshow(land, 'FaceColor', [1 0.92 0.8]);
geoshow(rivers,'Color','blue');
hold on

scatterm(wave_lat(Continent==8),wave_lon(Continent==8),30*(Hs_rcp85(Continent==8)./Hs(Continent==8)),(Hs_rcp85(Continent==8)./Hs(Continent==8)),'filled','MarkerEdgeColor','k')
c = colorbar;
c.Label.String = 'Q_{Wave},RCP8.5 (2081-2100) / Q_{Wave},(1979-2009)';