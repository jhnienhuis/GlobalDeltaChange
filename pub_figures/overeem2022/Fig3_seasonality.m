function Fig3_seasonality

%this is from: https://github.com/jhnienhuis/GlobalDeltaChange
d = load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat');

%river ice filename
f_riverice = "Fig3_Breakup_4Rivers1985_2015.csv";
f_seaice = "Fig3_SeaiceOpenWaterSTATS_selectedeltas.xlsx";

opts = delimitedTextImportOptions("NumVariables", 4);
opts.DataLines = [2, Inf];
opts.VariableNames = ["DeltaName", "BreakupMin", "BreakupMean", "BreakupMax"];
opts.VariableTypes = ["string", "double", "double", "double"];
breakup = readtable(f_riverice, opts);

[seaice,seaice_name] = xlsread(f_seaice,'A6:D9');

figure
set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 18.3, 10]);

%% colville
subplot(2,2,4)
d_name = 'Colville';
de = find(d.delta_name==d_name);
out.t_day = (1:366)';
out.QTide_day = ones(366,1).*d.QTide(de);
[out.QWave_day, out.energy_day] = get_hs_timeseries(d.wave_lat(de),d.wave_lon(de));


[out.Qs_day, out.discharge_day] = get_qs_timeseries(d.MouthLat(de),d.MouthLon(de),d.BasinArea(de));
plot(out.t_day,[out.Qs_day out.QWave_day out.QTide_day])
hold on
opts = delimitedTextImportOptions("NumVariables", 2);
opts.VariableTypes = ["double", "double"];
ColvilleQs = table2array(readtable("Fig3_ColvilleQs_arnborg1962.csv", opts));
plot(ColvilleQs(:,1),max(0,ColvilleQs(:,2)))


plot(table2array(breakup(strcmp(breakup.DeltaName, d_name),3)).*ones(2,1),get(gca,'YLim'),'--k')
plot(table2array(breakup(strcmp(breakup.DeltaName, d_name),[2 4])),mean(get(gca,'YLim')).*ones(2,1),'--k')

plot(seaice(strcmp(seaice_name,d_name),2).*ones(2,1),get(gca,'YLim'),'-b')
plot(seaice(strcmp(seaice_name,d_name),3).*ones(2,1),get(gca,'YLim'),'-b')

title(d.delta_name(de))
datetick('x','mmm')
ylabel('Sediment flux (kg/s)')
%% lena
subplot(2,2,1)
d_name = 'Lena';
de = find(d.delta_name==d_name);
out.t_day = (1:366)';
out.QTide_day = ones(366,1).*d.QTide(de);
[out.QWave_day, out.energy_day] = get_hs_timeseries(d.wave_lat(de),d.wave_lon(de));


[out.Qs_day, out.discharge_day] = get_qs_timeseries(d.MouthLat(de),d.MouthLon(de),d.BasinArea(de));
l = plot(out.t_day,[out.Qs_day out.QWave_day out.QTide_day]);
hold on

opts = delimitedTextImportOptions("NumVariables", 2);
opts.VariableTypes = ["double", "double"];
LenaQs = table2array(readtable("Fig3_LenaMonthlyObservedQs_fromHolmes.csv", opts));
plot(LenaQs(:,1),LenaQs(:,2))


plot(table2array(breakup(strcmp(breakup.DeltaName, d_name),3)).*ones(2,1),get(gca,'YLim'),'--k')
plot(table2array(breakup(strcmp(breakup.DeltaName, d_name),[2 4])),mean(get(gca,'YLim')).*ones(2,1),'--k')

plot(seaice(strcmp(seaice_name,d_name),2).*ones(2,1),get(gca,'YLim'),'-b')
plot(seaice(strcmp(seaice_name,d_name),3).*ones(2,1),get(gca,'YLim'),'-b')

title(d.delta_name(de))
legend(l,{'QRiver','QWave','QTide'})
datetick('x','mmm')
ylabel('Sediment flux (kg/s)')

%% yukon
subplot(2,2,2)
d_name = 'Yukon';
de = find(d.delta_name==d_name);
out.t_day = (1:366)';
out.QTide_day = ones(366,1).*d.QTide(de);
[out.QWave_day, out.energy_day] = get_hs_timeseries(d.wave_lat(de),d.wave_lon(de));


[out.Qs_day, out.discharge_day] = get_qs_timeseries(d.MouthLat(de),d.MouthLon(de),d.BasinArea(de));
plot(out.t_day,[out.Qs_day out.QWave_day out.QTide_day])
hold on

opts = delimitedTextImportOptions("NumVariables", 2);
opts.VariableTypes = ["double", "double"];
YukonQs = table2array(readtable("Fig3_YukonMonthlyObservedQs_fromHolmes.csv", opts));
plot(YukonQs(:,1),YukonQs(:,2))

plot(table2array(breakup(strcmp(breakup.DeltaName, d_name),3)).*ones(2,1),get(gca,'YLim'),'--k')
plot(table2array(breakup(strcmp(breakup.DeltaName, d_name),[2 4])),mean(get(gca,'YLim')).*ones(2,1),'--k')

plot(seaice(strcmp(seaice_name,d_name),2).*ones(2,1),get(gca,'YLim'),'-b')
plot(seaice(strcmp(seaice_name,d_name),3).*ones(2,1),get(gca,'YLim'),'-b')

title(d.delta_name(de))
datetick('x','mmm')
ylabel('Sediment flux (kg/s)')

%% jago
d.delta_name(10728) = 'Jago';
subplot(2,2,3)
d_name = 'Jago';
de = find(d.delta_name==d_name);
out.t_day = (1:366)';
out.QTide_day = ones(366,1).*d.QTide(de);
[out.QWave_day, out.energy_day] = get_hs_timeseries(d.wave_lat(de),d.wave_lon(de));

[out.Qs_day, out.discharge_day] = get_qs_timeseries(d.MouthLat(de),d.MouthLon(de),d.BasinArea(de));
plot(out.t_day,[out.Qs_day out.QWave_day out.QTide_day])
hold on
plot(table2array(breakup(strcmp(breakup.DeltaName, d_name),3)).*ones(2,1),get(gca,'YLim'),'--k')
plot(table2array(breakup(strcmp(breakup.DeltaName, d_name),[2 4])),mean(get(gca,'YLim')).*ones(2,1),'--k')

plot(seaice(strcmp(seaice_name,d_name),2).*ones(2,1),get(gca,'YLim'),'-b')
plot(seaice(strcmp(seaice_name,d_name),3).*ones(2,1),get(gca,'YLim'),'-b')

title(d.delta_name(de))
datetick('x','mmm')
ylabel('Sediment flux (kg/s)')

set(findall(gcf,'type','Axes'), 'FontSize', 7)
set(findall(gcf,'Type','text'), 'FontSize', 7)

saveas(gcf,'Fig3_seasonality.svg')

end

function [qs_day, discharge_day] = get_qs_timeseries(rm_lat,rm_lon,basinarea)

rm_lon(rm_lon>180) = rm_lon(rm_lon>180)-360;

%these files are large, let me know if you're interested and i will share them
a = geotiffread('D:\OneDrive - Universiteit Utrecht\WBMSed\bqart_sed-flow_dir\bqart_a.tif');
a(a<0) = nan; a = flipud(a);
d = matfile(['D:\OneDrive - Universiteit Utrecht\WBMSed\qs_timeseries']);
wbm_lat = d.rm_lat;
wbm_lon = d.rm_lon;
lon_grid = d.lon_grid;
lat_grid = d.lat_grid;
wbm_basinarea = a(sub2ind(size(a),wbm_lat,wbm_lon));

delta_coor = rm_lat+1i.*rm_lon;
wbm_coor = lat_grid(wbm_lat)+1i.*lon_grid(wbm_lon);
blub = abs(wbm_coor-delta_coor) +...
    (abs((wbm_basinarea-basinarea)./basinarea.*log10(basinarea)))+...
    (abs(wbm_coor-delta_coor)>8).*10;

[~,idx] = min(blub);
fac = wbm_basinarea(idx)./basinarea;
%mackenzie = [459,1595], colville = [296,1605], lena=[3035,1632]

sed_series = d.sed_series(:,idx).*fac;
discharge_series = d.discharge_series(:,idx).*fac;

t = d.t;
[yr_t,~] = datevec(t);

t_day = int64(t-datenum(yr_t,1,1))+1;
t_day(t_day>366) = 366;

qs_day = double(accumarray(t_day,sed_series,[366 1],@mean));
discharge_day = double(accumarray(t_day,discharge_series,[366 1],@mean));

end

function [QWave_day, energy_day] = get_hs_timeseries(wave_lat,wave_lon)
%these files are large, let me know if you're interested and i will share them
d = load('D:\OneDrive - Universiteit Utrecht\Waves\Software\DirMapCoast.mat','grd_lon','grd_lat','hs');

idx = find(abs(wave_lat-d.grd_lat)<1.5 & abs(wave_lon-d.grd_lon)<1.5);

if isempty(idx),
    [~,idx] = min((wave_lat-d.grd_lat').^2+(wave_lon-d.grd_lon').^2,[],2);
else,
    [~,idx2] = max(d.hs(1,idx)); idx = idx(idx2);
end
%these files are large, let me know if you're interested and i will share them
f = ['D:\OneDrive - Universiteit Utrecht\Waves\WaveWatch_timeseries\WaveWatch_timeseries_' num2str(idx) '.nc'];

t = ncread(f,'time');
dp = floor(double(ncread(f,'dp'))./100);
dp(dp==360)=0;
hs = double(ncread(f,'hs'))./1000; 
tp = double(ncread(f,'tp'))./1000;
[y,~] = datevec(double(t));
t_day = int64(t-datenum(y,1,1))+1;
t_day(t_day>366) = 366;

%hs_day = accumarray(t_day,hs,[],@mean);
xx = int64(dp.*1000)+t_day;

energy_day = accumarray(t_day,(1/8).*1030.*9.81.*(hs.^2),[366 1],@sum)./length(hs).*365.24;
energy = accumarray(xx,hs.^2.4.*tp.^0.2,[359366 1],@sum)./length(dp).*365.24;

[idx_day,idx_dp] = meshgrid(1:366,0:359);
energy = energy(int64(idx_dp.*1000+idx_day));

%hold on,
%plot(sum(energy,2)),
   
AngArray = linspace(-0.5*pi,0.5*pi,180);
g = 9.81;
k = 5.3e-6*1050*(g^1.5)*(0.5^1.2)*(sqrt(g*0.78)/(2*pi))^0.2;
sedtrans = 2650 * 0.6 * k .*  (cos(AngArray).^1.2) .* sin(AngArray);

QWave_day = zeros(size(energy_day));
for ii=1:366,
    sedconv = conv(energy(:,ii),sedtrans,'full');
    QWave_day(ii) = (max(sedconv)-min(sedconv));
end

end