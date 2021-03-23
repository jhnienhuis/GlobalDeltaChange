function [QWave_day, energy_day,num_day_50pct] = get_Qwave_timeseries(wave_lat,wave_lon)
d = load('D:\OneDrive - Universiteit Utrecht\Waves\Software\DirMapCoast.mat','grd_lon','grd_lat','hs');

idx = find(abs(wave_lat-d.grd_lat)<1.5 & abs(wave_lon-d.grd_lon)<1.5);

if isempty(idx),
    [~,idx] = min((wave_lat-d.grd_lat').^2+(wave_lon-d.grd_lon').^2,[],2);
else,
    [~,idx2] = max(d.hs(1,idx)); idx = idx(idx2);
end



%scatter(int64(grd_lon(idx).*30),int64(grd_lat(idx).*30),50)
%hold on
%scatter(int64(wave_lon.*30),int64(wave_lat.*30),'filled')
%m = matfile('D:\OneDrive - Universiteit Utrecht\Waves\WaveWatch_timeseries\WaveWatch_timeseries.mat');

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

num_day_50pct = find(cumsum(sort(QWave_day,'descend'))>(0.5*sum(QWave_day)),1);