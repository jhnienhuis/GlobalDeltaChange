%% get QWave
function get_Qwave
%load global map of directional wave energy spectra
load(['D:\OneDrive - Universiteit Utrecht\Waves\Software\DirMapGlo.mat'])
hs = flipud(hs); hs(hs<0.01) = nan;
energy = flipud(energy);
tp = flipud(tp); tp(tp<0.1) = nan;
AngArray = linspace(-0.5*pi,0.5*pi,180);
g = 9.81;
k = 5.3e-6*1050*(g^1.5)*(0.5^1.2)*(sqrt(g*0.78)/(2*pi))^0.2;
sedtrans = 2650 * (1-0.4) * k .*  (cos(AngArray).^1.2) .* sin(AngArray);
wave_multiplier = 3; %winnowing of fine sediments by waves
load('GlobalDeltaData.mat')
QWave = zeros(length(MouthLat),1); Hs = QWave; Tp = Hs;
res = 2;
remfun = @(lon) (rem(res*360-1+lon,res*360)+1);

wave_lat = floor(res*(91+MouthLat));
wave_lon = floor(remfun(2+res*(MouthLon)));

%figure, imagesc(hs), hold on, scatter(wave_lon,wave_lat,'r')
sbox = 1; mid_idx = ceil(((2*sbox+1)^2)/2);

for ii=1:numel(MouthLat),
    
    %find largest but 
    [~,x,y] = max2d(hs(wave_lat(ii)+(-sbox:sbox),remfun(wave_lon(ii)+(-sbox:sbox))));
    Hs(ii) = nanmedian(hs(wave_lat(ii)+(-sbox:sbox),remfun(wave_lon(ii)+(-sbox:sbox))),'all');
    Tp(ii) = nanmedian(tp(wave_lat(ii)+(-sbox:sbox),remfun(wave_lon(ii)+(-sbox:sbox))),'all');
    
    
    if isnan(Hs(ii)),
        wave_lat(ii) = nan;
        wave_lon(ii) = nan;
        QWave(ii) = 0;
        Hs(ii) = 0;
        Tp(ii) = 0;
    else,
        wave_lat(ii) = wave_lat(ii)-sbox-1+x;
        wave_lon(ii) = remfun(wave_lon(ii)-sbox-1+y);
        
        sedconv = conv(fliplr(tp(wave_lat(ii),wave_lon(ii)).^0.2.*squeeze(energy(wave_lat(ii),wave_lon(ii),:))),sedtrans,'full'); % sedconv in kg/s
        QWave(ii) = wave_multiplier.*(max(sedconv)-min(sedconv)).^1.2;
    end

    
end

wave_lat = floor(res*(90+MouthLat));
wave_lon = floor(remfun(1+res*(MouthLon)));

wave_lat = wave_lat/res-90;
wave_lon = wave_lon/res;

%figure, scatter(wave_lon,wave_lat,10,QWave==0,'filled')
%[~,b] = max(energy,[],3);
%imagesc(b), axis xy

save GlobalDeltaData Hs Tp QWave -append 

