%% get QWave
%load global map of directional wave energy spectra
load([dropbox filesep 'WorldDeltas\QWave\DirMapGlo.mat'])
hs = flipud(hs); hs(hs<0.01) = nan;
energy = flipud(energy);
tp = flipud(tp);
AngArray = linspace(-0.5*pi,0.5*pi,180);
g = 9.81;
k = 5.3e-6*1050*(g^1.5)*(0.5^1.2)*(sqrt(g*0.78)/(2*pi))^0.2;
sedtrans = 2650 * (1-0.4) * k .*  (cos(AngArray).^1.2) .* sin(AngArray);

load('GlobalDeltaData.mat')
QWave = zeros(length(MouthLat),1); Hs = QWave;
res = 2;
remfun = @(lon) (rem(res*360-1+lon,res*360)+1);

wave_lat = floor(res*(91+MouthLat));
wave_lon = floor(remfun(2+res*(MouthLon)));

%figure, imagesc(hs), hold on, scatter(wave_lon,wave_lat,'r')
sbox = 1; mid_idx = ceil(((2*sbox+1)^2)/2);

for ii=1:numel(MouthLat),
    
    %find smallest but 
    [h,x,y] = max2d(-1*hs(wave_lat(ii)+(-sbox:sbox),remfun(wave_lon(ii)+(-sbox:sbox))));
    
    Hs(ii) = -1*h;
    
    if isnan(h),
        wave_lat(ii) = nan;
        wave_lon(ii) = nan;
        QWave(ii) = nan;
        
    else,
        wave_lat(ii) = wave_lat(ii)-sbox-1+x;
        wave_lon(ii) = remfun(wave_lon(ii)-sbox-1+y);
        
        sedconv = conv(fliplr(tp(wave_lat(ii),wave_lon(ii)).^0.2.*squeeze(energy(wave_lat(ii),wave_lon(ii),:))),sedtrans,'full'); % sedconv in kg/s
        QWave(ii) = (max(sedconv)-min(sedconv)).^1.2;
    end

    
end

wave_lat = floor(res*(90+MouthLat));
wave_lon = floor(remfun(1+res*(MouthLon)));

wave_lat = wave_lat/res-90;
wave_lon = wave_lon/res;

%figure, scatter(wave_lon,wave_lat,10,QWave==0,'filled')
%[~,b] = max(energy,[],3);
%imagesc(b), axis xy

save GlobalDeltaData Hs QWave -append 