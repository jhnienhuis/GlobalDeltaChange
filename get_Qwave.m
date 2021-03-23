%% get QWave
function get_Qwave
%load global map of directional wave energy spectra
load(['D:\OneDrive - Universiteit Utrecht\Waves\Software\DirMapGlo.mat'],'hs','energy','tp','wind')

load([dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep 'GlobalDeltaData.mat'],'MouthLon','MouthLat','BasinID2','delta_name','QRiver_prist');

load('D:\Dropbox\GlobalFetch\GlobalFetch.mat','FetchAll','lon','lat')

hs = flipud(hs); hs(hs<0.01) = nan;
wind = flipud(wind); wind(isnan(hs)) = nan;
energy = flipud(energy);
tp = flipud(tp); tp(isnan(hs)) = nan;
AngArray = linspace(-0.5*pi,0.5*pi,180);
g = 9.81;
k = 5.3e-6*1050*(g^1.5)*(0.5^1.2)*(sqrt(g*0.78)/(2*pi))^0.2;
sedtrans = 2650 * (1-0.4) * k .*  (cos(AngArray).^1.2) .* sin(AngArray);
wave_multiplier = 3; %winnowing of fine sediments by waves

%local waves, unlimited depth because easy, from https://www.helpdeskwater.nl/publish/pages/157091/1209433-000-hye-0013-v4-r-inputdatabaseforthebretschneiderwavecalculationsfornarrowriverareas.pdf
Hsl = @(u,F) (0.283*tanh(0.0125*(F.*g./(u.^2)).^0.42).*(u.^2)./g);
Tpl = @(u,F) (2.4*pi*tanh(0.077*(F.*g./(u.^2)).^0.25).*u./g.*1.08);


QWave = zeros(length(MouthLat),1); Hs = QWave; Tp = Hs;
res = 2;
remfun = @(lon) (rem(res*360-1+lon,res*360)+1);

wave_lat = floor(res*(91+MouthLat));
wave_lon = floor(remfun(2+res*(MouthLon)));

%figure, imagesc(hs), hold on, scatter(wave_lon,wave_lat,'r')


%get fetch
FetchDeltas = get_Qwave_fetch(MouthLon,MouthLat,lat,lon,FetchAll,QRiver_prist,delta_name);

for ii=1:numel(MouthLat),

    sbox = max(1,ceil(log10(QRiver_prist(ii)))-2); 
    
    %find fetch, convert to nautical degrees
    fe = FetchDeltas(ii,[271:360,1:270]);
    
    [~,x,y] = max2d(hs(wave_lat(ii)+(-sbox:sbox),remfun(wave_lon(ii)+(-sbox:sbox))));
    
    wave_lat(ii) = wave_lat(ii)-sbox-1+x;
    wave_lon(ii) = remfun(wave_lon(ii)-sbox-1+y);
    
    if sum(fe==127)==0 || isnan(hs(wave_lat(ii),wave_lon(ii))), %include rare cases where no wave data was found nearby, assume 10m/s wind.
                
        %only local waves
        fe_mean = mean(max(0,fe));
        
        w = wind(wave_lat(ii),wave_lon(ii));
        w(isnan(w)) = 10;
                 
        Hs(ii) = Hsl(w./4,fe_mean.*1e3);
        Tp(ii) = Tpl(w./4,fe_mean.*1e3);
        sedconv = conv(Tp(ii).^0.2.*Hs(ii).^2.4,sedtrans,'full'); % sedconv in kg/s
        QWave(ii) = wave_multiplier.*(max(sedconv)-min(sedconv));
        
    else,

        
        %only get energy directed at delta.
        energy_in = squeeze(energy(wave_lat(ii),wave_lon(ii),:)).*(fe==127)';
        
        Hs(ii) = sqrt(sum(energy_in));
        Tp(ii) = nanmedian(tp(wave_lat(ii)+(-sbox:sbox),remfun(wave_lon(ii)+(-sbox:sbox))),'all');
                
        sedconv = conv(fliplr(tp(wave_lat(ii),wave_lon(ii)).^0.2.*energy_in),sedtrans,'full'); % sedconv in kg/s
        QWave(ii) = wave_multiplier.*(max(sedconv)-min(sedconv)).^1.2;
    end
       
end
%{
wave_lat = floor(res*(90+MouthLat));
wave_lon = floor(remfun(1+res*(MouthLon)));

wave_lat = wave_lat/res-90;
wave_lon = wave_lon/res;

%figure, scatter(wave_lon,wave_lat,10,QWave==0,'filled')
%[~,b] = max(energy,[],3);
%imagesc(b), axis xy
%}
save GlobalDeltaData Hs Tp QWave -append 

end

function [FetchDeltas] = get_Qwave_fetch(MouthLon,MouthLat,lat,lon,FetchAll,QRiver_prist,delta_name)

lat = [lat{:}];
lon = [lon{:}];
lon2 = mod(lon-1,360)+1;

%find big landmasses first, avoids small islands.
isl = cumsum(isnan(lon2));
isl = isl./max(isl);

fetch = cell2mat(cellfun(@rot90,FetchAll,'UniformOutput',false));

FetchDeltas = zeros(length(MouthLon),360);

%make the coastline length for fetch calculation depend on delta size
rr = 3+max(0,round(QRiver_prist.^0.6));

for ii=1:length(MouthLon),
    [~,idx0] = min((MouthLat(ii)-lat).^2+(MouthLon(ii)-lon2).^2+(isl*(rr(ii)-3)./2));
    
    idx2 = find(isnan(lat(idx0+(1:rr(ii)))),1,'first');
    idx2(isempty(idx2)) = rr(ii);
    idx3 = find(isnan(lat(idx0+(-rr(ii):-1))),1,'last')-rr(ii);
    idx3(isempty(idx3)) = -rr(ii);
    
    idx1 = max(1,min(length(lat),idx0+(idx3:idx2)));
    
    FetchDeltas(ii,1:360) = max(fetch(:,idx1),[],2);


end
%figure,plot(FetchDeltas(ii,:))


%{
%plot
ii=165;

[~,idx0] = min((MouthLat(ii)-lat).^2+(MouthLon(ii)-lon2).^2+(isl*(rr(ii)-3)./2));

idx2 = find(isnan(lat(idx0+(1:rr(ii)))),1,'first');
idx2(isempty(idx2)) = rr(ii);
idx3 = find(isnan(lat(idx0+(-rr(ii):-1))),1,'last')-rr(ii);
idx3(isempty(idx3)) = -rr(ii);

idx1 = max(1,min(length(lat),idx0+(idx3:idx2)));

subplot(2,1,1)
cla
plot(lon2(idx1),lat(idx1)), hold on
plot(lon2(idx0),lat(idx0),'r'), hold on
scatter(MouthLon(ii),MouthLat(ii))
set(gca,'DataAspectRatio',[1 1 1])
subplot(2,1,2)
cla
%polarplot(deg2rad(-(0:359)),max(0,max(fetch(:,idxx+rr),[],2)),'r'), hold on,
polarplot(deg2rad(-(0:359)),max(0,FetchDeltas(ii,:)),'r'), hold on,
%figure,plot(max(fetch(:,idxx+rr),[],2)')
%}
%


%{

blub = num2cell([[1:length(MouthLon)]' FetchDeltasplot])';
for ii=1:length(blub),
    desc{ii} = sprintf('%1.0i: %1.0i - %1.0i - %1.0i - %1.0i',blub{:,ii});
end

kmlwritepoint(['blub.kml'],MouthLat(1:1000),MouthLon(1:1000),'name',desc(1:1000))

%save FetchDeltas FetchDeltas
%}



end