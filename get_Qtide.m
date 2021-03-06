%% Get tidal fluxes
load([dropbox filesep 'WorldDeltas' filesep 'scripts' filesep 'GlobalDeltaData.mat'])

direc = 'D:\OneDrive - Universiteit Utrecht\Tides\';
lat_tide = ncread([direc 'hf.m2_tpxo8_atlas_30c_v1.nc'],'lat_z');
lon_tide = ncread([direc 'hf.m2_tpxo8_atlas_30c_v1.nc'],'lon_z');

m2 = int32(abs(single(ncread([direc 'hf.m2_tpxo8_atlas_30c_v1.nc'],'hRe'))+1i*single(ncread([direc 'hf.m2_tpxo8_atlas_30c_v1.nc'],'hIm'))));
m4 = int32(abs(single(ncread([direc 'hf.m4_tpxo8_atlas_30c_v1.nc'],'hRe'))+1i*single(ncread([direc 'hf.m4_tpxo8_atlas_30c_v1.nc'],'hIm'))));
s2 = int32(abs(single(ncread([direc 'hf.s2_tpxo8_atlas_30c_v1.nc'],'hRe'))+1i*single(ncread([direc 'hf.s2_tpxo8_atlas_30c_v1.nc'],'hIm'))));
k1 = int32(abs(single(ncread([direc 'hf.k1_tpxo8_atlas_30c_v1.nc'],'hRe'))+1i*single(ncread([direc 'hf.k1_tpxo8_atlas_30c_v1.nc'],'hIm'))));
k2 = int32(abs(single(ncread([direc 'hf.k2_tpxo8_atlas_30c_v1.nc'],'hRe'))+1i*single(ncread([direc 'hf.k2_tpxo8_atlas_30c_v1.nc'],'hIm'))));
o1 = int32(abs(single(ncread([direc 'hf.o1_tpxo8_atlas_30c_v1.nc'],'hRe'))+1i*single(ncread([direc 'hf.o1_tpxo8_atlas_30c_v1.nc'],'hIm'))));
p1 = int32(abs(single(ncread([direc 'hf.p1_tpxo8_atlas_30c_v1.nc'],'hRe'))+1i*single(ncread([direc 'hf.p1_tpxo8_atlas_30c_v1.nc'],'hIm'))));
q1 = int32(abs(single(ncread([direc 'hf.q1_tpxo8_atlas_30c_v1.nc'],'hRe'))+1i*single(ncread([direc 'hf.q1_tpxo8_atlas_30c_v1.nc'],'hIm'))));
n2 = int32(abs(single(ncread([direc 'hf.n2_tpxo8_atlas_30c_v1.nc'],'hRe'))+1i*single(ncread([direc 'hf.n2_tpxo8_atlas_30c_v1.nc'],'hIm'))));



%hamax = reshape(max([m2(:),s2(:),k1(:),o1(:)],[],2),size(s2));
hamax = (m2+m4+s2+k1+k2+o1+p1+q1+n2)/2;

%clearvars m2 m4 s2 k1 k2 o1 p1 q1 n2

F = ((k1 + o1)./(s2 + m2))>1.5;

hamax(hamax<1) = 0; %in mm!
res = 30;

remfun = @(lon) (rem(res*360-1+lon,res*360)+1);

mouth_lat = floor(1+res*(90+MouthLat));
mouth_lon = floor(remfun(res*MouthLon));

tide_a = zeros(size(MouthLon));

%plot for error check
%imagesc(hamax,[0 4000]), hold on, scatter(mouth_lon,mouth_lat,'or','filled'), axis xy
%
for ii=1:numel(MouthLat),
    tide_a(ii) = hamax(mouth_lat(ii),mouth_lon(ii));
    
    if tide_a(ii)==0,
        
            sbox = 10; mid_idx = 11;
            [loc_y,loc_x] = find(hamax(mouth_lat(ii)+(-sbox:sbox),remfun(mouth_lon(ii)+(-sbox:sbox)))>0);
            
        
        if isempty(loc_y),
            sbox = 100; mid_idx = 101;
            [loc_y,loc_x] = find(hamax(mouth_lat(ii)+(-sbox:sbox),remfun(mouth_lon(ii)+(-sbox:sbox)))>0);
            
            if isempty(loc_y),
                tide_a(ii) = nan;
            end
            
        end
    
        if ~isempty(loc_y),
            
            [~,tideloc] = min(abs((loc_x-mid_idx).^2+(loc_y-mid_idx).^2));
             
            mouth_lat(ii) = mouth_lat(ii)-mid_idx+loc_y(tideloc);
            mouth_lon(ii) = remfun(mouth_lon(ii)-mid_idx+loc_x(tideloc));
            
            tide_a(ii) = hamax(mouth_lat(ii),remfun(mouth_lon(ii)));
            
        end
    else,
        tide_a(ii) = hamax(mouth_lat(ii),remfun(mouth_lon(ii)));
    end

end

tide_omega = single(F(sub2ind(size(F),mouth_lat,mouth_lon)));
tide_omega(tide_omega==0) = 1.4e-4;
tide_omega(tide_omega==1) = 0.7e-4;
tide_a = tide_a./1000;
Discharge_prist(Discharge_prist<=0) = 1e-10;
QRiver_prist(QRiver_prist<=0) = 1e-10;
ChannelSlope(ChannelSlope==0) = 1e-4;

w_upstream = 10*Discharge_prist.^0.5;
d_upstream = 0.3*Discharge_prist.^0.35;
beta = w_upstream./d_upstream;
t_length = max(2*tide_a,d_upstream)./ChannelSlope;
t_length = min(t_length,pi./tide_omega.*sqrt(d_upstream.*10)); %put limiter on tidal intrusion distance based on tidal wave celerity

k = tide_omega./(pi*sqrt(0.2*1e-4)*1.65*55);

Discharge_tide = double(0.5*tide_omega.*k.*tide_a.^2.*t_length.^2.*beta.*(1+2.*ChannelSlope./(k.*tide_a))); %Qtide in m3/s of water %corrected!
QTide = Discharge_tide.*min(1,QRiver_prist./Discharge_prist); %Qtide in kg/s of sediment assuming same sediment concentration
TidalAmp = tide_a;

MouthLon(MouthLon<0) = MouthLon(MouthLon<0) + 360;

save('GlobalDeltaData.mat','QTide','Discharge_tide','ChannelSlope','QRiver_prist','TidalAmp','Discharge_prist','-append');