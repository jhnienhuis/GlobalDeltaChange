clr
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaProfile','bed_h')
load('D:\Drive\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaArea.mat','delta_area','src')
load('D:\Drive\github\GlobalDeltaChange\GlobalDeltaData.mat','BasinID2','BasinArea','MouthLon','MouthLat','delta_name');


%remove caspian sea and other non deltas and deltas > 10km2
idx = (~isnan(bed_h) & bed_h<0) & ~(MouthLon>46 & MouthLon<56 & MouthLat>34 & MouthLat<48) & delta_area>1e7 &src==1;

MouthLat = MouthLat(idx)';
MouthLon = MouthLon(idx)';

MouthLon(MouthLon>180) = MouthLon(MouthLon>180)-360;

a = geotiffread('D:\OneDrive - Universiteit Utrecht\WBMSed\bqart_sed-flow_dir\bqart_a.tif');
a(a<0) = nan; a = flipud(a);
d = matfile(['D:\OneDrive - Universiteit Utrecht\WBMSed\qs_timeseries']);
wbm_lat = d.rm_lat;
wbm_lon = d.rm_lon;
lon_grid = d.lon_grid;
lat_grid = d.lat_grid;
wbm_basinarea = a(sub2ind(size(a),wbm_lat,wbm_lon));

delta_coor = MouthLat+1i.*MouthLon;
wbm_coor = lat_grid(wbm_lat)+1i.*lon_grid(wbm_lon);

blub = abs(wbm_coor-delta_coor) +...
    (abs((wbm_basinarea-BasinArea(idx)')./BasinArea(idx)'.*log10(BasinArea(idx)')))+...
    (abs(wbm_coor-delta_coor)>8).*10;

[~,idx2] = min(blub);
fac = wbm_basinarea(idx2)./BasinArea(idx);
%mackenzie = [459,1595], colville = [296,1605], lena=[3035,1632]

discharge_series = d.discharge_series;
discharge_series = discharge_series(:,idx2).*fac';


Discharge50 = prctile(discharge_series,50,1)';
Discharge90 = prctile(discharge_series,90,1)';
Discharge99 = prctile(discharge_series,99,1)';
Discharge999 = prctile(discharge_series,99.9,1)';




BasinID2 = BasinID2(idx);
delta_name = delta_name(idx);


t = table(BasinID2,delta_name,Discharge50,Discharge90,Discharge99,Discharge999);

writetable(t,'GlobalDeltaDischargeExtremes.xlsx')



%[yr_t,~] = datevec(d.t);
%yr_t = yr_t-1979;
%rm_lat = (d.rm_lat-900)./10;
%mean annual maximum
% pk = zeros(size(rm_lat));
% m = zeros(size(rm_lat));
% p10 = zeros(size(rm_lat));
% p90 = zeros(size(rm_lat));
% 
% for idx=1:length(rm_lat)
%     if rem(idx,100)==1, disp(idx),end
%     x = d.discharge_series(:,idx);
%     m(idx) = mean(x);
%     pk(idx) = mean(accumarray(yr_t,x,[],@max));
%     p10(idx) = mean(accumarray(yr_t,x,[],@(x) (prctile(x,10))));
%     p90(idx) = mean(accumarray(yr_t,x,[],@(x) (prctile(x,90))));
% end






