function [sed_series,discharge_series,t] = get_Qriver_timeseries(rm_lat,rm_lon,basinarea)

rm_lon(rm_lon>180) = rm_lon(rm_lon>180)-360;
%
%colville
%lena
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

%qs_day = double(accumarray(t_day,sed_series,[366 1],@mean));
%discharge_day = double(accumarray(t_day,discharge_series,[366 1],@mean));
%num_day_50pct = find(cumsum(sort(qs_day,'descend'))>(0.5*sum(qs_day)),1);
%if isempty(num_day_50pct), num_day_50pct = nan; end

end