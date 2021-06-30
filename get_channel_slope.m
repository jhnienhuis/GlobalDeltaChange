%% Get channel slope
function get_channel_slope

%load located river mouths
load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','MouthLat','MouthLon','BasinArea','Continent','RiverID_ATLAS','delta_name')
%load('D:\OneDrive - Universiteit Utrecht\GlobalDEM\CoastalZone_int8.mat','cz');
load('D:\OneDrive - Universiteit Utrecht\GlobalDEM\SRTM15plus_int8.mat','cz');
sz = size(cz);
cz(cz<0) = 0;
res=240;
ele_max = 30; %maximum elevation to search for
Continent_name = {'na','af','ca','sa','eu','as','au'};
remfun = @(lon) (rem(res*360-1+lon,res*360)+1);

channel_len = zeros(length(MouthLat),ele_max);
channel_len_lat = zeros(length(MouthLat),ele_max);
channel_len_lon = zeros(length(MouthLat),ele_max);


for jj=1:7
    %clearvars -except jj cz MouthLon MouthLat BasinArea Continent res ele_max continents remfun channel_len channel_len_lat channel_len_lon
    
    jj
    
    
    idxc = find(Continent==jj);
    
    BasinAreac = single(BasinArea(idxc));
    BasinLat = round((MouthLat(idxc)+90)*res);
    BasinLon = round(remfun(MouthLon(idxc)*res));
    
    [len, len_lat, len_lon] = get_channel_slope_f(cz,BasinLon,BasinLat,BasinAreac,Continent_name{jj},res,ele_max,remfun);
    
    channel_len(idxc,:) = len;
    channel_len_lat(idxc,:) = len_lat;
    channel_len_lon(idxc,:) = len_lon;
    
    
end

channel_len_lat(channel_len_lat==0) = nan;
channel_len_lon(channel_len_lon==0) = nan;


save('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','channel_len_lat','channel_len_lon','channel_len','-append')
load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','channel_len_lat','channel_len_lon','channel_len')

%do other function for channels that were not found:

Continent_name = {'ar','eu','na','si','au','as','af','sa_north','sa_south'};
res=240;


for jj=1:9
    
    jj
    
    %load shapefile
    x = shaperead(['D:\OneDrive - Universiteit Utrecht\HydroSheds\RiverATLAS_Data_v10_shp\RiverATLAS_v10_shp\RiverATLAS_v10_' Continent_name{jj} '_filt2.shp']);
    
    
    shape_ID = [x.MAIN_RIV];
    
    %find all rivers on continent and %find all rivers where previous method didn't work
    [idx] = find(ismember(RiverID_ATLAS',shape_ID)& (channel_len(:,ele_max)==0) ); %
    
    
    %loop through rivers/deltas
    for ii=1:length(idx),
        
        %find all river sections that belong to delta
        idxr = find(shape_ID==RiverID_ATLAS(idx(ii)));
        
        %sort from mouth to upstream
        [~,idxrs] = sort([x(idxr).DIST_DN_KM],'descend');
        
        %extract lat/lon, put from mouth to upstream
        RIVlat = fliplr([x(idxr(idxrs)).Y]);
        RIVlon = fliplr([x(idxr(idxrs)).X]);
        
        idnan = isnan(RIVlat);
        
        RIVlat(idnan) = [];
        RIVlon(idnan) = [];
        
        %do distance between points
        di = [(deg2km(distance(RIVlat(1:end-1),RIVlon(1:end-1),RIVlat(2:end),RIVlon(2:end))))];
        
        %idxint = di==0;
        %RIVlat(idxint) = [];
        %RIVlon(idxint) = [];
        %di(idxint) = [];
        
        %do cumulative length of river
        di = cumsum([0 di]);
        
        %convert channel coordinates to grid indices
        xx = round(remfun(res*RIVlon+2));
        yy = round((res.*(90+RIVlat))+1);
        
        
        %plot(xx(~idnan),yy(~idnan))
        
        %put into linear index for path extraction
        id = sub2ind(sz,yy,xx);
        
        %retrieve elevation from DEM, rivers can only flow downstream,
        %hence the cummax or the reverse cumulative cummin
        %el = cummax((double(cz(id))));
        el = cummax((double(cz(id)))); %,'reverse');
        el(1) = 0; el(2) = min(el(2),ele_max); el(end) = max(el(end),1);
        %el(~idnan) = el;
        %el(idnan) = nan;
        %imagesc(cz(min(yy):max(yy),min(xx):max(xx))), hold on, plot(xx-(min(xx)),yy-min(yy),'or')
        %scatter(RIVlon,RIVlat,30,el,'filled')
        
        %find unique elevations up to max (default max is 50 m above sea level)
        [uniq,uniqi] = unique(el(el<(ele_max)));
        if length(uniq)==1, uniq(2) = 1; uniqi(2) = 2; end
        %get channel distance from mouth for these unique elevations
        channel_len(idx(ii),uniq+1) = cummax(di(uniqi).*1000);
        
        %get coordinates from profile
        channel_len_lat(idx(ii),uniq+1) = (RIVlat(uniqi));
        channel_len_lon(idx(ii),uniq+1) = (RIVlon(uniqi));
        
    end
    %
end
%}
%{
%four odd rivers not found?!?!
idx_notfound = find(max(abs(channel_len_lat),[],2)==0);
channel_len_lat(idx_notfound,:) = repmat(MouthLat(idx_notfound),[1 ele_max+1]);
channel_len_lon(idx_notfound,:) = repmat(MouthLon(idx_notfound),[1 ele_max+1]);
%}



maxa = max(channel_len,[],2);
channel_slope = zeros(size(maxa));
for ii=1:length(channel_len),
    if maxa(ii)<1,
        channel_slope(ii) = 1e-3;
    else,
        p=polyfit(channel_len(ii,:),log(1:ele_max),1);
        channel_slope(ii) = exp(p(2)).*p(1);
    end
end

channel_slope(isnan(channel_slope) | channel_slope<2e-5) = 2e-5;

%do mod to squeeze longitudes on a 0-360 grid.
for ii=1:length(channel_len),
    x = find(~isnan(channel_len_lon(ii,:)),1);
    if isempty(x),
        channel_len_lon(ii,:) = MouthLon(ii);
        channel_len_lat(ii,:) = MouthLat(ii);
    else,
    x = (mod(channel_len_lon(ii,x)-1,360)+1)-channel_len_lon(ii,x);
    channel_len_lon(ii,:) = channel_len_lon(ii,:)+x;
    end
    
end



save('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','channel_slope','channel_len_lat','channel_len_lon','channel_len','-append')


end


function [len, len_lat, len_lon] = get_channel_slope_f(cz,BasinLon,BasinLat,BasinAreac,Continent_name,res,ele_max,remfun)

%load accumulated drainage area (# cells)
d = ['D:\OneDrive - Universiteit Utrecht\HydroSheds\' Continent_name '_acc_15s_bil\' Continent_name '_acc_15s'];
fileID = fopen([d '.hdr'],'r');
hdr = textscan(fileID, '%s%s%[^\n\r]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true); fclose(fileID);
hdr = str2double(hdr{2});
ulX = round(res*rem(hdr(11)+360,360));
ulY = round(res*(90+hdr(12)));
a = int32(multibandread([d '.bil'],[hdr(3) hdr(4) hdr(5)],'int32',0,'bil','ieee-le'));
a(a<10) = 0;
disp('Loaded drainage accumulation')
imregmax = a.*int32(imregionalmax(a)); %find river mouths again
disp('Calculated regional maxima')
x = size(a,1); y = size(a,2);

%trim file
cz2 = cz(ulY-(0:(x-1)),remfun(ulX+(0:(y-1))));



%calculate drainage area (# cells) for the located river mouths
areapercell = 6371.^2.*2*pi/360/res*(sin(deg2rad((1/res:1/res:180)-90))-sin(deg2rad((1/res:1/res:180)-90-1/res)))';
BasinAreaCell = BasinAreac./areapercell(BasinLat);

%change grid to match .bil file
BasinLatSm = min(x,max(1,ulY-BasinLat));
BasinLonSm = min(y,max(1,remfun(BasinLon-ulX)));

BasinLat_acc = zeros(size(BasinLat));
BasinLon_acc = zeros(size(BasinLon));
upstreamLon = BasinLonSm;
upstreamLat = BasinLatSm;

%search parameters
dis = [sqrt(2) 1 sqrt(2); 1 0 1; sqrt(2) 1 sqrt(2)];
len = zeros(length(upstreamLat),ele_max);
len_lat = zeros(length(upstreamLat),ele_max);
len_lon = zeros(length(upstreamLat),ele_max);

for ii=1:length(upstreamLat),
    
    if mod(ii,200)==1, ii, end
    
    rlook = max(10,round(sqrt(BasinAreac(ii))/10));
    
    %find more appropriate river mouth location using the accumulation cells
    [~,idxa,idxb] = max2d(-abs(imregmax(min(x,max(1,upstreamLat(ii)+(-rlook:rlook))),min(y,max(1,upstreamLon(ii)+(-rlook:rlook))))-int32(BasinAreaCell(ii))));
    
    %somehow a very big river was not found
    
    upstreamLat(ii) = upstreamLat(ii)-rlook-1+idxa;
    upstreamLon(ii) = upstreamLon(ii)-rlook-1+idxb;
    
    %save that variable
    BasinLat_acc(ii) = upstreamLat(ii);
    BasinLon_acc(ii) = upstreamLon(ii);
    
    %if nothing found, return len==-1; this happens most often if the
    %channel has been found already.
    if a(upstreamLat(ii),upstreamLon(ii))<(0.5*BasinAreaCell(ii)),
        len(ii,:) = nan;
        len_lat(ii,:) = nan;
        len_lon(ii,:) = nan;
        continue,
    end
    
    %else, start search up to elevation = 15m
    k=0;
    while cz2(upstreamLat(ii),upstreamLon(ii))<ele_max,
        
        k=k+1;
        %set current cell to -999 to prevent it to be found in max2d
        a(upstreamLat(ii),upstreamLon(ii)) = -999;
        
        %find next maximum
        [maxa,idxa,idxb] = max2d(a(min(x,max(1,upstreamLat(ii)+(-1:1))),min(y,max(1,upstreamLon(ii)+(-1:1)))));
        
        %if less than 10 drainage cells, stop
        if maxa<10 || k==1000,
            break,
        end
        
        %set cell to new maximum
        upstreamLat(ii) = upstreamLat(ii)-2+idxa;
        upstreamLon(ii) = upstreamLon(ii)-2+idxb;
        
        %calculate elevation of cell
        ele = cz2(upstreamLat(ii),upstreamLon(ii))+1;
        len(ii,ele:end) = len(ii,ele:end)+dis(idxa,idxb);
        len_lat(ii,ele:end) = upstreamLat(ii);
        len_lon(ii,ele:end) = upstreamLon(ii);
    end
    
end

len = 1000*bsxfun(@times,sqrt(areapercell(BasinLat)),len); %m

len_lat(len_lat==0) = nan;
len_lon(len_lon==0) = nan;

len_lat = (ulY-len_lat)./res-90;
len_lon = (ulX+len_lon)./res;
%make image of this
%{
imagesc(cz2)
blub = bwboundaries(a==-999);
hold on
for ii=1:length(blub),
    plot(blub{ii}(:,2),blub{ii}(:,1),'w')
end
axis ij, axis tight
scatter(upstreamLon,upstreamLat,'r','filled','MarkerEdgeColor','w')
scatter(BasinLon_acc, BasinLat_acc,'g','filled','MarkerEdgeColor','w')
set(gca,'DataAspectRatio',[1 1 1],'XLim',([264 267]*240)-ulX,'YLim',ulY-([109 107.5]*240))
%}
end