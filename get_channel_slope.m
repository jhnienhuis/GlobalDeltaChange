%% Get channel slope
load('D:\GlobalDatasets\GlobalDEM\CoastalZone_int8.mat');
for jj=1:7
    clearvars -except jj cz
    res=240;
    ele_max = 20; %maximum elevation to search for
    continents = {'na','af','ca','sa','eu','as','au'};
    remfun = @(lon) (rem(res*360-1+lon,res*360)+1);
    jj
    
%load accumulated drainage area (# cells)
d = ['D:\GlobalDatasets\HydroSheds\' continents{jj} '_acc_15s_bil\' continents{jj} '_acc_15s'];
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

%load located river mouths
load([dropbox filesep 'WorldDeltas' filesep 'scripts' filesep 'RiverMouth' continents{jj} '.mat'],'MouthLat','MouthLon','BasinArea')
BasinArea = single(BasinArea);
BasinLat = round((MouthLat+90)*res);
BasinLon = round(remfun(MouthLon*res));

%calculate drainage area (# cells) for the located river mouths
areapercell = 6371.^2.*2*pi/360/res*(sin(deg2rad((1/res:1/res:180)-90))-sin(deg2rad((1/res:1/res:180)-90-1/res)))';
BasinAreaCell = BasinArea'./areapercell(BasinLat);

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
    
    rlook = max(10,round(sqrt(BasinArea(ii))/10));
    
    %find more appropriate river mouth location using the accumulation cells
    [maxa,idxa,idxb] = max2d(-abs(imregmax(min(x,max(1,upstreamLat(ii)+(-rlook:rlook))),min(y,max(1,upstreamLon(ii)+(-rlook:rlook))))-int32(BasinAreaCell(ii))));
    
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

%transform back
channel_len = 1000*bsxfun(@times,sqrt(areapercell(BasinLat)),len); %m
channel_slope = bsxfun(@rdivide,[1 1:(ele_max-1)],channel_len);
channel_slope(channel_slope>100) = 1e-3;
channel_slope = mean(channel_slope(:,1:20),2);

len_lat(len_lat==0) = nan;
len_lon(len_lon==0) = nan;

channel_len_lat = (ulY-len_lat)./res-90;
channel_len_lon = (ulX+len_lon)./res;


maxa = max(channel_len,[],2);
ChannelSlope = zeros(size(maxa));
for ii=1:length(channel_len),
    if maxa(ii)<1,
        ChannelSlope(ii) = 1e-3;
    else,
        p=polyfit(channel_len(ii,:),log(1:20),1); 
        ChannelSlope(ii) = exp(p(2)).*p(1);
    end
end

BasinLat_acc = (ulY - BasinLat_acc)/res;
BasinLon_acc = (BasinLon_acc + ulX)/res;

save([dropbox filesep 'WorldDeltas'  filesep 'scripts' filesep 'RiverMouth' continents{jj} '.mat'],'channel_slope','ChannelSlope','channel_len_lat','channel_len_lon','channel_len','-append')

end