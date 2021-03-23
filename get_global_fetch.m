function get_global_fetch
%get wind, annual average
%get shoreline section per delta
%max per angle for entire shoreline section
%if open ocean then get from wavewatch
%if not, do fetch from wind speed


%load shore_vec
%out = load('shore_vec_h.mat');
%noaa shoreline vector, available at https://www.ngdc.noaa.gov/mgg/shorelines/
s = shaperead('GSHHS_h_L1.shp');

%axis equal
%hold on
%35 & 47 are nova zembla
%long island: 157, north america: 4
%for ii=1:10,
%    plot(lon{ii},lat{ii})
%    text(double(lon{ii}(1)),double(lat{ii}(1)),num2str(ii))
%end
lon = {s(:).X};
lat = {s(:).Y};
lon_l = double([s(:).X]);
lat_l = double([s(:).Y]);
%split up in 1 degree boxes, only look at coasts within 1 degree
deg1 = floor(2.*lon_l+2i*lat_l);

%but loop through island?
islands = length(lat);
FetchAll = cell(size(islands));

for kk=9:islands
lon_isl = double([lon{kk}]);
lat_isl = double([lat{kk}]);
deg_isl = floor(2.*lon_isl+2i*lat_isl);
kk
%fetch map
n = length(lon_isl);
nspace = 1;
fetch = zeros(n,360,'int8');

for ii=1:nspace:n,
    %if mod(ii,1000)==1, progressbar(ii/n), end
    
    if isnan(lon_isl(ii)), continue, end
    
    %only look at coasts within +/- 0.5 degree
    sbox = deg_isl(ii)+(-1:1)+1i*(-1:1)';
    
    %select all coasts within 1 degree
    if ii==1 || deg_isl(ii)~=deg_isl(ii-nspace),
        idx = ismember(deg1,sbox);
        bridge = find(diff(find(idx))>1);
    end
    
    %convert to polar coordinates
    [theta,rho] = cart2pol(lon_l(idx)-lon_isl(ii),lat_l(idx)-lat_isl(ii));
    
    %convert angles to degrees
    theta = min(360,floor(mod((180*theta/pi)-1,360)+1));
    
    %convert distances to km (roughly)
    rho = 6371*rho*pi/180;
    
    %find selected coastline again, to block land
    idxx = find(rho==0,1);
    rho(idxx) = 999;
    
    %this is easy, just find the nearest coast per angle
    %unlimited fetch is 127km
    fetch(ii,:) = accumarray(theta',rho,[],@min);
    fetch(ii,fetch(ii,:)==0) = 127;
        
    %what are its neighbours
    idx1 = theta(max(1,idxx-1));
    idx2 = theta(min(length(theta),idxx+1));
    
    %use neighbours to block land
    if idx1>idx2,
        fetch(ii,[idx1:360, 1:idx2]) = -128;
    else,
        fetch(ii,idx1:idx2) = -128;
    end
    
    %now more difficult part, what if part of the coast is not covered by
    %the coastal node.
    %unwrap angles
    blob = round(180*unwrap(pi*theta/180)/pi);
    %figure, plot(blob,'-o'), hold on, scatter(idx1, blob(idx1))
    %get rid of unconnected parts of the shoreline (islands etc)
    blob([idxx bridge]) = nan;
    
    %find differences greater than 1 degree
    blobd = diff(blob);
    blobdf= find(abs(blobd)>1.1);
    blobds = sign(blobd(blobdf));
    
    %interpolate between these angles, give it the average distance
    for jj=1:size(blobdf,2),
        
        mm = floor(mod((blob(blobdf(jj)):blobds(jj):blob(blobdf(jj)+1))-1,360)+1);
        
        fetch(ii,mm) = min(fetch(ii,mm),(rho(blobdf(jj))+rho(blobdf(jj)+1))/2);
    end
    
end

%if coast is closer than 1 km, call it land.
fetch(fetch==0) = -128;

%put back in thing per island..
FetchAll{kk} = fetch;
end
save GlobalShorelineFetch FetchAll lon lat

%{
 save to shapefile
load GlobalFetch

Latitude = [lat_l{:}];
Longitude = [lon_l{:}];
Fetch = vertcat(FetchAll{:});
Fetch(Fetch<0) = 0;

name = [{'x';'y'}; cellstr([repmat('f',360,1)])];% num2str((1:360)','%-i')])];

csvwrite('GlobalFetch3.csv',name')
csvwrite('GlobalFetch3.csv',[Longitude(1)',Latitude(1)',Fetch(1,:)],1,0);

% shapefile is annoying
load GlobalFetch
s = struct;
s.Latitude = [lat_l{:}];
s.Longitude = [lon_l{:}];
Fetch = double(vertcat(FetchAll{:}));
Fetch(Fetch<0) = 0;

for ii=1:360,
    s.(['f' num2str(ii)]) = Fetch(:,ii);
end
s = geopoint(s);
s.Geometry = 'point';
dbf = makedbfspec(s);

for ii=1:360,
    dbf.(['f' num2str(ii)]).FieldDecimalCount = 0;
    dbf.(['f' num2str(ii)]).FieldLength = 1;
end

shapewrite(s,'GlobalFetch.shp','DbfSpec',dbf)
%}
%%plotting
%{
ii=4001; jj=10;
figure
plot([lon{:}],[lat{:}],'-')
hold on
axis equal

scatter([lon{jj}],[lat{jj}],20,max([FetchAll{jj}],[],2),'filled');
set(gca,'CLim',[0 127]);
scatter(lon{jj}(ii),lat{jj}(ii),100,'g','filled');

figure
polarplot(max(0,FetchAll{jj}(ii,:)))
%}