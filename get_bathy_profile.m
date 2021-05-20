function get_bathy_profile
%load delta database
load('GlobalDeltaData.mat','MouthLon','MouthLat');

%get continental shelf depth

%load world bathymetry
%[cz, refvec] = etopo('D:\OneDrive - Universiteit Utrecht\GlobalDEM\',5,[-90 90]);
%cz = cz(:,[(180*refvec(1)+1):end 1:(180*refvec(1))]);

load('D:\OneDrive - Universiteit Utrecht\GlobalDEM\SRTM15plus_int8.mat','cz');
cz(cz>0) = 0;
cz = double(cz);
refvec = [240,90,-180];

rsample = 20; cz = cz(1:rsample:end,1:rsample:end); refvec = [60/(0.25*rsample) 90 -180];



%hoy many km per cell latitude (simple!)
kmperlatcell = deg2km(1)./refvec(1);

%how many km per cell longitude (depends on latitude)
kmperloncell = kmperlatcell*cos(linspace(-0.5*pi,0.5*pi,size(cz,1)));

%find the coastline contour (read the contourc to see what the output looks like
%C = contourc(Z,[0 0]);

%uncomment to plot 0 and 200m contour line
%find all instances of the coastline (all "islands")
%I = find(C(1,:)==200);

%hold on
%for ii=1:length(I)-1,
%    plot(C(1,(1+I(ii)):(I(ii+1)-1)),C(2,(1+I(ii)):(I(ii+1)-1)))
%end



%turn delta locations into imag
delta_loc = (refvec(1)*MouthLon)+1i.*(refvec(1).*(MouthLat+90));

%plot(delta_loc,'or')
%plot 0 m contour
%C0 = contourc(Z,[0 0]);
%I = find(C0(1,:)==0);
%hold on
%for ii=1:length(I)-1,
%    plot(C0(1,(1+I(ii)):(I(ii+1)-1)),C0(2,(1+I(ii)):(I(ii+1)-1)),'b')
%end


%do profiles for these contour lines
shelf_lines = 0:-5:-100;

%pre-allocate matrices
shelf_len = zeros(length(delta_loc),length(shelf_lines));
shelf_len_lat = zeros(length(delta_loc),length(shelf_lines));
shelf_len_lon = zeros(length(delta_loc),length(shelf_lines));

shelf_len_lat(:,1) = imag(delta_loc);
shelf_len_lon(:,1) = real(delta_loc);



%find minimum distance to the next contour line for every contour
for jj=1:length(shelf_lines),
    jj
    S = contourc(cz,[shelf_lines(jj) shelf_lines(jj)]);
    S2 = S(1,S(1,:)~=shelf_lines(jj))+1i*S(2,S(1,:)~=shelf_lines(jj));
    
    %idx = zeros(length(delta_loc),1);
    %for all the shoreline positions
    for ii=1:length(delta_loc),
        
        %distance to contour in km (simple minimum of the complex magnitude
        %of the difference the two contour lines
        [shelf_len(ii,jj),idx] = min(abs((kmperloncell(uint16(imag(delta_loc(ii))))+1i*kmperlatcell).*(S2-delta_loc(ii))));
        
        %use new contour locations as the next place to start finding minima.
        delta_loc(ii) = S2(idx);
        shelf_len_lat(ii,jj) = imag(S2(idx));
        shelf_len_lon(ii,jj) = real(S2(idx));
        
    end
    

    
end
%snap deltas to coastline.
shelf_len(:,1) = 0;

shelf_len_lat(shelf_len_lat==0) = nan;
shelf_len_lon(shelf_len_lon==0) = nan;

%turn matrix indices into latitude and longitude
shelf_len_lat = (shelf_len_lat./refvec(1))-90;
shelf_len_lon = shelf_len_lon./refvec(1);

shelf_len = cumsum(shelf_len,2).*1000; %turn to m

%remove distances away from shoreline larger than 400km? (great lakes etc)
shelf_len_lat(shelf_len>5e5) = nan;
shelf_len_lon(shelf_len>5e5) = nan;
shelf_len(shelf_len>5e5) = nan;

% plot a particular transect
%figure
%ii=500;
%plot(shelf_len(ii,:),c_lines,'-o')
%xlabel('Offshore Distance (m)')
%ylabel('Depth (m)')
%hold on


% add shelf to global delta dataset
save GlobalDeltaData.mat shelf_len shelf_len_lon shelf_len_lat shelf_lines -append