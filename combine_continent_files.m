%% combine continent files

%merge all files into one mat file
continents = {'na','af','ca','sa','eu','as','au','ao'};
remfun = @(lon) (rem(360-1+lon,360)+1);

k=0;
for jj = 1:8,
    load([dropbox filesep 'WorldDeltas\scripts\Rivermouth' continents{jj} '.mat'],'BasinArea');
    sz = length(BasinArea);
    Continent(k+(1:sz)) = jj;
    k=k+sz;
end
Continent = Continent';
BasinArea = zeros(k,1);
MouthLat = zeros(k,1);
MouthLon = zeros(k,1);
ChannelSlope = zeros(k,1);
QRiver_dist = zeros(k,1);
QRiver_prist = zeros(k,1);
Discharge_dist = zeros(k,1);
Discharge_prist = zeros(k,1);
BasinID = zeros(k,1);

for jj = 1:8,
    out = load([dropbox filesep 'WorldDeltas\scripts\Rivermouth' continents{jj} '.mat']);
    
    if jj==8,
        fact = 1;
        MouthLon(Continent==jj) = out.MouthLon;
        MouthLat(Continent==jj) = out.MouthLat;
        BasinID(Continent==jj) = 1:length(out.MouthLat);
    else,
        fact = double(out.BasinArea')./out.cellareabasins;
        MouthLat(Continent==jj) = out.MouthLat;
        MouthLon(Continent==jj) = out.MouthLon;
        ChannelSlope(Continent==jj) = out.channel_slope;
        
       
        BasinID(Continent==jj) = out.BasinID;
        nnz(unique(out.BasinID))
        nnz(out.BasinID)
    end
        
    BasinArea(Continent==jj) = out.BasinArea;
    QRiver_dist(Continent==jj) = out.QRiver_dist.*fact;
    QRiver_prist(Continent==jj) = out.QRiver_prist.*fact;
    Discharge_dist(Continent==jj) = out.Discharge_dist.*fact;
    Discharge_prist(Continent==jj) = out.Discharge_prist.*fact;
    
end

idx = Discharge_prist>1;
idx = idx & (QRiver_prist>0.05 | QRiver_dist>0.05 | Discharge_prist>20);
QRiver_prist = max(0,QRiver_prist);
QRiver_dist = max(0,QRiver_dist);
QRiver_prist = QRiver_prist(idx);
QRiver_dist = QRiver_dist(idx);
BasinID = BasinID(idx);
BasinArea = BasinArea(idx);
MouthLat = MouthLat(idx);
MouthLon = remfun(MouthLon(idx));
Discharge_dist = Discharge_dist(idx);
Discharge_prist = Discharge_prist(idx);
ChannelSlope = ChannelSlope(idx);
Continent = Continent(idx);

% add continent flag
diva = shaperead('D:\OneDrive - Universiteit Utrecht\DIVA\cls_p18_2.shp','UseGeoCoords',true);
remfun = @(lon) (rem(360-1+lon,360)+1);
diva_lon = remfun([diva(:).Lon]);
diva_lat = ([diva(:).Lat]);
diva_size = arrayfun(@(x) length(x.Lon),diva);
diva_type = cell2mat(arrayfun(@(x,y) (x.*ones(1,y)),[diva(:).CPC],diva_size','UniformOutput',false));
diva_cont = cell2mat(arrayfun(@(x,y) (x.*ones(1,y)),[diva(:).GVARID],diva_size','UniformOutput',false));
[~,ia] = unique(int32(diva_lat) + 1000*int32(diva_lon));
diva_lon = diva_lon(ia); diva_lat = diva_lat(ia); diva_type = diva_type(ia); diva_cont = diva_cont(ia);
%scatter(diva_lon,diva_lat,40,diva_cont)

diva_both = diva_lon + 1i*diva_lat;
delta_both = MouthLon + 1i*MouthLat;
[~,x] = min(abs(repmat(delta_both,1,length(diva_both))-diva_both),[],2);
Region = int32(diva_cont(x))';
%scatter(MouthLon,MouthLat,30,cont)
hold on,
l = unique(diva_cont);
for ii=1:length(l),
    
    scatter(MouthLon(Region==l(ii)),MouthLat(Region==l(ii)),30,Region(Region==l(ii)),'filled')
end

Region(Region==0) = 3;
Region(Region==23) = 18;
Region(Region>14) = Region(Region>14)-1;
Region_str = {'East Africa', 'South Asia','West Africa','Baltic','Eastern Central America','Western Central America','Carribean','Russia','East Asia','Middle-East','Northern Africa','Eastern North America','Western North America','Northern Mediterranean','Western Europe','Oceania','Pacific Islands','Eastern South America','Western South America','Southeast Asia'};


% save files
BasinID2 = int64(BasinID.*10+Continent);


save([dropbox filesep 'WorldDeltas\scripts\GlobalDeltaData.mat'],'BasinArea','MouthLat','MouthLon','ChannelSlope','channel_len','channel_len_lat','channel_len_lon',...
    'QRiver_dist','QRiver_prist','Discharge_dist','Discharge_prist','BasinID','BasinID2','Continent','Region','Region_str','-append');

%kmlwrite('GlobalDeltaData',MouthLat,MouthLon,'Name',cellstr(num2str(BasinArea,'%1.0e')))
%kmlwrite('GlobalDeltaBasinID',MouthLat,MouthLon,'Name',cellstr(num2str(BasinID,'%1.0f')),'Description',cellstr(num2str(BasinArea,'%1.0e')))