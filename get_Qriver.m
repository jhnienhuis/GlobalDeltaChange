function get_Qriver
% get flux per river mouth, for all continents except arctic
load([dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep 'GlobalDeltaData.mat'],'Continent','MouthLat','MouthLon','BasinArea','BasinID','QRiver_prist','Discharge_prist','QRiver_dist','Discharge_dist');

c_list = {'na','af','ca','sa','eu','as','au'};

%resolution is 240 per degree (15s) (test)
res = 240;
%resolution WBMSed is 10 per degree (6 min)
res_wbm = 10;
remfun = @(lon) (rem(res_wbm*360-1+lon,res_wbm*360)+1);

PristYield = load('D:\GlobalDatasets\WBMSed\Prist_Yield.mat');
DistYield = load('D:\GlobalDatasets\WBMSed\Dist_Yield.mat');

%area per cell of WBMSed in km2
areapercell = repmat(6371.^2.*2*pi/360/res_wbm*(sin(deg2rad(-89.9:0.1:90))-sin(deg2rad(-90:0.1:89.9))),360*res_wbm,1);

for jj=1:length(c_list),
    
    %load river mouth info
    %load([dropbox filesep 'WorldDeltas\scripts\Rivermouth' continents{jj} '.mat'])
    jj
    %get all north american river basins > 50km2
    basins = shaperead(['D:\GlobalDatasets\HydroSheds\' c_list{jj} '_bas_15s_beta\' c_list{jj} '_bas_15s_beta.shp'],'Selector',{@(v1) (v1>40),'AREA_SQKM'});
    disp('Loaded basins...')
    basins_list = [basins(:).BASIN_ID];
    
    h=waitbar(0);
    for ii=find(Continent==jj),
        
        
        
        waitbar(ii/length(BasinID))
        
        kk=find(basins_list==BasinID(ii));
        
        %get yield
        a = false(res_wbm*360,res_wbm*180);
        a(sub2ind(size(a),floor(remfun(res_wbm*basins(kk).X(~isnan(basins(kk).Y)))),floor(res_wbm*(90+basins(kk).Y(~isnan(basins(kk).Y)))))) = 1;
        
        if any(a(1,:)),
            abasin = circshift(imfill(circshift(a,360*res_wbm/2),8,'holes'),-(360*res_wbm/2));
        else,
            abasin = imfill(a,8,'holes');
        end
        cellareabasins = sum(areapercell(abasin(:)));
        
        fact = BasinArea(ii)./cellareabasins;
        
        QRiver_prist(ii) = sum(sum(PristYield.WBMsed(abasin))).*fact;
        QRiver_dist(ii) = sum(sum(DistYield.WBMsed(abasin))).*fact;
        Discharge_prist(ii) = sum(sum(PristYield.WBMdis(abasin))).*fact;
        Discharge_dist(ii) = sum(sum(DistYield.WBMdis(abasin))).*fact;
    end
    close(h)
    
end


%limit number of deltas
idx = Discharge_prist>1;
idx = idx & (QRiver_prist>0.05 | QRiver_dist>0.05 | Discharge_prist>20);

QRiver_prist = max(0,QRiver_prist(idx));
QRiver_dist = max(0,QRiver_dist(idx));
BasinID = BasinID(idx);
BasinArea = BasinArea(idx);
MouthLat = MouthLat(idx);
MouthLon = remfun(MouthLon(idx));
Discharge_dist = Discharge_dist(idx);
Discharge_prist = Discharge_prist(idx);
Continent = Continent(idx);

save([dropbox filesep 'github' filesep 'GlobalDeltaChange' filesep 'GlobalDeltaData.mat'],'Continent','MouthLat','MouthLon','BasinArea','BasinID','QRiver_prist','Discharge_prist','QRiver_dist','Discharge_dist','-append')