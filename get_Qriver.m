%% get flux per river mouth
continents = {'na','af','ca','sa','eu','as','au'};

%resolution is 240 per degree (15s)
res = 240;
%resolution WBMSed is 10 per degree (6 min)
res_wbm = 10;
remfun = @(lon) (rem(res_wbm*360-1+lon,res_wbm*360)+1);

PristYield = load('D:\GlobalDatasets\WBMSed\Prist_Yield.mat');
DistYield = load('D:\GlobalDatasets\WBMSed\Dist_Yield.mat');

%area per cell of WBMSed in km2
areapercell = repmat(6371.^2.*2*pi/360/res_wbm*(sin(deg2rad(-89.9:0.1:90))-sin(deg2rad(-90:0.1:89.9))),360*res_wbm,1);

for jj=1:length(continents),
    
    %load river mouth info
    load([dropbox filesep 'WorldDeltas\scripts\Rivermouth' continents{jj} '.mat'])
    jj
    %get all north american river basins > 50km2
    basins = shaperead(['D:\GlobalDatasets\HydroSheds\' continents{jj} '_bas_15s_beta\' continents{jj} '_bas_15s_beta.shp'],'Selector',{@(v1) (v1>40),'AREA_SQKM'});
    disp('Loaded basins...')
    basins_list = [basins(:).BASIN_ID];
    
    %preallocate flux vectors
    cellareabasins = zeros(length(BasinID),1);
    QRiver_prist = zeros(length(BasinID),1);
    Discharge_prist = zeros(length(BasinID),1);
    QRiver_dist = zeros(length(BasinID),1);
    Discharge_dist = zeros(length(BasinID),1);
    h=waitbar(0);
    for ii=1:length(BasinID),
        
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
        cellareabasins(ii) = sum(areapercell(abasin(:)));
        
        QRiver_prist(ii) = sum(sum(PristYield.WBMsed(abasin)));
        QRiver_dist(ii) = sum(sum(DistYield.WBMsed(abasin)));
        Discharge_prist(ii) = sum(sum(PristYield.WBMdis(abasin)));
        Discharge_dist(ii) = sum(sum(DistYield.WBMdis(abasin)));
    end
    close(h)
    save([dropbox filesep 'WorldDeltas\scripts\Rivermouth' continents{jj} '.mat'],'cellareabasins','QRiver_prist','Discharge_prist','QRiver_dist','Discharge_dist','-append')
end
