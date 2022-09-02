%arctic deltas
%this is from: https://github.com/jhnienhuis/GlobalDeltaChange
%more information is here: https://www.nature.com/articles/s41586-019-1905-9
d = load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','QRiver_dist','Continent','QRiver_prist');

arctic = d.Continent==8;
ee = load('D:\Dropbox\github\GlobalDeltaChange\land_area_change\GlobalDeltaData_AreaChange.mat');

% delta area change
%median delta gross changes (dep+ero) in m2/yr for ~arctic and arctic deltas
median(ee.dep_aqua(~arctic & d.QRiver_dist>1)-ee.ero_aqua(~arctic & d.QRiver_dist>1))*1e6
median(ee.dep_aqua(arctic & d.QRiver_dist>1)-ee.ero_aqua(arctic & d.QRiver_dist>1))*1e6

%mean net changes for ~arctc and non arctic deltas
mean(ee.net_aqua(~arctic & d.QRiver_dist>1))*1e6
mean(ee.net_aqua(arctic & d.QRiver_dist>1))*1e6

%efficiency of land gain for growing arctic deltas. comparing land gain
%(m2/yr) to fluvial sediment flux (from kg/s to m3/yr assuming ~1600kg/m3 bulk density)
1e6.*sum(ee.net_aqua(~arctic &ee.net_aqua>0))/sum(365*24*3600./1600.*d.QRiver_prist(~arctic &ee.net_aqua>0))
1e6.*sum(ee.net_aqua(arctic &ee.net_aqua>0))/sum(365*24*3600./1600.*d.QRiver_prist(arctic &ee.net_aqua>0))
