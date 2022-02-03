function create_xls

load('D:\Drive\github\GlobalDeltaChange\land_area_change\GlobalDeltaData_AreaChange.mat','net_aqua')
load('D:\Drive\github\GlobalDeltaChange\GlobalDeltaData.mat','BasinID2','delta_name','MouthLon','MouthLat','BasinArea','Discharge_prist','QRiver_prist','QRiver_dist','QTide','QWave');


t = table(BasinID2,delta_name,MouthLon,MouthLat,BasinArea,Discharge_prist,QRiver_prist,QRiver_dist,QTide,QWave,net_aqua);

writetable(t,'GlobalDeltaData.xlsx')
end