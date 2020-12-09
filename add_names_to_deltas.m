load([dropbox filesep 'WorldDeltas' filesep 'scripts' filesep 'GlobalDeltaData.mat'])

%do slope adjustment of some famous deltas
[~,~,data] = xlsread([dropbox filesep 'WorldDeltas' filesep 'FamousDeltaData.xlsx'],'A1:M68');
[~,idx] = ismember([data{3:end,1}]+imag([data{3:end,4}]),BasinID+imag(Continent));
delta_name = data(3:end,2);
delta_name_id = [data{3:end,1}];
%
if any(idx==0), disp('some deltas not found in xls'); end
ChannelSlope(idx) = [data{3:end,7}];
Discharge_prist(idx) = [data{3:end,5}];
QRiver_prist(idx) = [data{3:end,6}];

QRiver_dist(idx) = [data{3:end,11}];
Discharge_dist(idx) = [data{3:end,10}];
QWave(idx) = [data{3:end,13}];

QWave = QWave.*5;

%local changes to pekel and aquamonitor
[~,idx] = ismember(delta_name_id,BasinID);
src = [2,1,1,1,1,1,2,2,1,2,2,2,2,2,1,2,1,1,2,1,1,2,2,2,2,2,2,1,2,1,1,2,1,2,1,2,2,2,1,2,2,2,1,2,1,2,2,2,2,2,2,1,1,1,2,2,1,2,1,1,2,2,1,2,1,2];

ee.net_aqua(idx(src==1)) = ee.net_pekel(idx(src==1));
ee.ero_aqua(idx(src==1)) = ee.net_aqua(idx(src==1)) - ee.dep_pekel(idx(src==1));

ee.ero_aqua(idx(src==1)) = min(0,ee.ero_aqua(idx(src==1))+rand(sum(src==1),1)-0.5);
ee.net_aqua(idx(src==1)) = ee.dep_aqua(idx(src==1))+ ee.ero_aqua(idx(src==1));

save('GlobalDeltaData.mat','ChannelSlope','QRiver_prist','QRiver_dist','QWave','QTide','ee','Discharge_dist','Discharge_prist','delta_name','delta_name_id','-append');
