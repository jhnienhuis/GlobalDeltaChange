load('GlobalDeltaData.mat')

[~,mor_pred] = max([QWave,QRiver_prist,QTide],[],2);
icon(mor_pred == 3) = {'http://maps.google.com/mapfiles/kml/paddle/red-blank.png'};
icon(mor_pred == 2) = {'http://maps.google.com/mapfiles/kml/paddle/grn-blank.png'};
icon(mor_pred == 1) = {'http://maps.google.com/mapfiles/kml/paddle/blu-blank.png'};
[~,ic] = ismember(delta_name_id,BasinID);
name = cell(size(BasinID));
name(ic) = delta_name;

blub = num2cell([BasinID QRiver_prist QRiver_dist QWave QTide])';

for ii=1:length(blub),
    desc{ii} = sprintf('BasinID: %1.0i</br>QRiver_prist: %1.0e kg/s</br>QRiver_dist: %1.0e kg/s</br>QWave: %1.0e kg/s</br>QTide: %1.0e kg/s</br>',blub{:,ii});
end

kmlwritepoint('GlobalDeltaData.kml',MouthLat,MouthLon,'name',name,'icon',icon,'description',desc)
