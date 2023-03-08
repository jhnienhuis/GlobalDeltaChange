function [QRiverLog,QWaveLog,QTideLog] = DeltaLogMaker(QRiver,QWave,QTide)

ter_map_lin = linspace(0,1,1001);
map_own = [0 0.03 0.1 0.2 0.5 0.80 0.9 0.97 1];
%which is approximately this: 0.5+ 0.2T -0.1T*log10(1-T^2) where T=(2x-1)
ter_map = interp1(linspace(0,1,length(map_own)),map_own,ter_map_lin,'pchip');


QTotal = QTide+QRiver+QWave;
QWaveLog = interp1(ter_map,ter_map_lin,QWave./QTotal,'spline');
QRiverLog = interp1(ter_map,ter_map_lin,QRiver./QTotal,'spline');
QTideLog = interp1(ter_map,ter_map_lin,QTide./QTotal,'spline');

end
