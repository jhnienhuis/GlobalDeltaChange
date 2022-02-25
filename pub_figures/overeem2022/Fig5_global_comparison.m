function Fig5_global_comparison

%this is from: https://github.com/jhnienhuis/GlobalDeltaChange
d = load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat');

%this if from: https://github.com/jhnienhuis/GlobalDeltaChange/tree/master/land_area_change
ee = load('D:\Dropbox\github\GlobalDeltaChange\land_area_change\GlobalDeltaData_AreaChange.mat');

%this is from: https://github.com/jhnienhuis/GlobalDeltaSeaLevel/tree/master/export_data
s = load('D:\Dropbox\github\GlobalDeltaSeaLevel\export_data\GlobalDeltaSeaLevelData.mat','DeltaSLR');

p = load('Fig5_DischargePeakiness.mat');

arctic = d.Continent==8;

subplot(2,5,1)
histogram(log10(abs(d.shelf_slope(~arctic & d.QRiver_dist>1))),[-4:0.25:0],'Normalization','probability'), hold on
histogram(log10(abs(d.shelf_slope(arctic & d.QRiver_dist>1))),[-4:0.25:0],'Normalization','probability')
set(gca,'XTick',[-4:0],'XTickLabels',{'10-4','10-3','10-2','10-1','100'})
xlabel('shelf slope (m/m)')
legend('Temperate Deltas','Arctic Deltas')
set(gca,'YLim',[0 0.25],'YTick',[0 0.1 0.2])

subplot(2,5,2)
gro = 1e6.*(ee.dep_aqua-ee.ero_aqua);
histogram(log10(abs(gro(~arctic  & d.QRiver_dist>1))),[0:0.5:7],'Normalization','probability'), hold on
histogram(log10(abs(gro(arctic  & d.QRiver_dist>1))),[0:0.5:7],'Normalization','probability')
set(gca,'XTick',[0:2:7],'XTickLabels',{'100','102','104','106'})
xlabel('gross changes (m2yr-1)')
set(gca,'YLim',[0 0.25],'YTick',[0 0.1 0.2])

subplot(2,5,7)
eff = 1e6.*ee.net_aqua./(365*24*3600./1600.*d.QRiver_prist);
histogram(log10(abs(eff(~arctic &ee.net_aqua>0))),[-6:0.5:2],'Normalization','probability'), hold on
histogram(log10(abs(eff(arctic &ee.net_aqua>0))),[-6:0.5:2],'Normalization','probability')
set(gca,'XTick',[-6:2:2],'XTickLabels',{'10-6','10-4','10-2','100','102'})
xlabel('eff (m-1)')
set(gca,'YLim',[0 0.25],'YTick',[0 0.1 0.2])

subplot(2,5,6)

histogram(s.DeltaSLR(~arctic).*1e3,[-15:1:5],'Normalization','probability'), hold on
histogram(s.DeltaSLR(arctic).*1e3,[-15:1:5],'Normalization','probability')
xlabel('RSLR (mm/yr)')
set(gca,'YLim',[0 0.5],'YTick',[0 0.2 0.4])

subplot(2,5,3)
R = d.QRiver_dist./d.QWave;
histogram(log10(abs(R(~arctic & d.QRiver_dist>1))),[-4:0.5:4],'Normalization','probability'), hold on
histogram(log10(abs(R(arctic & d.QRiver_dist>1))),[-4:0.5:4],'Normalization','probability')
xlabel('QRiver / QWave')
set(gca,'XTick',[-4:2:4],'XTickLabels',{'10-4','10-2','100','102','104'})
set(gca,'YLim',[0 0.25],'YTick',[0 0.1 0.2])

subplot(2,5,8)
T = d.Discharge_tide./d.Discharge_prist;
histogram(log10(abs(T(~arctic & d.QRiver_dist>1))),[-4:0.5:4],'Normalization','probability'), hold on
histogram(log10(abs(T(arctic & d.QRiver_dist>1))),[-4:0.5:4],'Normalization','probability')
xlabel('QTide / QRiver')
set(gca,'XTick',[-4:2:4],'XTickLabels',{'10-4','10-2','100','102','104'})
set(gca,'YLim',[0 0.25],'YTick',[0 0.1 0.2])

subplot(2,5,4)
p90p10 = p.p90./p.p10;
histogram(p90p10(p.rm_lat<60 & p.m>2),[0:10:200 max(p90p10)],'Normalization','probability'), hold on
histogram(p90p10(p.rm_lat>60 & p.m>2),[0:10:200 max(p90p10)],'Normalization','probability')
xlabel('p90 / p10')
set(gca,'XTick',[0:50:200])
set(gca,'YLim',[0 0.5],'YTick',[0 0.2 0.4])

%
% make triangle
[QRiver_prist_log,QWave_prist_log,QTide_prist_log] = DeltaLogMaker(d.QRiver_prist(arctic),d.QWave(arctic),d.QTide(arctic));

subplot(2,5,[5 10])

%cbrewer can be found at the matlab fileexchange: https://nl.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab
colormap(flipud(cbrewer('div', 'RdYlBu', 64)))

%ternplot can be found at the matlab fileexchange: https://nl.mathworks.com/matlabcentral/fileexchange/2299-alchemyst-ternplot
[~,x0,y0] = ternplot(QTide_prist_log,QRiver_prist_log,QWave_prist_log,'scatter',3+(d.QRiver_prist(arctic)).^0.7,log10(d.QRiver_prist(arctic)),'filled');

%add names
d.delta_name(d.BasinID2==2728) = "Jago";
arc_name = ["Lena","Colville","MacKenzie","Yukon","Jago","Ob","Yenisei","Kolyma"];
[~,idx_name] = ismember(arc_name,d.delta_name(arctic));
text(x0(idx_name),y0(idx_name),arc_name)


caxis([0,4])
set(gca,'Layer','top')
set([findall(gcf,'String','  25'); findall(gcf,'String','25')],'String','0.1')
set([findall(gcf,'String','  50'); findall(gcf,'String','50')],'String','0.5')
set([findall(gcf,'String','  75'); findall(gcf,'String','75')],'String','0.9')
h = colorbar('location','east');
set(h,'YTick',0:5,'YTickLabel',cellstr(num2str((10.^(0:5))','%1.0f'))); ylabel(h,'QRiver (kgs-1)')

%
set(findall(gcf,'type','Axes'), 'FontSize', 7)
set(findall(gcf,'Type','text'), 'FontSize', 7)
set(gcf, 'Units', 'Centimeters', 'OuterPosition', [0, 0, 18.3, 10]);
saveas(gcf,'Fig5_global_comparison.svg')
end


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