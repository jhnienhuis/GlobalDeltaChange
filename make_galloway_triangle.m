
out = load('GlobalDeltaData.mat');
%addpath([dropbox filesep 'WorldDeltas'])

%some key statistics
[~,mor] = max([out.QWave,out.QRiver_prist,out.QTide],[],2);
histcounts(mor)./length(mor)

%out.QTide(out.QTide<=0 | isnan(out.QTide)) = 1;
%out.QWave(out.QWave<=0 | isnan(out.QWave)) = 1;
%out.QRiver_prist(out.QRiver_prist<=0 | isnan(out.QRiver_prist)) = 1;
%out.QRiver_dist(out.QRiver_dist<=0 | isnan(out.QRiver_dist)) = 1;
%out.QWave = out.QWave.*5;
[~,~,x] = xlsread([dropbox filesep 'WorldDeltas' filesep 'FamousDeltaData.xlsx'],'A1:I69');
x_Name = x(3:end,2);
x_CN_prist = [x{3:end,8}];
x_CN_dist = [x{3:end,9}];

[~,idx] = ismember([x{3:end,1}]+1i*[x{3:end,4}],out.BasinID+1i*out.Continent); idx = idx(idx>0);
%table((1:length(idx))',x_Name,out.QRiver_prist(idx),out.QRiver_dist(idx),out.QWave(idx),out.QTide(idx))
idx2 = [1,5,6,12,13,15,16,17,18,19,20,23,25,26,27,29,30,32,33,35,37,40,41,45,49,54,62,67,36];
idx = idx(idx2);
%out.QWave(idx(1)) = 1008;

t = table(idx',x_Name(idx2),out.Discharge_prist(idx),out.QRiver_prist(idx),out.QRiver_dist(idx),out.QWave(idx),out.QTide(idx),out.ee.net_aqua(idx));
%t = table((1:length(idx))',x_Name(idx2),out.QWave(idx),out.QTide(idx));
sortrows(t,2)

%name change plot
figure

[QRiver_prist_log,QWave_prist_log,QTide_prist_log] = DeltaLogMaker(out.QRiver_prist(idx),out.QWave(idx),out.QTide(idx));

[~,x0,y0] = ternplot(QTide_prist_log,QRiver_prist_log,QWave_prist_log,'scatter','filled');

[QRiver_dist_log,QWave_dist_log,QTide_dist_log] = DeltaLogMaker(out.QRiver_dist(idx),out.QWave(idx),out.QTide(idx));

[~,x1,y1] = ternplot(QTide_dist_log,QRiver_dist_log,QWave_dist_log,'scatter');
hold on
[X0, Y0] = ds2nfu(x0, y0);
[X1, Y1] = ds2nfu(x1, y1);

for ii=1:length(x1),
    if distance(X0(ii),Y0(ii),X1(ii), Y1(ii))>0.01,
            annotation('arrow',[X0(ii) X1(ii)]',[Y0(ii) Y1(ii)]')
    else,
        scatter(x1(ii),y1(ii),'.k')
    end
end

text(x1+0.02, y1, x_Name(idx2));
set(gca,'Layer','top')

set([findall(gcf,'String','  25'); findall(gcf,'String','25')],'String','0.1')
set([findall(gcf,'String','  50'); findall(gcf,'String','50')],'String','0.5')
set([findall(gcf,'String','  75'); findall(gcf,'String','75')],'String','0.9')