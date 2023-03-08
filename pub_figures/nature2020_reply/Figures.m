%% Fig. 1
addpath('D:\Drive\github\GlobalDeltaSeaLevel')
v1 = load('GlobalDeltaData_v1.mat','QRiver_dist','QRiver_prist','ee');
load('D:\Drive\github\GlobalDeltaChange\GlobalDeltaData.mat','BasinID2','delta_name','MouthLon');
f = [gdrive filesep 'github' filesep 'GlobalDeltaChange' filesep];
ed_BasinID2 = get_edmonds_data(BasinID2);
ed_idx = ismember(BasinID2,ed_BasinID2);

land = v1.ee.net_aqua;
change = v1.QRiver_dist./v1.QRiver_prist;
changediff = v1.QRiver_dist-v1.QRiver_prist;

[~,~,b] = histcounts(change,[0 0.5 1.5 inf]);

a(1) = subplot(2,3,1);
%all data
[r(1),sp(1)] = corr(change,land,'Type','Spearman');
[rd(1),spd(1)] = corr(changediff,land,'Type','Spearman');
[h(1),p(1)] = kstest2(land(b==1),land(b==3));
m = accumarray(b,land,[],@mean);
s = accumarray(b,land,[],@(x) (std(x)./sqrt(numel(x))));
errorbar([1 3],m([1 3]),s([1 3]),'ob','LineWidth',2,'MarkerFaceColor','b'), hold on
errorbar(2,m(2),s(2),'ob')
title(['(a) all data'])
ylabel('Delta Area Change (km^2/yr)'),

a(2) = subplot(2,3,2);
%using caldwell deltas
[r(2),sp(2)] = corr(change(ed_idx==1),land(ed_idx==1),'Type','Spearman');
[rd(2),spd(2)] = corr(changediff(ed_idx==1),land(ed_idx==1),'Type','Spearman');
[h(2),p(2)] = kstest2(land(b==1 & ed_idx==1),land(b==3& ed_idx==1));
m = accumarray(b(ed_idx),land(ed_idx),[],@mean);
s = accumarray(b(ed_idx),land(ed_idx),[],@(x) (std(x)./sqrt(numel(x))));
errorbar([1 3],m([1 3]),s([1 3]),'ob','LineWidth',2,'MarkerFaceColor','b'), hold on
errorbar(2,m(2),s(2),'ob')
title(['(b) Caldwell data'])


a(3) = subplot(2,3,3);
%excl yellow
[r(3),sp(3)] = corr(change(delta_name~="Yellow"),land(delta_name~="Yellow"),'Type','Spearman');
[rd(3),spd(3)] = corr(changediff(delta_name~="Yellow"),land(delta_name~="Yellow"),'Type','Spearman');
[h(3),p(3)] = kstest2(land(b==1 & delta_name~="Yellow"),land(b==3 & delta_name~="Yellow"));
m = accumarray(b(delta_name~="Yellow"),land(delta_name~="Yellow"),[],@mean);
s = accumarray(b(delta_name~="Yellow"),land(delta_name~="Yellow"),[],@(x) (std(x)./sqrt(numel(x))));
errorbar([1 3],m([1 3]),s([1 3]),'ob','LineWidth',2,'MarkerFaceColor','b'), hold on
errorbar(2,m(2),s(2),'ob')
title(['(c) all data excl. Yellow River'])
set(gcf, 'Units', 'Inches', 'OuterPosition', [0, 0, 7, 3]);

a(4) = subplot(2,3,4);
%smaller buffer size:

fileID = fopen([f 'land_area_change' filesep 'GlobalDeltaChange_p50.csv'],'r');
data = textscan(fileID, '%q%f%f%f%f%f%f%f%f%q%[^\n\r]', 'Delimiter', ',', 'HeaderLines' ,1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
data = cell2mat(data(2:9)); %basin ID, deltaArea, change, deposition, erosion
[~,idx] = ismember(BasinID2,data(:,1));
idx(idx==0) = 2;
land = data(idx,4)/20; %aquamonitor change deposition erosion per year (1984_240_2013_48 => 20 years)
[h(4),p(4)] = kstest2(land(b==1),land(b==3));
[r(4),sp(4)] = corr(change,land,'Type','Spearman');
[rd(4),spd(4)] = corr(changediff,land,'Type','Spearman');
m = accumarray(b,land,[],@mean);
s = accumarray(b,land,[],@(x) (std(x)./sqrt(numel(x))));
errorbar([1 3],m([1 3]),s([1 3]),'ob','LineWidth',2,'MarkerFaceColor','b'), hold on
errorbar(2,m(2),s(2),'ob')
ylabel('Delta Area Change (km^2/yr)'),
title(["(d) 50% buffer size"])


a(5) = subplot(2,3,5);
%larger buffer size:

fileID = fopen([f 'land_area_change' filesep 'GlobalDeltaChange_p150.csv'],'r');
data = textscan(fileID, '%q%f%f%f%f%f%f%f%f%q%[^\n\r]', 'Delimiter', ',', 'HeaderLines' ,1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
data = cell2mat(data(2:9)); %basin ID, deltaArea, change, deposition, erosion
[~,idx] = ismember(BasinID2,data(:,1));
idx(idx==0) = 2;
land = data(idx,4)/20; %aquamonitor change deposition erosion per year (1984_240_2013_48 => 20 years)
[h(5),p(5)] = kstest2(land(b==1),land(b==3));
[r(5),sp(5)] = corr(change,land,'Type','Spearman');
[rd(5),spd(5)] = corr(changediff,land,'Type','Spearman');
m = accumarray(b,land,[],@mean);
s = accumarray(b,land,[],@(x) (std(x)./sqrt(numel(x))));
errorbar([1 3],m([1 3]),s([1 3]),'ob','LineWidth',2,'MarkerFaceColor','b'), hold on
errorbar(2,m(2),s(2),'ob')

title(["(e) 150% buffer size"])
xlabel('Human-induced Sediment Flux Change (Q^d_{river} / Q^p_{river})')

a(6) = subplot(2,3,6);
%altered land cover mask #2:

fileID = fopen([f 'land_area_change' filesep 'GlobalDeltaChange_lc_c_w2.csv'],'r');
data = textscan(fileID, '%q%f%f%f%f%f%f%f%f%q%[^\n\r]', 'Delimiter', ',', 'HeaderLines' ,1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
data = cell2mat(data(2:9)); %basin ID, deltaArea, change, deposition, erosion
[~,idx] = ismember(BasinID2,data(:,1));
idx(idx==0) = 2;
land = data(idx,4)/20; %aquamonitor change deposition erosion per year (1984_240_2013_48 => 20 years)
[h(6),p(6)] = kstest2(land(b==1),land(b==3));
[r(6),sp(6)] = corr(change,land,'Type','Spearman');
[rd(6),spd(6)] = corr(changediff,land,'Type','Spearman');
m = accumarray(b,land,[],@mean);
s = accumarray(b,land,[],@(x) (std(x)./sqrt(numel(x))));
errorbar([1 3],m([1 3]),s([1 3]),'ob','LineWidth',2,'MarkerFaceColor','b'), hold on
errorbar(2,m(2),s(2),'ob')
xlim([0.5 3.5])
title("(f) altered land cover mask")

%a = get(gcf);
for ii=1:length(a),
    a(ii).XTickLabel = {'<-50%','-50% .. +50%','>+50%'};
    a(ii).XLim = [0.5 3.5];
    a(ii).XTick = [1 2 3];
    text(a(ii),1.5,0,[join(["KS-test, p = " num2str(p(ii),'%1.1e')]),...
        join(["Spearmans, p = " num2str(sp(ii),'%1.1e')])])
    
end



saveas(gcf,'Fig1.png')
saveas(gcf,'Fig1.fig')
saveas(gcf,'Fig1.svg')

%% Table S1

table((1:6)',r',sp',rd',spd')


%% Figure S1
clr

v1 = load('GlobalDeltaData_v1.mat','QRiver_dist','QRiver_prist','ee');

land = v1.ee.net_aqua;
change = v1.QRiver_dist./v1.QRiver_prist;

%raw data from paper:
subplot(1,2,1),
b = [-1.1E-2,3.9E-3,7.6E-3,9.9E-3,6.6E-3,5.1E-3,2.0E-4,4.4E-3,3.0E-3,1.0E-2,7.4E-3,1.2E-2,2.8E-2,1.5E-2,2.9E-2,2.5E-2,3.0E-2,2.6E-2,5.6E-2,3.3E-2];
scatter(0.05:0.1:1.95,b,'k','filled'), hold on, grid on
scatter(change,land,1,[0.4 0.4 0.4],'filled')
set(gca,'YLim',[-0.02 0.08]), box on, set(gca,'XLim',[0 2])
plot([0 2],[-0.013 0.037],'r');
title('Published figure')
xlabel('Q_{dist}/Q_{prist}')
ylabel('Land Area Change (km2/yr)')


subplot(1,2,2)

%idx yellow = 6302;

edges = linspace(0,2,21);
%change(6302) = NaN;
[xbin] = discretize(change,edges);
%ybin = accumarray(xbin(~isnan(xbin) & v1.QRiver_prist>10),land(~isnan(xbin) & v1.QRiver_prist>10),[],@(x) (trimmean(x,15)));

ybin = accumarray(xbin(~isnan(xbin)),land(~isnan(xbin)),[],@(x) (mean(x)));
ybin = max(-0.02,min(ybin,0.08));
scatter(edges(2:end)-0.1,ybin,40,'k','filled'), grid on, hold on
%scatter(edges(2)-0.1,-0.02,40,'k','filled'); text(edges(2)-0.02,-0.018,num2str(round(ybin(1),1)));
scatter(change,land,1,[0.4 0.4 0.4],'filled')
set(gca,'YLim',[-0.02 0.08]), box on, set(gca,'XLim',[0 2])
plot([0 2],polyval(polyfit(change(change>edges(1)&change<edges(end)),land(change>edges(1)&change<edges(end)),1),[0 2]),'r');

title('v1 Data')
xlabel('Q_{dist}/Q_{prist}')
ylabel('Land Area Change (km2/yr)')


%saveas(gcf,'FigC1.fig')
%saveas(gcf,'FigC1.svg')


%% Figure S3
clr
addpath('D:\Dropbox\github\GlobalDeltaSeaLevel')
v1 = load('GlobalDeltaData_v1.mat','QRiver_dist','QRiver_prist','ee');
load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','BasinID2','delta_name','MouthLon');

land = v1.ee.net_aqua;
change = v1.QRiver_dist./v1.QRiver_prist;
[~,~,b] = histcounts(change,[0 0.5 1.5 inf]);

a(1) = subplot(1,3,1);
%all data
[r(1),sp(1)] = corr(change,land,'Type','Spearman');
[h(1),p(1)] = kstest2(land(b==1),land(b==3));
m = accumarray(b,land,[],@mean);
s = accumarray(b,land,[],@(x) (std(x)./sqrt(numel(x))));
errorbar([1 3],m([1 3]),s([1 3]),'ob','LineWidth',2,'MarkerFaceColor','b'), hold on
errorbar(2,m(2),s(2),'ob')
title(['(a) v1 dataset'])
ylabel('Delta Area Change (km^2/yr)'),

a(2) = subplot(1,3,2);
%all data incl besset
besset = readtable('D:\Dropbox\2020\2020 GlobalDeltas Nature\reply\besset_delta_change.csv');
[~,besset_idx] = ismember(besset.Deltas,delta_name);

land = v1.ee.net_aqua;
land(besset_idx) = besset.DeltaChange_km2_yr_;

[r(2),sp(2)] = corr(change,land,'Type','Spearman');
[h(2),p(2)] = kstest2(land(b==1),land(b==3));
m = accumarray(b,land,[],@mean);
s = accumarray(b,land,[],@(x) (std(x)./sqrt(numel(x))));
errorbar([1 3],m([1 3]),s([1 3]),'ob','LineWidth',2,'MarkerFaceColor','b'), hold on
errorbar(2,m(2),s(2),'ob')
title(['(b) incl besset et al'])
xlabel('Human-induced Sediment Flux Change (Q^d_{river} / Q^p_{river})')

%v2 dataset
v2 = load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','QRiver_dist','QRiver_prist');
v2.ee = load('D:\Dropbox\github\GlobalDeltaChange\land_area_change\GlobalDeltaData_AreaChange.mat');

land2 = v2.ee.net_aqua;
change2 = v2.QRiver_dist./v2.QRiver_prist;
[~,~,b2] = histcounts(change2,[0 0.5 1.5 inf]);

a(3) = subplot(1,3,3);
[h(3),p(3)] = kstest2(land2(b2==1),land2(b2==3));
[r(3),sp(3)] = corr(change2,land2,'Type','Spearman');
m = accumarray(b2,land2,[],@mean);
s = accumarray(b2,land2,[],@(x) (std(x)./sqrt(numel(x))));
errorbar([1 3],m([1 3]),s([1 3]),'ob','LineWidth',2,'MarkerFaceColor','b'), hold on
errorbar(2,m(2),s(2),'ob')
xlim([0.5 3.5])
title("(c) v2 dataset")

%a = get(gcf);
for ii=1:length(a),
    a(ii).XTickLabel = {'<-50%','-50% .. +50%','>+50%'};
    a(ii).XLim = [0.5 3.5];
    a(ii).XTick = [1 2 3];
    text(a(ii),1.5,0,[join(["KS-test, p = " num2str(p(ii),'%1.1e')]),...
        join(["Spearmans, p = " num2str(sp(ii),'%1.1e')])])
    
end
saveas(gcf,'FigS3.png')
saveas(gcf,'FigS3.fig')

%% Figure S4
clr
addpath('D:\Dropbox\github\GlobalDeltaSeaLevel')
v1 = load('GlobalDeltaData_v1.mat','QRiver_dist','QRiver_prist','ee');
load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','BasinID2','delta_name','MouthLon');

land = v1.ee.net_aqua;
change = v1.QRiver_dist./v1.QRiver_prist;
[~,~,b] = histcounts(change,[0 0.5 1.5 inf]);

a(1) = subplot(1,3,1);
%all data
[r(1),sp(1)] = corr(change,land,'Type','Spearman');
[h(1),p(1)] = kstest2(land(b==1),land(b==3));
m = accumarray(b,land,[],@mean);
s = accumarray(b,land,[],@(x) (std(x)./sqrt(numel(x))));
errorbar([1 3],m([1 3]),s([1 3]),'ob','LineWidth',2,'MarkerFaceColor','b'), hold on
errorbar(2,m(2),s(2),'ob')
title(['(a) v1 dataset'])
ylabel('Delta Area Change (km^2/yr)'),

a(2) = subplot(1,3,2);
%removing 10% largest change
idx = isoutlier(land,'percentiles',[5 95]);
sum(v1.QRiver_prist(idx))./sum(v1.QRiver_prist)

[r(2),sp(2)] = corr(change(idx~=1),land(idx~=1),'Type','Spearman');
[h(2),p(2)] = kstest2(land(b==1 &(idx~=1)),land(b==3&(idx~=1)));
m = accumarray(b(idx~=1),land(idx~=1),[],@mean);
s = accumarray(b(idx~=1),land(idx~=1),[],@(x) (std(x)./sqrt(numel(x))));
errorbar([1 3],m([1 3]),s([1 3]),'ob','LineWidth',2,'MarkerFaceColor','b'), hold on
errorbar(2,m(2),s(2),'ob')
title(['(b) removing 5% lowest&highest change'])
xlabel('Human-induced Sediment Flux Change (Q^d_{river} / Q^p_{river})')


a(3) = subplot(1,3,3);

ll = [0.01:0.01:0.1 0.2:0.1:1];
p = zeros(1000,1);
ps = zeros(1000,1);
pm = zeros(length(ll),1);
psm = zeros(length(ll),1);

for jj=1:length(ll),
    
    for ii=1:1000,
        idx = randsample(length(land),round(ll(jj)*length(land)));
        idxb = false(size(land));
        idxb(idx) = true;
        
        [~,p(ii)] = kstest2(land(b==1 &idxb),land(b==3&idxb));
        %[~,ps(ii)] = corr(change(idxb),land(idxb),'Type','Spearman');
    end
    
    pm(jj) = mean(p);
    %psm(jj) = mean(ps);  
end

plot(ll*100,pm,'-o')
title("(c) expected p-value")
xlabel("Sample Size (% of all potential deltas)")
ylabel("Expected p-value");


%a = get(gcf);
for ii=1:length(a)-1,
    a(ii).XTickLabel = {'<-50%','-50% .. +50%','>+50%'};
    a(ii).XLim = [0.5 3.5];
    a(ii).XTick = [1 2 3];
    text(a(ii),1.5,0,[join(["KS-test, p = " num2str(p(ii),'%1.1e')]),...
        join(["Spearmans, p = " num2str(sp(ii),'%1.1e')])])
    
end
saveas(gcf,'FigS4.png')
saveas(gcf,'FigS4.fig')


