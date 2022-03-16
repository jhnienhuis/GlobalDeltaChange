function [t] = global_delta_validation

%load prediction

f = [gdrive filesep 'github' filesep 'GlobalDeltaChange' filesep];
load([f 'GlobalDeltaData.mat'],'QRiver_prist','QTide','QWave','BasinID2','delta_name')

%load validation
v = readtable('global_delta_validation.xlsx');

%mapped delta morphology accuracy
[~,idx] = ismember(v.BasinID2,BasinID2);

%detecting river mouths
c_rm = [sum(v.rm>0) v.missed_deltas(1)*2; sum(v.rm==0) inf]./((v.missed_deltas(1).*2)+numel(v.rm));

%global uncertainty in number of deltas
d_unc = 1.96*numel(BasinID2)*sqrt(c_rm(1)*(1-c_rm(1))./((v.missed_deltas(1)*2)+numel(v.rm)-1));

%predicted morphology
[~,mor_pred] = max([QWave,QRiver_prist,QTide],[],2);
morphology = ["Wave dominated","River dominated","Tide dominated"]';
d1 = mor_pred(idx);
d2 = v.morphology;

c = confusionmat(d1,d2);

tc = table(morphology,c,'VariableNames',["Predicted","Observed"]);

%extra bootstrap
c_rel = c./numel(d2);

c_unc_btstrp = std(bootstrp(1000,@(d1,d2) (sum(diag(confusionmat(d1,d2)/length(idx)))),d1,d2));

%percentage mispredicted, individual prediction accuracy of the classes
c_unc = (diag(c_rel)'./sum(c_rel))-c_unc_btstrp;

%total accuracy of the mapping
s_a = (1-c_unc).*histcounts(mor_pred);

%rNames = ["Number of Deltas";morphology];
rvalue = [[numel(BasinID2),histcounts(mor_pred)]',[round([d_unc,s_a]')]];    

%ta = table(morphology,round(100*c_unc)','VariableNames',["Morphology","Prediction accuracy (%)"]);
%tb = table(rNames,rvalue,'VariableNames',["Morphology","Global number +/- Uncertainty"]);



%% do remote sensing accuracy :-s
f = [gdrive filesep 'github' filesep 'GlobalDeltaChange' filesep];
load([f 'land_area_change' filesep 'GlobalDeltaData_AreaChange.mat'],'delta_name','BasinID2','net_aqua','net_pekel');

v = readtable('global_delta_validation.xlsx');
v.land_area_change = cellfun(@str2num,v.land_area_change);
[~,idx] = ismember(v.BasinID2(~isnan(v.land_area_change)),BasinID2);

%relative mean absolute errors (MAE), as fraction of mean change

% 1% error in water vs. sand (see pekel paper) per year
se(1) = 0.01;
mae(1) = 0.01;

% 6% error between the two maps per year
se(2) = var(net_pekel-net_aqua)./sqrt(numel(net_aqua));
mae(2) = mean(abs(net_pekel-net_aqua))./mean(abs(net_aqua))

% error in the delta region mapping
x = v.land_area_change(~isnan(v.land_area_change));
y = net_aqua(idx);

se(3) = var(x-y)./sqrt(numel(x));
mae(3) = mean(abs(x-y))./mean(abs(y));

%standard error confidence interval of the sum
all = sum(net_aqua);
se_all = sqrt(sum(se).*numel(BasinID2));

unc_matrix = [[100.*mae' mae'.*sum(net_aqua)];[sum(100*mae) sum(mae)'.*mean(net_aqua)];[100.*se_all./all se_all]];

%% adjust readme file

txt = fileread([f 'readme.rst']);
% Split into lines
txt = string(regexp(txt, '\r\n', 'split'));
% Find lines starting with certain letters with confusion matrix
w(1) = find(startsWith(txt, "|           | Wave       |"));
w(2) = find(startsWith(txt, "| Predicted | River      |"));
w(3) = find(startsWith(txt, "|           | Tide       |"));

for ii=1:3,
    str = strfind(txt(w(ii)),digitsPattern(3));
    for jj=1:3,
        txt(w(ii)) = replaceBetween(txt(w(ii)),str(jj),str(jj)+2,string(num2str(c(ii,jj)','%03.0f')));
    end
    
end

% Find lines starting with certain letters with individual predictions
w = find(contains(txt,"Prediction accuracy (%)"));
w(1:3) = w+[2 3 4];

for ii=1:3,    
    txt(w(ii)) = replace(txt(w(ii)),extract(txt(w(ii)),digitsPattern(2))',cellstr(num2str(round(100*c_unc(ii)),'%02.0f'))');
end

%add global uncertainty of morphology
w = find(contains(txt,"Uncertainty (+/- 1std)"));
w = w+[2 3 4 5];

for ii=1:4,    
    txt(w(ii)) = replace(txt(w(ii)),extract(txt(w(ii)),digitsPattern(4,5))',[{num2str(rvalue(ii,1),'%05.0f')},{num2str(rvalue(ii,2),'%04.0f')}]);
end

%add global uncertainty of morphology
w = find(contains(txt,"Selection            Percentage of      Expressed in "));
w = w+[3 4 5 7 8];

for ii=1:length(w),    
    txt(w(ii)) = replace(txt(w(ii)),extract(txt(w(ii)),digitsPattern(3)+"%"),strcat(num2str(unc_matrix(ii,1),'%03.0f'), "%"));
    txt(w(ii)) = replace(txt(w(ii)),extract(txt(w(ii)),digitsPattern(3)+"."+digitsPattern(2)),num2str(unc_matrix(ii,1),'%06.2f'));
end




% Write to new file
fid = fopen([f 'readme.rst'], 'wt');
fprintf(fid, '%s\n', txt{:});
fclose(fid);





