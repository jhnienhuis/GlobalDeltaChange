function [t] = global_delta_validation

%load prediction

f = [gdrive filesep 'github' filesep 'GlobalDeltaChange' filesep];
load([f 'GlobalDeltaData.mat'],'QRiver_prist','QTide','QWave','BasinID')
load([f 'land_area_change' filesep 'GlobalDeltaData_AreaChange.mat'],'net_aqua','net_pekel');

%load validation
v = load('global_delta_validation.mat');

%mapped delta morphology accuracy
[~,idx] = ismember(v.BasinID,BasinID);

%detecting river mouths
c_rm = [sum(v.rm>0) 2*v.missed_deltas; sum(v.rm==0) inf]./((v.missed_deltas*2)+numel(v.rm));

%global uncertainty in number of deltas
d_unc = 1.96*numel(BasinID)*sqrt(c_rm(1)*(1-c_rm(1))./((v.missed_deltas*2)+numel(v.rm)-1));

%predicted morphology
[~,mor_pred] = max([QWave,QRiver_prist,QTide],[],2);
morphology = ["Wave dominated","River dominated","Tide dominated"]';
d1 = mor_pred(idx(idx~=0));
d2 = v.mor(idx~=0);

c = confusionmat(d1,d2);

tc = table(morphology,c,'VariableNames',["Predicted","Observed"]);



%extra bootstrap
c_rel = c./numel(d2);

c_unc_btstrp = std(bootstrp(1000,@(d1,d2) (sum(diag(confusionmat(d1,d2)/length(idx)))),d1,d2));

%percentage mispredicted, individual prediction accuracy of the classes
c_unc = (diag(c_rel)'./sum(c_rel))-c_unc_btstrp;

%total accuracy of the mapping
s_a = (1-c_unc).*histcounts(mor_pred);

rNames = ["Number of Deltas";morphology];
rvalue = [[numel(BasinID),histcounts(mor_pred)]',[round([d_unc,s_a]')]];    

ta = table(morphology,round(100*c_unc)','VariableNames',["Morphology","Prediction accuracy (%)"]);
tb = table(rNames,rvalue,'VariableNames',["Morphology","Global number +/- Uncertainty"]);

%% adjust readme file

txt = fileread([f 'readme.rst']);
% Split into lines
txt = regexp(txt, '\r\n', 'split');
% Find lines starting with certain letters with confusion matrix
w(1) = find(startsWith(txt, "|           | Wave       |"));
w(2) = find(startsWith(txt, "| Predicted | River      |"));
w(3) = find(startsWith(txt, "|           | Tide       |"));

for ii=1:3,    
    txt(w(ii)) = replace(txt(w(ii)),extract(txt(w(ii)),digitsPattern(3))',cellstr(num2str(c(ii,:)','%03.0f'))');
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

% Write to new file
fid = fopen([f 'readme.rst'], 'wt');
fprintf(fid, '%s\n', txt{:});
fclose(fid);



%% do remote sensing accuracy :-s
v = load('DeltaValidation.mat');

load('GlobalDeltaData.mat')
[~,idx] = ismember(v.BasinID,BasinID);

%number of validation points = 40;
nr = 40;

% 1% error in water vs. sand (see pekel paper) per year
se_1 = 0.01;
% 6% error between the two maps per year
se_2 = var(ee.dep_pekel-ee.dep_aqua)./sqrt(nr);

% error in the delta region mapping
x = v.net(~isnan(v.net) & idx~=0)./31;
y = ee.net_pekel(idx(~isnan(v.net) & idx~=0));

se_3 = var(x-y)./sqrt(numel(x));

%standard error confidence interval of the sum
all = sum(ee.net_aqua)*28
se_all = sqrt(28*(se_1+se_2+se_3))*sqrt(numel(BasinID))

se_all_net = sqrt(2*(sqrt(28*(se_1+se_2+se_3))*sqrt(numel(BasinID))).^2)



tb = table(rNames,rvalue,'VariableNames',["Morphology","Global number +/- Uncertainty"]);

