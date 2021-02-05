continents = [1:7 5];

%retrieved from FAO aquamaps river names: http://www.fao.org/nr/water/aquamaps/
fname = {'noram_37341','africa_37333','centam_37249','samerica_37330','europe_37253','asia_37331','australia_37252','neareast_37340'};

out = shaperead('D:\Dropbox\WorldDeltas\scripts\GlobalDeltaData.shp');

load('D:\Dropbox\WorldDeltas\scripts\GlobalDeltaData.mat','Continent'); %,'MouthLat','MouthLon'
%scatter(out_mat.MouthLon,out_mat.MouthLat,20,out_mat.Continent,'filled')

maj_name = cell(size(out));
sub_name = cell(size(out));

for kk=1:length(continents),
    
    
    names = shaperead(['D:\OneDrive - Universiteit Utrecht\HydroSheds\river_names\rivers_' fname{kk} '.shp'],'Selector',{@(v1) (v1 == -999),'TOBAS_ID'},'Attributes',{'MAJ_NAME','SUB_NAME'});
        
    f = find(cellfun(@length,maj_name)<2 & Continent==continents(kk))';
    
    for jj=f,
        if mod(jj,100)==1, jj, end
        for ii=1:length(names),
            xi = floor(length(names(ii).X)/2);
            b = inpolygon(names(ii).X(xi),names(ii).Y(xi),out(jj).X,out(jj).Y);
            if b==1, break, end
        end
        
        if b==1,
            maj_name{jj} = names(ii).MAJ_NAME;
            sub_name{jj} = names(ii).SUB_NAME;
        else,
            maj_name{jj} = "";
            sub_name{jj} = "";
        end
        
        
    end
end
maj_name(cellfun(@isempty,maj_name)) = cellstr(strings(sum(cellfun(@isempty,maj_name)),1));
sub_name(cellfun(@isempty,sub_name)) = cellstr(strings(sum(cellfun(@isempty,sub_name)),1));
save GlobalNames sub_name maj_name


%%
clr
load GlobalNames sub_name maj_name delta_name_id delta_name_continent delta_name
out = load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat');
maj_name = string(maj_name);
sub_name = string(sub_name);

%find all unique major names and give those to the largest drainage basin
r = contains(maj_name,"Coast"); 
maj_name(r) = "";

[blub,ida,idc] = unique(maj_name);

for ii=2:length(ida),
    idx = find(idc==ii);
    [~,idy] = max(out.BasinArea(idx));
    sub_name(idx(idy)) = maj_name(idx(idy));
    
end


%first give all deltas subname
out.delta_name = strings(size(sub_name));

[~,ida,idc] = unique(sub_name);

for ii=2:length(ida),
    idx = find(idc==ii);
    [~,idy] = max(out.BasinArea(idx));
    out.delta_name(idx(idy)) = sub_name(idx(idy));
    
end




[~,idx] = ismember([delta_name_id*10 + delta_name_continent'],out.BasinID2);
out.delta_name(idx) = delta_name;
table(idx',delta_name,out.delta_name(idx),int64(out.BasinID(idx)))



%add some additional ones
idx = [1948,2818,1498,3558,408,1338,2088,3628,1548,3658,638,3888,3758];
n = {'Ob','Yenisei','Dvina','Hatanga','Neva','Kunwak','Pur','Olenek','Anadyr','Shuangtaizi','Kuskokwim','Reka Nizhnyaya','Pyasina'};
[~,idx] = ismember(idx,out.BasinID2);
out.delta_name(idx) = n;
out.delta_name(10749) = 'Colville';
%out = rmfield(out,{'delta_name','delta_name_id','delta_name_continent'});
%out.delta_name = delta_name;

t = table(out.delta_name,out.BasinArea,int64(out.BasinID2));
t = sortrows(t,2,'descend');
t(1:100,:);

save GlobalDeltaData -struct out