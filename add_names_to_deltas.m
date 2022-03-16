function add_names_to_deltas
continents = [1:7 5];

%retrieved from FAO aquamaps river names: http://www.fao.org/nr/water/aquamaps/
fname = {'noram_37341','africa_37333','centam_37249','samerica_37330','europe_37253','asia_37331','australia_37252','neareast_37340'};

sname = 'GlobalDeltaBasins';
unzip(['D:\Dropbox\github\GlobalDeltaChange\export_data\' sname '_shp.zip']);
out = shaperead(['D:\Dropbox\github\GlobalDeltaChange\' sname]);
delete([sname '.dbf'],[sname '.shx'],[sname '.shp'])

load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','Continent'); %,'MouthLat','MouthLon'
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

load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','BasinArea','BasinID2');
maj_name = string(maj_name);
sub_name = string(sub_name);

%find all unique major names and give those to the largest drainage basin
r = contains(maj_name,"Coast"); 
maj_name(r) = "";

[blub,ida,idc] = unique(maj_name);

for ii=2:length(ida),
    idx = find(idc==ii);
    [~,idy] = max(BasinArea(idx));
    sub_name(idx(idy)) = maj_name(idx(idy));
    
end


%first give all deltas subname
delta_name = strings(size(sub_name));

[~,ida,idc] = unique(sub_name);

for ii=2:length(ida),
    idx = find(idc==ii);
    [~,idy] = max(BasinArea(idx));
    delta_name(idx(idy)) = sub_name(idx(idy));
    
end

additional_names = ["Mississippi",    "Brazos",    "Colorado, TX",    "Grijalva",    "Colorado, MX",    "Eel, CA",    "Columbia, WA",    "Fraser",    "Tinajones, CO",    "Orinoco",    "Essequibo",    "Amazon",    "Sao Francisco",    "Uruguay",    "Parana",    "Mekong",    "Yellow",    "Yangtze",    "Ganges-Brahmaputra",    "Nile",    "Irriwaddy",    "Fly",    "Schelde",    "Themes",    "Rhine-Meuse",    "Danube",    "Niger",    "Weser",    "Elbe",    "Po",    "Ombrone",    "Rhone",    "Ebro",    "Congo",    "Volta",    "Senegal",    "Ceyan",    "Chao",    "Copper",    "Godavari",    "Homathko",    "Indigirka",    "Indus",    "Klamath",    "Klinaklini",    "Kolyma",    "Krishna",    "Lena",    "Limpopo",    "MacKenzie",    "Magdalena",    "Mahanadi",    "Orange",    "Pechora",    "Pescara",    "Song Hong",    "Squamish",    "Tigris",    "Var",    "Vistula",    "Volga",    "Waipaoa",    "Yana",    "Yukon",    "Zhujiang","Ob",    "Yenisei",    "Dvina",    "Hatanga",    "Neva",    "Kunwak",    "Pur",    "Olenek",    "Anadyr",    "Shuangtaizi",    "Kuskokwim",    "Reka Nizhnyaya",    "Pyasina","Colville","Yellow, US"];
additional_basinid2 = [     4267691,     4281751,     4301321,      995653,     4056471,     3646801,     3173581,     2619721,      214804,      268554,      325964,      540974,      841614,     1187004,     1210764,     4050686,     1005716,     1780166,     2259146,      200412,     3032116,     7126506,     1353513,     1330295,     1248635,     1832055,      908572,     1131775,     1099655,     1873515,     2205845,     2098785,     2433835,     1086772,      811472,      472622,     3173945,     3279946,         548,     2996486,     2155461,        3338,     2180336,     3576671,     2131071,        2678,     3117246,        3648,     1548922,        2528,       94474,     2763826,     1596892,        2178,     2224895,     2765116,     2457481,     3355105,     2045595,      962135,     1743925,     1766287,        3318,        1258,     2372006,1948,2818,1498,3558,408,1338,2088,3628,1548,3658,638,3888,3758,2938,4076891];

%n = {'Ob','Yenisei','Dvina','Hatanga','Neva','Kunwak','Pur','Olenek','Anadyr','Shuangtaizi','Kuskokwim','Reka Nizhnyaya','Pyasina'};
[~,idx] = ismember(additional_basinid2,BasinID2);
delta_name(idx) = additional_names;
%out.delta_name(10749) = 'Colville';
%out = rmfield(out,{'delta_name','delta_name_id','delta_name_continent'});
%out.delta_name = delta_name;

t = table(delta_name,BasinArea,int64(BasinID2));
t = sortrows(t,2,'descend');
t(1:100,:);

save GlobalDeltaData delta_name -append