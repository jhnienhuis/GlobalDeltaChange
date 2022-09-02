%% Figures Paper Delta Morphodynamics
% Figure 9 - Prediction delta morphology from field-based data

clr

load('D:\drive\github\GlobalDeltaChange\GlobalDeltaData.mat')

% Find indices within latitudinal extents
indexLat_PacificCol = find(round(MouthLat*10000)/10000>=1.6417...
    & round(MouthLat*1000)/1000<=7.108);
indexLat_CaribbeanCol = find(round(MouthLat*1000)/1000>=7.929...
    & round(MouthLat*100)/100<=11.77);
% Find indices within longitudinal extents that already fit latitudes
index_PacificCol0 = indexLat_PacificCol(round(MouthLon(indexLat_PacificCol)*10)/10>=281 ...
    & round(MouthLon(indexLat_PacificCol)*10)/10<=282.9);
index_CaribbeanCol0 = indexLat_CaribbeanCol(round(MouthLon(indexLat_CaribbeanCol)*10)/10>=282.7 ...
    & round(MouthLon(indexLat_CaribbeanCol)*10)/10<=287.5);

MouthLat_Col = [MouthLat(index_PacificCol0); MouthLat(index_CaribbeanCol0)];
MouthLon_Col = [MouthLon(index_PacificCol0); MouthLon(index_CaribbeanCol0)];

% Find Basin IDs from latitudinal-longitudinal indices
basinID_CaribbeanCol = BasinID(index_CaribbeanCol0);
basinID_PacificCol = BasinID(index_PacificCol0);

% World deltas
% Galloway (1975, Fig. 3)
index_Mississippi = find(BasinID==426769); % 100
index_SaoFrancisco = find(BasinID==84161); % 101
index_Copper = find(BasinID==54); % 102
index_Danube = find(BasinID==183205); % 103
index_Mahakam = find(BasinID==549139); % 104
index_Fly = find(BasinID==712650); % 105

% Nienhuis et al. (2020, Fig. 2)
index_Amazon = find(BasinID==54097); % 106
index_Eel = find(BasinID==364680); % 107
index_Yangtze = find(BasinID==178016); % 108
index_Columbia = find(BasinID==317358); % 109
index_Elbe = find(BasinID==109965); % 110
index_Grijalva = find(BasinID==99565); % 111
index_Lena = find(BasinID==364); % 112
index_Nile = find(BasinID==20041); % 113

% Best (2019)
index_GangesBrahmaputra = find(BasinID==225914); % 114
index_Irrawaddy = find(BasinID==303211); % 115
index_Parana = find(BasinID==121076); % 116
index_Orinoco = find(BasinID==26855); % 117
index_Ob = find(BasinID==194); % 118
index_Yenisei = find(BasinID==281); % 119, it has two indices, chosen the second
     index_Yenisei = index_Yenisei(2,1);
index_MacKenzie = find(BasinID==252); % 120
index_Congo = find(BasinID==108677); % 121 
index_Niger = find(BasinID==90857); % 122 
index_Zambezi = find(BasinID==141681); % 123
index_Limpopo = find(BasinID==154892); % 124


% Indices within compiled dataset (Pacific)
% South to North
index_Mira = find(BasinID==41423);
index_Patia = find(BasinID==39739); % At Sanquianga
index_SanJuan = find(BasinID==37364); % At Boca San Juan
index_Atrato = find(BasinID==29299); % Atrato at Boca El Roto
index_Magdalena = find(BasinID==9447); % Magdalena at Bocas de Ceniza
index_Sinu = find(BasinID==21480); % Sinu at Tinajones




index_NorthAmerica = [index_Mississippi; index_Grijalva; index_Eel; ...
    index_Columbia; index_Copper; index_MacKenzie];

index_SouthAmerica = [index_Amazon; index_Orinoco; index_SaoFrancisco; ...
    index_Parana; index_Magdalena; index_Sinu; index_Atrato; index_SanJuan; ...
    index_Patia; index_Mira];

index_AfricaEurope = [index_Congo; index_Niger; index_Nile; index_Zambezi; ...
    index_Limpopo; index_Elbe; index_Danube];

index_AsiaOceania = [index_Irrawaddy; index_GangesBrahmaputra; index_Ob; ...
    index_Yenisei; index_Lena; index_Yangtze; index_Mahakam; index_Fly];




index_total = [index_NorthAmerica; index_SouthAmerica; index_AfricaEurope; index_AsiaOceania];
    index_WorldEx_total_NA = (1:length(index_NorthAmerica))';
    index_WorldEx_total_SA = (1:length(index_SouthAmerica))' + index_WorldEx_total_NA(end);
    index_WorldEx_total_AE = (1:length(index_AfricaEurope))' + index_WorldEx_total_SA(end);
    index_WorldEx_total_AO = (1:length(index_AsiaOceania))' + index_WorldEx_total_AE(end);

    
    
Deltanamefile_NorthAmerica = {'Mississippi';...
    'Grijalva';...
    'Eel';...
    'Columbia';...
    'Copper';...
    'MacKenzie'};

Deltanamefile_SouthAmerica = {'Amazon';...
    'Orinoco';...
    'SaoFrancisco';...
    'Parana';...
    'Magdalena';...
    'Sinu';...
    'Atrato';...
    'SanJuan';...
    'Patia';...
    'Mira'};

Deltanamefile_AfricaEurope = {'Congo';...
    'Niger';...
    'Nile';...
    'Zambezi';...
    'Limpopo';...
    'Elbe';...
    'Danube'};

Deltanamefile_AsiaOceania = {'Irrawaddy';...
    'GangesBrahmaputra';...
    'Ob';...
    'Yenisei';...
    'Lena';...
    'Yangtze';...
    'Mahakam';...
    'Fly'}; 



Deltaname_NorthAmerica = {'Mississippi';...
    'Grijalva';...
    'Eel';...
    'Columbia';...
    'Copper';...
    'MacKenzie'};

Deltaname_SouthAmerica = {'Amazon';...
    'Orinoco';...
    'Sao Francisco';...
    'Parana';...
    'Magdalena';...
    'Sinu';...
    'Atrato';...
    'San Juan';...
    'Patia';...
    'Mira'};

Deltaname_AfricaEurope = {'Congo';...
    'Niger';...
    'Nile';...
    'Zambezi';...
    'Limpopo';...
    'Elbe';...
    'Danube'};

Deltaname_AsiaOceania = {'Irrawaddy';...
    'Ganges-Brahmaputra';...
    'Ob';...
    'Yenisei';...
    'Lena';...
    'Yangtze';...
    'Mahakam';...
    'Fly'};











% North America
QRiver_NorthAmerica = QRiver_prist(index_NorthAmerica);
qRiver_NorthAmerica = Discharge_prist(index_NorthAmerica);
QTide_NorthAmerica = QTide(index_NorthAmerica);
QWave_NorthAmerica = QWave(index_NorthAmerica);
    QWave_NorthAmerica(isnan(QWave_NorthAmerica)==1) = 0;

Sum_Q_NorthAmerica = QRiver_NorthAmerica...
    + QWave_NorthAmerica...
    + QTide_NorthAmerica;


% South America
QRiver_SouthAmerica = QRiver_prist(index_SouthAmerica);
qRiver_SouthAmerica = Discharge_prist(index_SouthAmerica);
QTide_SouthAmerica = QTide(index_SouthAmerica);
QWave_SouthAmerica = QWave(index_SouthAmerica);
    QWave_SouthAmerica(isnan(QWave_SouthAmerica)==1) = 0;

Sum_Q_SouthAmerica = QRiver_SouthAmerica...
    + QWave_SouthAmerica...
    + QTide_SouthAmerica;


% Africa and Europe
QRiver_AfricaEurope = QRiver_prist(index_AfricaEurope);
qRiver_AfricaEurope = Discharge_prist(index_AfricaEurope);
QTide_AfricaEurope = QTide(index_AfricaEurope);
QWave_AfricaEurope = QWave(index_AfricaEurope);
    QWave_AfricaEurope(isnan(QWave_AfricaEurope)==1) = 0;

Sum_Q_AfricaEurope = QRiver_AfricaEurope...
    + QWave_AfricaEurope...
    + QTide_AfricaEurope;


% Asia and Oceania
QRiver_AsiaOceania = QRiver_prist(index_AsiaOceania);
qRiver_AsiaOceania = Discharge_prist(index_AsiaOceania);
QTide_AsiaOceania = QTide(index_AsiaOceania);
QWave_AsiaOceania = QWave(index_AsiaOceania);
    QWave_AsiaOceania(isnan(QWave_AsiaOceania)==1) = 0;

Sum_Q_AsiaOceania = QRiver_AsiaOceania...
    + QWave_AsiaOceania...
    + QTide_AsiaOceania;






% r_x values
rRiver_NorthAmerica = QRiver_NorthAmerica./Sum_Q_NorthAmerica;
rWave_NorthAmerica = QWave_NorthAmerica./Sum_Q_NorthAmerica;
rTide_NorthAmerica = QTide_NorthAmerica./Sum_Q_NorthAmerica;
rRiver_SouthAmerica = QRiver_SouthAmerica./Sum_Q_SouthAmerica;
rWave_SouthAmerica = QWave_SouthAmerica./Sum_Q_SouthAmerica;
rTide_SouthAmerica = QTide_SouthAmerica./Sum_Q_SouthAmerica;
rRiver_AfricaEurope = QRiver_AfricaEurope./Sum_Q_AfricaEurope;
rWave_AfricaEurope = QWave_AfricaEurope./Sum_Q_AfricaEurope;
rTide_AfricaEurope = QTide_AfricaEurope./Sum_Q_AfricaEurope;
rRiver_AsiaOceania = QRiver_AsiaOceania./Sum_Q_AsiaOceania;
rWave_AsiaOceania = QWave_AsiaOceania./Sum_Q_AsiaOceania;
rTide_AsiaOceania = QTide_AsiaOceania./Sum_Q_AsiaOceania;

[rRiver_NorthAmerica_log,rWave_NorthAmerica_log,rTide_NorthAmerica_log] = DeltaLogMaker(rRiver_NorthAmerica, rWave_NorthAmerica, rTide_NorthAmerica);
[rRiver_SouthAmerica_log,rWave_SouthAmerica_log,rTide_SouthAmerica_log] = DeltaLogMaker(rRiver_SouthAmerica, rWave_SouthAmerica, rTide_SouthAmerica);
[rRiver_AfricaEurope_log,rWave_AfricaEurope_log,rTide_AfricaEurope_log] = DeltaLogMaker(rRiver_AfricaEurope, rWave_AfricaEurope, rTide_AfricaEurope);
[rRiver_AsiaOceania_log,rWave_AsiaOceania_log,rTide_AsiaOceania_log] = DeltaLogMaker(rRiver_AsiaOceania, rWave_AsiaOceania, rTide_AsiaOceania);












% Field-based ratios
% Multiple distributaries (e.g., Patia) coalesce to a single delta, thus
% the Qriver must be modified to account for this
GEdatapath = [cd filesep 'data' filesep];
WW3datapath = [cd filesep 'data' filesep];
Tidesdatapath = [cd filesep 'data' filesep];


fileNameGE_NorthAmerica = {'Wor_Mississippi';...1
    'Wor_Grijalva';...12
    'Wor_Eel';...8
    'Wor_Columbia';...10
    'Wor_Copper';...3
    'Wor_MacKenzie'};% 21

fileNameGE_SouthAmerica = {'Wor_Amazon';...7
    'Wor_Orinoco';...18
    'Wor_SaoFrancisco';...2
    'Wor_Parana';...17
    'Car_MagdalenaAtBocasDeCeniza';...26
    'Car_SinuAtTinajones'
    'Car_Atrato';...27
    'Pac_SanJuan';...29
    'Pac_Patia';...30
    'Pac_Mira'};% 

fileNameGE_AfricaEurope = {'Wor_Congo';...22
    'Wor_Niger';...23
    'Wor_Nile';...
    'Wor_Zambezi';...24
    'Wor_Limpopo';...
    'Wor_Elbe';...11
    'Wor_Danube'};%4
    
fileNameGE_AsiaOceania = {'Wor_Irrawaddy';...16
    'Wor_GBM';...15
    'Wor_Ob';...19
    'Wor_Yenisei';...20
    'Wor_Lena';...13
    'Wor_Yangtze';...9
    'Wor_Mahakam';...5
    'Wor_Fly'};% 6    






fileNameWW3_NorthAmerica = {'WaveWatch_timeseries_213_Wor_Mississippi.nc';...15
    'WaveWatch_timeseries_15382_Wor_Grijalva.nc';...26
    'WaveWatch_timeseries_9221_Wor_Eel.nc';...22
    'WaveWatch_timeseries_9252_Wor_Columbia.nc';...24
    'WaveWatch_timeseries_14531_Wor_Copper.nc';...17
	'WaveWatch_timeseries_14650_Wor_MacKenzie.nc'}; %35

fileNameWW3_SouthAmerica = {'WaveWatch_timeseries_17090_Wor_Amazon.nc';...21
    'WaveWatch_timeseries_16854_Wor_Orinoco.nc';...32
    'WaveWatch_timeseries_17250_Wor_SaoFrancisco.nc';...16
    'WaveWatch_timeseries_16925_Wor_Parana.nc';...31
    'WaveWatch_timeseries_16179_Car_MagdalenaAtBocasDeCenizaCO.nc';...4    
    'WaveWatch_timeseries_16130_Car_SinuAtTinajonesCO.nc';... % 7
    'WaveWatch_timeseries_16085_Car_AtratoCO.nc';...1
    'WaveWatch_timeseries_16056_Pac_SanJuanCO.nc';...
    'WaveWatch_timeseries_15982_Pac_PatiaCO.nc';...12
    'WaveWatch_timeseries_15943_Pac_MiraCO.nc'};

fileNameWW3_AfricaEurope = {'WaveWatch_timeseries_9668_Wor_Congo.nc';...36
	'WaveWatch_timeseries_9597_Wor_Niger.nc';...37
    'WaveWatch_timeseries_2909_Wor_Nile.nc';...
	'WaveWatch_timeseries_7662_Wor_Zambezi.nc';...38
	'WaveWatch_timeseries_9972_Wor_Limpopo.nc';...
    'WaveWatch_timeseries_5506_Wor_Elbe.nc';...25
    'WaveWatch_timeseries_2882_Wor_Danube.nc'}; %18
    
fileNameWW3_AsiaOceania = {'WaveWatch_timeseries_10833_Wor_Irrawaddy.nc';...30
    'WaveWatch_timeseries_10753_Wor_GBM.nc';...29
	'WaveWatch_timeseries_10462_Wor_Ob.nc';...33
	'WaveWatch_timeseries_10675_Wor_Yenisei.nc';...34
    'WaveWatch_timeseries_12247_Wor_Lena.nc';...27
    'WaveWatch_timeseries_11794_Wor_Yangtze.nc';...23
    'WaveWatch_timeseries_11499_Wor_Mahakam.nc';...19
    'WaveWatch_timeseries_12886_Wor_Fly.nc'}; % 20
	





fileNameTides_NorthAmerica = {'Wor_Mississippi_tides.csv';...15
    'Wor_Grijalva_tides.csv';...26
    'Wor_Eel_tides.csv';...22
    'Wor_Columbia_tides.csv';...24
    'Wor_Copper_tides.csv';...17
	'Wor_MacKenzie_tides.csv'}; %35

fileNameTides_SouthAmerica = {'Wor_Amazon_tides.csv';...21
    'Wor_Orinoco_tides.csv';...32
    'Wor_SaoFrancisco_tides.csv';...16
    'Wor_Parana_tides.csv';...31
    'Car_MagdalenaAtBocasDeCeniza_tides.csv';...4
    'Car_SinuAtTinajones_tides.csv';... % 7
    'Car_Atrato_tides.csv';...1
    'Pac_SanJuan_tides.csv';...
    'Pac_Patia_tides.csv';...12
    'Pac_Mira_tides.csv'};

fileNameTides_AfricaEurope = {'Wor_Congo_tides.csv';...36
	'Wor_Niger_tides.csv';...37
    'Wor_Nile_tides.csv';...
	'Wor_Zambezi_tides.csv';...38
	'Wor_Limpopo_tides.csv';...
    'Wor_Elbe_tides.csv';...25
    'Wor_Danube_tides.csv'}; %18

fileNameTides_AsiaOceania = {'Wor_Irrawaddy_tides.csv';...30
    'Wor_GBM_tides.csv';...29
	'Wor_Ob_tides.csv';...33
	'Wor_Yenisei_tides.csv';...34
    'Wor_Lena_tides.csv';...27
    'Wor_Yangtze_tides.csv';...23
    'Wor_Mahakam_tides.csv';...19
    'Wor_Fly_tides.csv'}; % 20





% Profile data (in the order of 'kmlFile' struct')
% Difference in elevation (in m)
diff_elev_NorthAmerica = [4;... Mississippi 15
    2;... Grijalva 26
    10;... Eel 22
    3;... Columbia 24
    142;... Copper 17
	13]; % MacKenzie 35
    
diff_elev_SouthAmerica =[4;... Amazon 21
    12;... Orinoco 32
    8;... Sao Francisco 16
    4;... Parana 31
    2; ... Magdalena at Bocas de Ceniza 4
    2; ... Sinu at Tinajones 7    
    1; ... Atrato 1
    4; ... San Juan at San Juan 13
    41; ... Patia 12
    9]; %  Mira 11    

diff_elev_AfricaEurope = [5;... Congo 36
	15;... Niger 37
    3;... Nile 28
	6;... Zambezi 38
	2;... Limpopo
    5;... Elbe 25
    3]; % Danube 18
    
diff_elev_AsiaOceania = [5;... Irrawaddy 30
    8;... Ganges-Brahmaputra-Meghna 29
	4;... Ob 33
	1;... Yenisei 34
    11;... Lena 27
    1;... Yangtze 23
    2;... Mahakam 19
    2]; % Fly 20







% Distance from upstream location to mouth (in km)
diff_dist_NorthAmerica = [(382-40);... Mississippi 15
    41.6;...Grijalva 26
    21.1;... Eel 22
    160;... Columbia 24
    157;... Copper 17
	389]; % MacKenzie 35
    
diff_dist_SouthAmerica = [1081;... Amazon 21
    375;... Orinoco 32
    146;... Sao Francisco 16
    296;... Parana 31
    59.5; ... Magdalena at Bocas de Ceniza 4
    (8.91-0);... Sinu at Tinajones 7
    48;... Atrato 1
    (58.6-3.0);... San Juan at San Juan 13
    (89.5-0); ... Patia 12
    (49.1-0)]; % Mira 11    
    
diff_dist_AfricaEurope = [113;... Congo 36
	221;... Niger 37
    211;... Nile 28
	90.9;... Zambezi 38
	40.4;... Limpopo
    159;... Elbe 25
    198]; % Danube 18
    
diff_dist_AsiaOceania = [245;... Irrawaddy 30
    318;... Ganges-Brahmaputra-Meghna 29
	332;... Ob 33
	307;... Yenisei 34
    359;... Lena 27
    470;... Yangtze 23    
    210;... Mahakam 19
    106]; % Fly 20


    
	







S_NorthAmerica = diff_elev_NorthAmerica./(diff_dist_NorthAmerica*1000);
S_SouthAmerica = diff_elev_SouthAmerica./(diff_dist_SouthAmerica*1000);
S_AfricaEurope = diff_elev_AfricaEurope./(diff_dist_AfricaEurope*1000);
S_AsiaOceania = diff_elev_AsiaOceania./(diff_dist_AsiaOceania*1000);


% North America
N_NorthAmerica = length(Deltanamefile_NorthAmerica);
Ndist_obs_NorthAmerica = nan(N_NorthAmerica,1);
Delta_thetaobs_NorthAmerica = nan(N_NorthAmerica,1);
theta_rightobs_NorthAmerica = nan(N_NorthAmerica,1);
theta_leftobs_NorthAmerica = nan(N_NorthAmerica,1);
R_obs_NorthAmerica = nan(N_NorthAmerica,1);
R_NorthAmerica = nan(N_NorthAmerica,1);
Ndist_pred_NorthAmerica = nan(N_NorthAmerica,1);
theta_pred_NorthAmerica = nan(N_NorthAmerica,1);
theta_maxQw_NorthAmerica = nan(N_NorthAmerica,1);
theta_minQw_NorthAmerica = nan(N_NorthAmerica,1);
Qw_obs_NorthAmerica = nan(N_NorthAmerica,1);
T_NorthAmerica = nan(N_NorthAmerica,1);
wu_NorthAmerica = nan(N_NorthAmerica,1);
wm_NorthAmerica = nan(N_NorthAmerica,1);
% % % QAngTide_NorthAmerica = nan(N_NorthAmerica,1);
wm_wu_pred_NorthAmerica = nan(N_NorthAmerica,1);
T_obs_NorthAmerica = nan(N_NorthAmerica,1);
QTide_obs_NorthAmerica = nan(N_NorthAmerica,1);
QRiver_ret_NorthAmerica = nan(N_NorthAmerica,1);
for jj = 1:N_NorthAmerica
    % r_Wave and r_River from morphology
    [Ndist_obs_NorthAmerica(jj,1), Delta_thetaobs_NorthAmerica(jj,1), theta_rightobs_NorthAmerica(jj,1), theta_leftobs_NorthAmerica(jj,1), ...
        R_obs_NorthAmerica(jj,1), Qw_obs_NorthAmerica(jj,1), ~, ~, ...
        R_NorthAmerica(jj,1), Ndist_pred_NorthAmerica(jj,1), theta_pred_NorthAmerica(jj,1), theta_maxQw_NorthAmerica(jj,1), theta_minQw_NorthAmerica(jj,1),...
        T_NorthAmerica(jj,1), wu_NorthAmerica(jj,1), wm_NorthAmerica(jj,1), wm_wu_pred_NorthAmerica(jj,1), T_obs_NorthAmerica(jj,1), QTide_obs_NorthAmerica(jj,1), ...
        QRiver_ret_NorthAmerica(jj,1)] = ...
        jfpa_MorphologyPred(QRiver_NorthAmerica(jj,1), qRiver_NorthAmerica(jj,1), QWave_NorthAmerica(jj,1), QTide_NorthAmerica(jj,1), ...
        S_NorthAmerica(jj,1), fileNameGE_NorthAmerica{jj}, fileNameWW3_NorthAmerica{jj}, fileNameTides_NorthAmerica{jj}, ...
        GEdatapath, WW3datapath, Tidesdatapath);
        
% % %     pause
% % %     
% % %     set(gcf,'PaperPositionMode','auto')
% % %     % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % %     print(gcf,'-dpng','-r800',['z_' fileNameGE{jj} '_WaveRose'])
% % %     close(gcf)
    
% % %     jj
% % %     fileNameGE_WorldEx_Gall{jj}
% % %     pause

end


% South America
N_SouthAmerica = length(Deltanamefile_SouthAmerica);
Ndist_obs_SouthAmerica = nan(N_SouthAmerica,1);
Delta_thetaobs_SouthAmerica = nan(N_SouthAmerica,1);
theta_rightobs_SouthAmerica = nan(N_SouthAmerica,1);
theta_leftobs_SouthAmerica = nan(N_SouthAmerica,1);
R_obs_SouthAmerica = nan(N_SouthAmerica,1);
R_SouthAmerica = nan(N_SouthAmerica,1);
Ndist_pred_SouthAmerica = nan(N_SouthAmerica,1);
theta_pred_SouthAmerica = nan(N_SouthAmerica,1);
theta_maxQw_SouthAmerica = nan(N_SouthAmerica,1);
theta_minQw_SouthAmerica = nan(N_SouthAmerica,1);
Qw_obs_SouthAmerica = nan(N_SouthAmerica,1);
T_SouthAmerica = nan(N_SouthAmerica,1);
wu_SouthAmerica = nan(N_SouthAmerica,1);
wm_SouthAmerica = nan(N_SouthAmerica,1);
% % % QAngTide_SouthAmerica = nan(N_SouthAmerica,1);
wm_wu_pred_SouthAmerica = nan(N_SouthAmerica,1);
T_obs_SouthAmerica = nan(N_SouthAmerica,1);
QTide_obs_SouthAmerica = nan(N_SouthAmerica,1);
QRiver_ret_SouthAmerica = nan(N_SouthAmerica,1);
for jj = 1:N_SouthAmerica
    % r_Wave and r_River from morphology
    [Ndist_obs_SouthAmerica(jj,1), Delta_thetaobs_SouthAmerica(jj,1), theta_rightobs_SouthAmerica(jj,1), theta_leftobs_SouthAmerica(jj,1), ...
        R_obs_SouthAmerica(jj,1), Qw_obs_SouthAmerica(jj,1), ~, ~, ...
        R_SouthAmerica(jj,1), Ndist_pred_SouthAmerica(jj,1), theta_pred_SouthAmerica(jj,1), theta_maxQw_SouthAmerica(jj,1), theta_minQw_SouthAmerica(jj,1),...
        T_SouthAmerica(jj,1), wu_SouthAmerica(jj,1), wm_SouthAmerica(jj,1), wm_wu_pred_SouthAmerica(jj,1), T_obs_SouthAmerica(jj,1), QTide_obs_SouthAmerica(jj,1), ...
        QRiver_ret_SouthAmerica(jj,1)] = ...
        jfpa_MorphologyPred(QRiver_SouthAmerica(jj,1), qRiver_SouthAmerica(jj,1), QWave_SouthAmerica(jj,1), QTide_SouthAmerica(jj,1), ...
        S_SouthAmerica(jj,1), fileNameGE_SouthAmerica{jj}, fileNameWW3_SouthAmerica{jj}, fileNameTides_SouthAmerica{jj}, ...
        GEdatapath, WW3datapath, Tidesdatapath);
        
% % %     pause
% % %     
% % %     set(gcf,'PaperPositionMode','auto')
% % %     % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % %     print(gcf,'-dpng','-r800',['z_' fileNameGE_WorldEx_Nien{jj} '_WaveRose'])
% % %     close(gcf)
    
% % %     jj
% % %     fileNameGE_WorldEx_Nien{jj}
% % %     pause

end


% Africa and Europe
N_AfricaEurope = length(Deltanamefile_AfricaEurope);
Ndist_obs_AfricaEurope = nan(N_AfricaEurope,1);
Delta_thetaobs_AfricaEurope = nan(N_AfricaEurope,1);
theta_rightobs_AfricaEurope = nan(N_AfricaEurope,1);
theta_leftobs_AfricaEurope = nan(N_AfricaEurope,1);
R_obs_AfricaEurope = nan(N_AfricaEurope,1);
R_AfricaEurope = nan(N_AfricaEurope,1);
Ndist_pred_AfricaEurope = nan(N_AfricaEurope,1);
theta_pred_AfricaEurope = nan(N_AfricaEurope,1);
theta_maxQw_AfricaEurope = nan(N_AfricaEurope,1);
theta_minQw_AfricaEurope = nan(N_AfricaEurope,1);
Qw_obs_AfricaEurope = nan(N_AfricaEurope,1);
T_AfricaEurope = nan(N_AfricaEurope,1);
wu_AfricaEurope = nan(N_AfricaEurope,1);
wm_AfricaEurope = nan(N_AfricaEurope,1);
% % % QAngTide_AfricaEurope = nan(N_AfricaEurope,1);
wm_wu_pred_AfricaEurope = nan(N_AfricaEurope,1);
T_obs_AfricaEurope = nan(N_AfricaEurope,1);
QTide_obs_AfricaEurope = nan(N_AfricaEurope,1);
QRiver_ret_AfricaEurope = nan(N_AfricaEurope,1);
for jj = 1:N_AfricaEurope
    % r_Wave and r_River from morphology
    [Ndist_obs_AfricaEurope(jj,1), Delta_thetaobs_AfricaEurope(jj,1), theta_rightobs_AfricaEurope(jj,1), theta_leftobs_AfricaEurope(jj,1), ...
        R_obs_AfricaEurope(jj,1), Qw_obs_AfricaEurope(jj,1), ~, ~, ...
        R_AfricaEurope(jj,1), Ndist_pred_AfricaEurope(jj,1), theta_pred_AfricaEurope(jj,1), theta_maxQw_AfricaEurope(jj,1), theta_minQw_AfricaEurope(jj,1),...
        T_AfricaEurope(jj,1), wu_AfricaEurope(jj,1), wm_AfricaEurope(jj,1), wm_wu_pred_AfricaEurope(jj,1), T_obs_AfricaEurope(jj,1), QTide_obs_AfricaEurope(jj,1), ...
        QRiver_ret_AfricaEurope(jj,1)] = ...
        jfpa_MorphologyPred(QRiver_AfricaEurope(jj,1), qRiver_AfricaEurope(jj,1), QWave_AfricaEurope(jj,1), QTide_AfricaEurope(jj,1), ...
        S_AfricaEurope(jj,1), fileNameGE_AfricaEurope{jj}, fileNameWW3_AfricaEurope{jj}, fileNameTides_AfricaEurope{jj}, ...
        GEdatapath, WW3datapath, Tidesdatapath);
        
% % %     pause
% % %     
% % %     set(gcf,'PaperPositionMode','auto')
% % %     % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % %     print(gcf,'-dpng','-r800',['z_' fileNameGE{jj} '_WaveRose'])
% % %     close(gcf)
    
% % %     jj
% % %     fileNameGE_WorldEx_Best{jj}
% % %     pause

end


% Asia and Oceania
N_AsiaOceania = length(Deltanamefile_AsiaOceania);
Ndist_obs_AsiaOceania = nan(N_AsiaOceania,1);
Delta_thetaobs_AsiaOceania = nan(N_AsiaOceania,1);
theta_rightobs_AsiaOceania = nan(N_AsiaOceania,1);
theta_leftobs_AsiaOceania = nan(N_AsiaOceania,1);
R_obs_AsiaOceania = nan(N_AsiaOceania,1);
R_AsiaOceania = nan(N_AsiaOceania,1);
Ndist_pred_AsiaOceania = nan(N_AsiaOceania,1);
theta_pred_AsiaOceania = nan(N_AsiaOceania,1);
theta_maxQw_AsiaOceania = nan(N_AsiaOceania,1);
theta_minQw_AsiaOceania = nan(N_AsiaOceania,1);
Qw_obs_AsiaOceania = nan(N_AsiaOceania,1);
T_AsiaOceania = nan(N_AsiaOceania,1);
wu_AsiaOceania = nan(N_AsiaOceania,1);
wm_AsiaOceania = nan(N_AsiaOceania,1);
% % % QAngTide_AsiaOceania = nan(N_AsiaOceania,1);
wm_wu_pred_AsiaOceania = nan(N_AsiaOceania,1);
T_obs_AsiaOceania = nan(N_AsiaOceania,1);
QTide_obs_AsiaOceania = nan(N_AsiaOceania,1);
QRiver_ret_AsiaOceania = nan(N_AsiaOceania,1);
for jj = 1:N_AsiaOceania
    [Ndist_obs_AsiaOceania(jj,1), Delta_thetaobs_AsiaOceania(jj,1), theta_rightobs_AsiaOceania(jj,1), theta_leftobs_AsiaOceania(jj,1), ...
        R_obs_AsiaOceania(jj,1), Qw_obs_AsiaOceania(jj,1), ~, ~, ...
        R_AsiaOceania(jj,1), Ndist_pred_AsiaOceania(jj,1), theta_pred_AsiaOceania(jj,1), theta_maxQw_AsiaOceania(jj,1), theta_minQw_AsiaOceania(jj,1),...
        T_AsiaOceania(jj,1), wu_AsiaOceania(jj,1), wm_AsiaOceania(jj,1), wm_wu_pred_AsiaOceania(jj,1), T_obs_AsiaOceania(jj,1), QTide_obs_AsiaOceania(jj,1), ...
        QRiver_ret_AsiaOceania(jj,1)] = ...
        jfpa_MorphologyPred(QRiver_AsiaOceania(jj,1), qRiver_AsiaOceania(jj,1), QWave_AsiaOceania(jj,1), QTide_AsiaOceania(jj,1), ...
        S_AsiaOceania(jj,1), ...
        fileNameGE_AsiaOceania{jj}, fileNameWW3_AsiaOceania{jj}, fileNameTides_AsiaOceania{jj},...
        GEdatapath, WW3datapath, Tidesdatapath);
    
% % %     pause
% % %     
% % %     set(gcf,'PaperPositionMode','auto')
% % %     % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % %     print(gcf,'-dpng','-r800',['z_' fileNameGE_ColCar_Rest{jj} '_WaveRose'])
% % %     close(gcf)

% % %     jj
% % %     fileNameGE_ColCar_Rest{jj}
    
end






% Calculate r_x from Q_wave^M, Q_river, and T^M
% North America
rRiver_obs_NorthAmerica = QRiver_ret_NorthAmerica./(QRiver_ret_NorthAmerica + QWave_NorthAmerica + QTide_obs_NorthAmerica);
rWave_obs_NorthAmerica = QWave_NorthAmerica./(QRiver_ret_NorthAmerica + QWave_NorthAmerica + QTide_obs_NorthAmerica);
rTide_obs_NorthAmerica = QTide_obs_NorthAmerica./(QRiver_ret_NorthAmerica + QWave_NorthAmerica + QTide_obs_NorthAmerica);

[rRiver_obs_NorthAmerica_log, rWave_obs_NorthAmerica_log, rTide_obs_NorthAmerica_log] = ...
    DeltaLogMaker(rRiver_obs_NorthAmerica, rWave_obs_NorthAmerica, rTide_obs_NorthAmerica);


% South America
rRiver_obs_SouthAmerica = QRiver_ret_SouthAmerica./(QRiver_ret_SouthAmerica + QWave_SouthAmerica + QTide_obs_SouthAmerica);
rWave_obs_SouthAmerica = QWave_SouthAmerica./(QRiver_ret_SouthAmerica + QWave_SouthAmerica + QTide_obs_SouthAmerica);
rTide_obs_SouthAmerica = QTide_obs_SouthAmerica./(QRiver_ret_SouthAmerica + QWave_SouthAmerica + QTide_obs_SouthAmerica);

[rRiver_obs_SouthAmerica_log, rWave_obs_SouthAmerica_log, rTide_obs_SouthAmerica_log] = ...
    DeltaLogMaker(rRiver_obs_SouthAmerica, rWave_obs_SouthAmerica, rTide_obs_SouthAmerica);



% Africa and Europe
rRiver_obs_AfricaEurope = QRiver_ret_AfricaEurope./(QRiver_ret_AfricaEurope + QWave_AfricaEurope + QTide_obs_AfricaEurope);
rWave_obs_AfricaEurope = QWave_AfricaEurope./(QRiver_ret_AfricaEurope + QWave_AfricaEurope + QTide_obs_AfricaEurope);
rTide_obs_AfricaEurope = QTide_obs_AfricaEurope./(QRiver_ret_AfricaEurope + QWave_AfricaEurope + QTide_obs_AfricaEurope);

[rRiver_obs_AfricaEurope_log, rWave_obs_AfricaEurope_log, rTide_obs_AfricaEurope_log] = ...
    DeltaLogMaker(rRiver_obs_AfricaEurope, rWave_obs_AfricaEurope, rTide_obs_AfricaEurope);



% Asia and Oceania
rRiver_obs_AsiaOceania = QRiver_ret_AsiaOceania./(QRiver_ret_AsiaOceania + QWave_AsiaOceania + QTide_obs_AsiaOceania);
rWave_obs_AsiaOceania = QWave_AsiaOceania./(QRiver_ret_AsiaOceania + QWave_AsiaOceania + QTide_obs_AsiaOceania);
rTide_obs_AsiaOceania = QTide_obs_AsiaOceania./(QRiver_ret_AsiaOceania + QWave_AsiaOceania + QTide_obs_AsiaOceania);

[rRiver_obs_AsiaOceania_log, rWave_obs_AsiaOceania_log, rTide_obs_AsiaOceania_log] = ...
    DeltaLogMaker(rRiver_obs_AsiaOceania, rWave_obs_AsiaOceania, rTide_obs_AsiaOceania);




% % % figure
% % % subplot(1,3,1)
% % % histogram(r_river,round(sqrt(length(r_river))))
% % % xlabel('$r_{river}$','Fontsize',16,'Interpreter','latex')
% % % ylabel('Counts','Fontsize',16)
% % % set(gca,'YLim',[0 100],'XLim',[0 1])
% % % 
% % % subplot(1,3,2)
% % % histogram(r_wave,round(sqrt(length(r_wave))))
% % % xlabel('$r_{wave}$','Fontsize',16,'Interpreter','latex')
% % % set(gca,'YLim',[0 100],'XLim',[0 1])
% % % 
% % % subplot(1,3,3)
% % % histogram(r_tide,round(sqrt(length(r_tide))))
% % % xlabel('$r_{tide}$','Fontsize',16,'Interpreter','latex')
% % % set(gca,'YLim',[0 100],'XLim',[0 1])
% % % % % % 
% % % % % % set(gcf,'PaperPositionMode','auto')
% % % % % % print(gcf,'-dpng','-r400','rFractions_histograms')
% % % % % % close(gcf)



%% Ternary diagrams
%{
% North American deltas
% Fluxes-based vs. Morphology-based Predictions
fontsiz = 16;
fontnam = 'Arial';
markersiz1 = 80;
markersiz2 = 170;
markercol1 = [0.6 0.6 0.6];
markeredgecol1 = [0 0 0];
markercol2 = [1 0 1];
markeredgecol2 = [0 0 0];
headsize = 0.00001;
headangle = 22.5;
linwidth = 1;
lincolor = [0.6 0.6 0.6];
tc_axis = [0.99 0.99 0.99];

figure
[~, ~, ~, ~, ~, ~] = ...
    jfpa_ternplotquiver(rTide_NorthAmerica_log,...
    rRiver_NorthAmerica_log,...
    rWave_NorthAmerica_log,...
    rTide_obs_NorthAmerica_log,...
    rRiver_obs_NorthAmerica_log,...
    rWave_obs_NorthAmerica_log,...
    fontsiz, 0.10, 0.00,...
    'scatter', markersiz1, markersiz2, markercol1, markeredgecol1, markercol2, markeredgecol2, ...
    headsize, headangle, linwidth, lincolor, tc_axis, 'p'); hold on
% Legend
xcoorleg = 0.03;
ycoorleg = 0.82;
dxcoorleg = 0.05;
dycoorleg = 0.1;
scatter(xcoorleg, ycoorleg, markersiz1, 'markerfacecolor', markercol1, 'markeredgecolor', markeredgecol1)
%quiver_tri(xcoorleg, ycoorleg, dxcoorleg, -dycoorleg, headsize, headangle, linwidth, lincolor)
scatter(xcoorleg+dxcoorleg, ycoorleg-dycoorleg, markersiz2, 'markerfacecolor', markercol2, 'markeredgecolor', 'k')
text(xcoorleg+0.02, ycoorleg,'Prediction','FontSize',fontsiz-5)
text(xcoorleg+dxcoorleg+0.02, ycoorleg-dycoorleg,'Observation','FontSize',fontsiz-5)

title('World Deltas (North America) - Prediction vs. Observation')
hlabel = ternlabel('Relative Q_t_i_d_e', 'Relative Q_r_i_v_e_r', 'Relative Q_w_a_v_e',...
    0.12, 0.05);
set(hlabel,'FontName',fontnam,'FontSize',fontsiz)

% % % set(gcf,'PaperPositionMode','auto')
% % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % print(gcf,'-dpng','-r600','y_Fig06_MorphologyDeltaPrediction_NorthAmerica')
% % % close(gcf)






% South American
% Fluxes-based vs. Morphology-based Predictions
fontsiz = 16;
fontnam = 'Arial';
markersiz1 = 80;
markersiz2 = 170;
markercol1 = [0.6 0.6 0.6];
markeredgecol1 = [0 0 0];
markercol2 = [1 0 1];
markeredgecol2 = [0 0 0];
headsize = 0.00001;
headangle = 22.5;
linwidth = 1;
lincolor = [0.6 0.6 0.6];
tc_axis = [0.99 0.99 0.99];

figure
[~, ~, ~, ~, ~, ~] = ...
    jfpa_ternplotquiver(rTide_SouthAmerica_log,...
    rRiver_SouthAmerica_log,...
    rWave_SouthAmerica_log,...
    rTide_obs_SouthAmerica_log,...
    rRiver_obs_SouthAmerica_log,...
    rWave_obs_SouthAmerica_log,...
    fontsiz, 0.10, 0.00,...
    'scatter', markersiz1, markersiz2, markercol1, markeredgecol1, markercol2, markeredgecol2, ...
    headsize, headangle, linwidth, lincolor, tc_axis); hold on
% Legend
xcoorleg = 0.03;
ycoorleg = 0.82;
dxcoorleg = 0.05;
dycoorleg = 0.1;
scatter(xcoorleg, ycoorleg, markersiz1, 'markerfacecolor', markercol1, 'markeredgecolor', markeredgecol1)
%quiver_tri(xcoorleg, ycoorleg, dxcoorleg, -dycoorleg, headsize, headangle, linwidth, lincolor)
scatter(xcoorleg+dxcoorleg, ycoorleg-dycoorleg, markersiz2, 'markerfacecolor', markercol2, 'markeredgecolor', 'k')
text(xcoorleg+0.02, ycoorleg,'Prediction','FontSize',fontsiz-5)
text(xcoorleg+dxcoorleg+0.02, ycoorleg-dycoorleg,'Observation','FontSize',fontsiz-5)

title('World Deltas (South America) - Prediction vs. Observation')
hlabel = ternlabel('Relative Q_t_i_d_e', 'Relative Q_r_i_v_e_r', 'Relative Q_w_a_v_e',...
    0.12, 0.05);
set(hlabel,'FontName',fontnam,'FontSize',fontsiz)

% % % set(gcf,'PaperPositionMode','auto')
% % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % print(gcf,'-dpng','-r600','y_Fig06_MorphologyDeltaPrediction_SouthAmerica')
% % % close(gcf)





% African and European deltas
% Fluxes-based vs. Morphology-based Predictions
fontsiz = 16;
fontnam = 'Arial';
markersiz1 = 80;
markersiz2 = 170;
markercol1 = [0.6 0.6 0.6];
markeredgecol1 = [0 0 0];
markercol2 = [1 0 1];
markeredgecol2 = [0 0 0];
headsize = 0.00001;
headangle = 22.5;
linwidth = 1;
lincolor = [0.6 0.6 0.6];
tc_axis = [0.99 0.99 0.99];

figure
[~, ~, ~, ~, ~, ~] = ...
    jfpa_ternplotquiver(rTide_AfricaEurope_log,...
    rRiver_AfricaEurope_log,...
    rWave_AfricaEurope_log,...
    rTide_obs_AfricaEurope_log,...
    rRiver_obs_AfricaEurope_log,...
    rWave_obs_AfricaEurope_log,...
    fontsiz, 0.10, 0.00,...
    'scatter', markersiz1, markersiz2, markercol1, markeredgecol1, markercol2, markeredgecol2, ...
    headsize, headangle, linwidth, lincolor, tc_axis); hold on
% Legend
xcoorleg = 0.03;
ycoorleg = 0.82;
dxcoorleg = 0.05;
dycoorleg = 0.1;
scatter(xcoorleg, ycoorleg, markersiz1, 'markerfacecolor', markercol1, 'markeredgecolor', markeredgecol1)
%quiver_tri(xcoorleg, ycoorleg, dxcoorleg, -dycoorleg, headsize, headangle, linwidth, lincolor)
scatter(xcoorleg+dxcoorleg, ycoorleg-dycoorleg, markersiz2, 'markerfacecolor', markercol2, 'markeredgecolor', 'k')
text(xcoorleg+0.02, ycoorleg,'Prediction','FontSize',fontsiz-5)
text(xcoorleg+dxcoorleg+0.02, ycoorleg-dycoorleg,'Observation','FontSize',fontsiz-5)

title('World Deltas (Africa and Europe) - Prediction vs. Observation')
hlabel = ternlabel('Relative Q_t_i_d_e', 'Relative Q_r_i_v_e_r', 'Relative Q_w_a_v_e',...
    0.12, 0.05);
set(hlabel,'FontName',fontnam,'FontSize',fontsiz)

% % % set(gcf,'PaperPositionMode','auto')
% % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % print(gcf,'-dpng','-r600','y_Fig06_MorphologyDeltaPrediction_AfricaEurope')
% % % close(gcf)






% Asia and Oceania
% Fluxes-based vs. Morphology-based Predictions
markersiz1 = 80;
markersiz2 = 170;
markercol1 = [0.6 0.6 0.6];
markeredgecol1 = [0 0 0];
markercol2 = [1 0 1];
markeredgecol2 = [0 0 0];
headsize = 0.00001;
headangle = 22.5;
linwidth = 1;
lincolor = [0.6 0.6 0.6];
tc_axis = [0.99 0.99 0.99];

figure
[~, ~, ~, ~, ~, ~] = ...
    jfpa_ternplotquiver(rTide_AsiaOceania_log,...
    rRiver_AsiaOceania_log,...
    rWave_AsiaOceania_log,...
    rTide_obs_AsiaOceania_log,...
    rRiver_obs_AsiaOceania_log,...
    rWave_obs_AsiaOceania_log,...
    fontsiz, 0.10, 0.00,...
    'scatter', markersiz1, markersiz2, markercol1, markeredgecol1, markercol2, markeredgecol2, ...
    headsize, headangle, linwidth, lincolor, tc_axis); hold on
% Legend
xcoorleg = 0.03;
ycoorleg = 0.82;
dxcoorleg = 0.05;
dycoorleg = 0.1;
scatter(xcoorleg, ycoorleg, markersiz1, 'markerfacecolor', markercol1, 'markeredgecolor', markeredgecol1)
%quiver_tri(xcoorleg, ycoorleg, dxcoorleg, -dycoorleg, headsize, headangle, linwidth, lincolor)
scatter(xcoorleg+dxcoorleg, ycoorleg-dycoorleg, markersiz2, 'markerfacecolor', markercol2, 'markeredgecolor', 'k')
text(xcoorleg+0.02, ycoorleg,'Prediction','FontSize',fontsiz-5)
text(xcoorleg+dxcoorleg+0.02, ycoorleg-dycoorleg,'Observation','FontSize',fontsiz-5)

title('World Deltas (Asia and Oceania) - Prediction vs. Observation')
hlabel = ternlabel('Relative Q_t_i_d_e', 'Relative Q_r_i_v_e_r', 'Relative Q_w_a_v_e',...
    0.12, 0.05);
set(hlabel,'FontName',fontnam,'FontSize',fontsiz)
%}
% % % set(gcf,'PaperPositionMode','auto')
% % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % print(gcf,'-dpng','-r600','y_Fig06_MorphologyDeltaPrediction_AsiaOceania')
% % % close(gcf)
%{
%% OLD - Figure 3 - Prediction delta morphology from field-based data

% % % clear all, close all, clc
% % % addpath(genpath('D:\Users\Jfpa\Dropbox\01_Research\03_Projects\08_Global-Deltas_NSF_2018-2020\0_Calculations'))
% % % loadpath = 'D:\Users\Jfpa\Dropbox\01_Research\03_Projects\08_Global-Deltas_NSF_2018-2020\';
% % % % % % loadfile = 'GlobalDeltaData_2020.mat';
% % % loadfile = 'GlobalDeltaChange-master\GlobalDeltaData.mat';
% % % load([loadpath loadfile])
% % % 
% % % % Find indices within latitudinal extents
% % % indexLat_PacificCol = find(round(MouthLat*10000)/10000>=1.6417...
% % %     & round(MouthLat*1000)/1000<=7.108);
% % % indexLat_CaribbeanCol = find(round(MouthLat*1000)/1000>=7.929...
% % %     & round(MouthLat*100)/100<=11.77);
% % % % Find indices within longitudinal extents that already fit latitudes
% % % index_PacificCol0 = indexLat_PacificCol(round(MouthLon(indexLat_PacificCol)*10)/10>=281 ...
% % %     & round(MouthLon(indexLat_PacificCol)*10)/10<=282.9);
% % % index_CaribbeanCol0 = indexLat_CaribbeanCol(round(MouthLon(indexLat_CaribbeanCol)*10)/10>=282.7 ...
% % %     & round(MouthLon(indexLat_CaribbeanCol)*10)/10<=287.5);
% % % 
% % % MouthLat_Col = [MouthLat(index_PacificCol0); MouthLat(index_CaribbeanCol0)];
% % % MouthLon_Col = [MouthLon(index_PacificCol0); MouthLon(index_CaribbeanCol0)];
% % % 
% % % % Find Basin IDs from latitudinal-longitudinal indices
% % % basinID_CaribbeanCol = BasinID(index_CaribbeanCol0);
% % % basinID_PacificCol = BasinID(index_PacificCol0);
% % % 
% % % 
% % % % % % figure
% % % % % % plot(MouthLon,MouthLat,'.','Color',[0.6 0.6 0.6]),hold on
% % % % % % plot(MouthLon(index_PacificCol0)...
% % % % % %     ,MouthLat(index_PacificCol0),'*r')
% % % % % % plot(MouthLon(index_CaribbeanCol0)...
% % % % % %     ,MouthLat(index_CaribbeanCol0),'*g')
% % % % % % xlabel('Longitude','Fontsize',20)
% % % % % % ylabel('Latitude','Fontsize',20)
% % % % % % axis equal
% % % % % % set(gca,'XLim',[279 289],'YLim',[0 13],'Fontsize',16)
% % % % % % grid on
% % % % % % 
% % % % % % % % % set(gcf,'PaperPositionMode','auto')
% % % % % % % % % print(gcf,'-dpng','-r400','map_ColombianDeltas')
% % % % % % % % % close(gcf)
% % % 
% % % % World deltas
% % % % Galloway (1975, Fig. 3)
% % % index_Mississippi = find(BasinID==426769); % 100
% % % index_SaoFrancisco = find(BasinID==84161); % 101
% % % index_Copper = find(BasinID==54); % 102
% % % index_Danube = find(BasinID==183205); % 103
% % % index_Mahakam = find(BasinID==549139); % 104
% % % index_Fly = find(BasinID==712650); % 105
% % % 
% % % % Nienhuis et al. (2020, Fig. 2)
% % % index_Amazon = find(BasinID==54097); % 106
% % % index_Eel = find(BasinID==364680); % 107
% % % index_Yangtze = find(BasinID==178016); % 108
% % % index_Columbia = find(BasinID==317358); % 109
% % % index_Elbe = find(BasinID==109965); % 110
% % % index_Grijalva = find(BasinID==99565); % 111
% % % index_Lena = find(BasinID==364); % 112
% % % index_Nile = find(BasinID==20041); % 113
% % % 
% % % % Best (2019)
% % % index_GangesBrahmaputra = find(BasinID==225914); % 114
% % % index_Irrawaddy = find(BasinID==303211); % 115
% % % index_Parana = find(BasinID==121076); % 116
% % % index_Orinoco = find(BasinID==26855); % 117
% % % index_Ob = find(BasinID==194); % 118
% % % index_Yenisei = find(BasinID==281); % 119, it has two indices, chosen the second
% % %      index_Yenisei = index_Yenisei(2,1);
% % % index_MacKenzie = find(BasinID==252); % 120
% % % index_Congo = find(BasinID==108677); % 121 
% % % index_Niger = find(BasinID==90857); % 122 
% % % index_Zambezi = find(BasinID==141681); % 123
% % % index_Limpopo = find(BasinID==154892); % 124
% % % 
% % % 
% % % index_WorldEx_Gall = [index_Mississippi; index_SaoFrancisco; index_Copper; ...
% % %     index_Danube; index_Mahakam; index_Fly];
% % % 
% % % index_WorldEx_Nien = [index_Amazon; index_Eel; index_Yangtze; ...
% % %     index_Columbia; index_Elbe; index_Grijalva; index_Lena; index_Nile];
% % % 
% % % index_WorldEx_Best = [index_GangesBrahmaputra; index_Irrawaddy; ...
% % %     index_Parana; index_Orinoco; index_Ob; index_Yenisei; ...
% % %     index_MacKenzie; index_Congo; index_Niger; index_Zambezi; index_Limpopo];
% % % 
% % % index_WorldEx_total = [index_WorldEx_Gall; index_WorldEx_Nien; index_WorldEx_Best];
% % %     index_WorldEx_total_Gall = (1:length(index_WorldEx_Gall))';
% % %     index_WorldEx_total_Nien = (1:length(index_WorldEx_Nien))' + index_WorldEx_total_Gall(end);
% % %     index_WorldEx_total_Best = (1:length(index_WorldEx_Best))' + index_WorldEx_total_Nien(end);
% % % 
% % % WorldEx_Deltanamefile_Gall = {'Mississippi';...
% % %     'SaoFrancisco';...
% % %     'Copper';...
% % %     'Danube';...
% % %     'Mahakam';...
% % %     'Fly'};
% % % 
% % % WorldEx_Deltanamefile_Nien = {'Amazon';...
% % %     'Eel';...
% % %     'Yangtze';...
% % %     'Columbia';...
% % %     'Elbe';...
% % %     'Grijalva';...
% % %     'Lena';...
% % %     'Nile'};
% % % 
% % % WorldEx_Deltanamefile_Best = {'GangesBrahmaputra';...
% % %     'Irrawaddy';...
% % %     'Parana';...
% % %     'Orinoco';...
% % %     'Ob';...
% % %     'Yenisei';...
% % %     'MacKenzie';...
% % %     'Congo';...
% % %     'Niger';...
% % %     'Zambezi';...
% % %     'Limpopo'};
% % % 
% % % 
% % % 
% % % WorldEx_Deltaname_Gall = {'Mississippi';...
% % %     'Sao Francisco';...
% % %     'Copper';...
% % %     'Danube';...
% % %     'Mahakam';...
% % %     'Fly'};
% % % 
% % % WorldEx_Deltaname_Nien = {'Amazon';...
% % %     'Eel';...
% % %     'Yangtze';...
% % %     'Columbia';...
% % %     'Elbe';...
% % %     'Grijalva';...
% % %     'Lena';...
% % %     'Nile'};
% % % 
% % % WorldEx_Deltaname_Best = {'Ganges-Brahmaputra';...
% % %     'Irrawaddy';...
% % %     'Parana';...
% % %     'Orinoco';...
% % %     'Ob';...
% % %     'Yenisei';...
% % %     'MacKenzie';...
% % %     'Congo';...
% % %     'Niger';...
% % %     'Zambezi';...
% % %     'Limpopo'};
% % % 
% % % 
% % % 
% % % 
% % % % Fraction 'r' of the total sediment flux contributed by waves, tides, and
% % % % river
% % % % [Nienhuis et al. 2020 Nature, Eq. 4; Google Earth dataset]
% % % BasinID_Col = [basinID_PacificCol; basinID_CaribbeanCol];
% % % 
% % % % Indices within compiled dataset (Pacific)
% % % % South to North
% % % index_Pac_Baudo = find(BasinID==35855);
% % % index_Pac_Caunapi = find(BasinID==41294);
% % % index_Pac_Mira = find(BasinID==41423);
% % % index_Pac_PatiaAtMajagual = find(BasinID==40136); % At Majagual
% % % index_Pac_PatiaAtSanquianga = find(BasinID==39739); % At Sanquianga
% % % index_Pac_PatiaAtIscuandeW = find(BasinID==39605); % At Iscuande West
% % % index_Pac_PatiaAtIscuandeE = find(BasinID==39648); % At Iscuande East
% % % index_Pac_SanJuanAtSanJuan = find(BasinID==37364); % At Boca San Juan
% % % index_Pac_SanJuanAtCharambira = find(BasinID==36980); % At Boca Charambira
% % % index_Pac_SanJuanAtTogoroma = find(BasinID==36812); % At Boca Togoroma
% % % index_Pac_Tola = find(BasinID==40996);
% % % 
% % % 
% % % 
% % % 
% % % % Indices within compiled list of indices from global dataset
% % % % 'basinID_PacificCol'
% % % % % % index_PacificCol = [index_Pac_Baudo;...
% % % % % %     index_Pac_Caunapi;...
% % % % % %     index_Pac_Mira;...
% % % % % %     index_Pac_PatiaAtMajagual;...
% % % % % %     index_Pac_PatiaAtSanquianga;... % At Sanquianga
% % % % % %     index_Pac_PatiaAtIscuandeE;...
% % % % % %     index_Pac_PatiaAtIscuandeW;...
% % % % % %     index_Pac_SanJuanAtSanJuan;... % At Boca San Juan
% % % % % %     index_Pac_SanJuanAtCharambira;...
% % % % % %     index_Pac_SanJuanAtTogoroma;...
% % % % % %     index_Pac_Tola];
% % % % % % 
% % % % % % indexM_PacificCol = [index_Pac_Baudo;...
% % % % % %     index_Pac_Caunapi;...
% % % % % %     index_Pac_Mira;...
% % % % % %     index_Pac_PatiaAtSanquianga;... % At Sanquianga
% % % % % %     index_Pac_SanJuanAtSanJuan;... % At Boca San Juan
% % % % % %     index_Pac_Tola];
% % % % % % 
% % % % % % 
% % % % % % PacificCol_Deltanamefile = {'Baudo';...
% % % % % %     'Caunapi';...
% % % % % %     'Mira';...
% % % % % %     'Patia';... % Patia At Sanquianga
% % % % % %     'SanJuan';... % At Boca San Juan
% % % % % %     'Tola'};
% % % % % % 
% % % % % % PacificCol_Deltaname = {'Baudo';...
% % % % % %     'Caunapi';...
% % % % % %     'Mira';...
% % % % % %     'Patia';... % Patia At Sanquianga
% % % % % %     'San Juan';... % At Boca San Juan
% % % % % %     'Tola'};
% % % 
% % % index_PacificCol = [index_Pac_SanJuanAtSanJuan;...
% % %     index_Pac_PatiaAtSanquianga;...
% % %     index_Pac_Mira];
% % % 
% % % indexM_PacificCol = [index_Pac_SanJuanAtSanJuan;...
% % %     index_Pac_PatiaAtSanquianga;...
% % %     index_Pac_Mira];
% % % 
% % % 
% % % PacificCol_Deltanamefile = {'SanJuan';...
% % %         'Patia';...
% % %         'Mira'};
% % % 
% % % PacificCol_Deltaname = {'San Juan';...
% % %         'Patia';... 
% % %         'Mira'};
% % % 
% % % 
% % % 
% % % % Indices within compiled dataset (Caribbean)
% % % % % % index_Car_AtratoAtTarena = find(BasinID_Col==28871); % Atrato at Boca Tarena
% % % % % % index_Car_AtratoAtElRoto = find(BasinID_Col==29299); % Atrato at Boca El Roto
% % % % % % index_Car_Buritaca = find(BasinID_Col==8520);
% % % % % % index_Car_Leon = find(BasinID_Col==29843);
% % % % % % index_Car_MagdalenaAtBocasDeCeniza = find(BasinID_Col==9447); % Magdalena at Bocas de Ceniza
% % % % % % index_Car_MagdalenaAtCorrea = find(BasinID_Col==18557); % Magdalena at Correa (Canal del Dique)
% % % % % % index_Car_MagdalenaAtMatunilla = find(BasinID_Col==17638); % Magdalena at Matunilla (Canal del Dique)
% % % % % % index_Car_SinuAtTinajones = find(BasinID_Col==21480); % Sinu at Tinajones
% % % % % % index_Car_Turbo = find(BasinID_Col==29238);
% % % % % % 
% % % % % % % Indices within compiled list of indices from global dataset
% % % % % % % 'basinID_CaribbeanCol'
% % % % % % index_CaribbeanCol = [index_Car_AtratoAtTarena;...
% % % % % %     index_Car_AtratoAtElRoto;...
% % % % % %     index_Car_Buritaca;...
% % % % % %     index_Car_Leon;...
% % % % % %     index_Car_MagdalenaAtBocasDeCeniza;...
% % % % % %     index_Car_MagdalenaAtCorrea;... % Magdalena at Correa (Canal del Dique)
% % % % % %     index_Car_MagdalenaAtMatunilla;... % Magdalena at Matunilla (Canal del Dique)
% % % % % %     index_Car_SinuAtTinajones;... % Sinu at Tinajones
% % % % % %     index_Car_Turbo];
% % % % % % 
% % % % % % indexM_CaribbeanCol = [index_Car_AtratoAtElRoto;...
% % % % % %     index_Car_Buritaca;...
% % % % % %     index_Car_Leon;...
% % % % % %     index_Car_MagdalenaAtBocasDeCeniza;...
% % % % % %     index_Car_MagdalenaAtCorrea;... % Magdalena at Correa (Canal del Dique)
% % % % % %     index_Car_MagdalenaAtMatunilla;... % Magdalena at Matunilla (Canal del Dique)
% % % % % %     index_Car_SinuAtTinajones;... % Sinu at Tinajones
% % % % % %     index_Car_Turbo];
% % % % % % 
% % % % % % CaribbeanCol_Deltanamefile = {'Atrato';...
% % % % % %     'Buritaca';...
% % % % % %     'Leon';...
% % % % % %     'MagdalenaAtBocasDeCeniza';... % Magdalena at Bocas de Ceniza
% % % % % %     'MagdalenaAtCorrea';... % Magdalena at Correa (Canal del Dique)
% % % % % %     'MagdalenaAtMatunilla';... % Magdalena at Matunilla (Canal del Dique)
% % % % % %     'SinuAtTinajones';... % Sinu at Tinajones
% % % % % %     'Turbo'}; 
% % % % % % 
% % % % % % CaribbeanCol_Deltaname = {'Atrato';...
% % % % % %     'Buritaca';...
% % % % % %     'Leon';...
% % % % % %     'Magdalena At Bocas De Ceniza';...
% % % % % %     'Magdalena At Correa';... % Magdalena at Correa (Canal del Dique)
% % % % % %     'Magdalena At Matunilla';... % Magdalena at Matunilla (Canal del Dique)
% % % % % %     'Sinu At Tinajones';... % Sinu at Tinajones
% % % % % %     'Turbo'};
% % % 
% % % 
% % % index_Car_AtratoAtTarena = find(BasinID==28871); % Atrato at Boca Tarena
% % % index_Car_AtratoAtElRoto = find(BasinID==29299); % Atrato at Boca El Roto
% % % index_Car_Buritaca = find(BasinID==8520);
% % % index_Car_Leon = find(BasinID==29843);
% % % index_Car_MagdalenaAtBocasDeCeniza = find(BasinID==9447); % Magdalena at Bocas de Ceniza
% % % index_Car_MagdalenaAtCorrea = find(BasinID==18557); % Magdalena at Correa (Canal del Dique)
% % % index_Car_MagdalenaAtMatunilla = find(BasinID==17638); % Magdalena at Matunilla (Canal del Dique)
% % % index_Car_SinuAtTinajones = find(BasinID==21480); % Sinu at Tinajones
% % % index_Car_Turbo = find(BasinID==29238);
% % % 
% % % % Indices within compiled list of indices from global dataset
% % % % 'basinID_CaribbeanCol'
% % % index_CaribbeanCol = [index_Car_MagdalenaAtBocasDeCeniza;...
% % %     index_Car_AtratoAtElRoto;...
% % %     index_Car_SinuAtTinajones];
% % % 
% % % indexM_CaribbeanCol = [index_Car_MagdalenaAtBocasDeCeniza;...
% % %     index_Car_AtratoAtElRoto;...
% % %     index_Car_SinuAtTinajones];
% % % 
% % % CaribbeanCol_Deltanamefile = {'MagdalenaAtBocasDeCeniza';...
% % %     'Atrato';...
% % %     'SinuAtTinajones'}; 
% % % 
% % % CaribbeanCol_Deltaname = {'Magdalena At Bocas De Ceniza';...
% % %     'Atrato';...
% % %     'Sinu At Tinajones'};
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % % Basin information
% % % % % % loadpathbasindata = 'D:\Users\Jfpa\Dropbox\01_Research\03_Projects\08_Global-Deltas_NSF_2018-2020\GlobalDeltaBasins_shp\';
% % % % % % loadfilebasindata = 'GlobalDeltaBasins.shp';
% % % % % % S = shaperead([loadpathbasindata loadfilebasindata]);
% % % % % % [~, msg] = fopen([loadpathbasindata loadfilebasindata]);
% % % % % % mapshow(S([index_PacificCol0; index_CaribbeanCol0]));
% % % % % % S_Col = S([index_PacificCol0; index_CaribbeanCol0]);
% % % % % % Discharge_prist_Col = Discharge_prist([index_PacificCol0; index_CaribbeanCol0]);
% % % % % % Qriver_prist_Col = QRiver_prist([index_PacificCol0; index_CaribbeanCol0]);
% % % 
% % % 
% % % 
% % % 
% % % % Pacific
% % % SumQ_ColPac = QRiver_prist(indexM_PacificCol)...
% % %     + QWave(indexM_PacificCol)...
% % %     + QTide(indexM_PacificCol);
% % % QRiver_ColPac = QRiver_prist(indexM_PacificCol);
% % % qRiver_ColPac = Discharge_prist(indexM_PacificCol);
% % % 
% % % 
% % % 
% % % 
% % % % Caribbean
% % % SumQ_ColCar = QRiver_prist(indexM_CaribbeanCol)...
% % %     + QWave(indexM_CaribbeanCol)...
% % %     + QTide(indexM_CaribbeanCol);
% % % QRiver_ColCar = QRiver_prist(indexM_CaribbeanCol);
% % % qRiver_ColCar = Discharge_prist(indexM_CaribbeanCol);
% % % 
% % % 
% % % 
% % % % World deltas (Galloway)
% % % QRiver_WorldEx_Gall = QRiver_prist(index_WorldEx_Gall);
% % % qRiver_WorldEx_Gall = Discharge_prist(index_WorldEx_Gall);
% % % QTide_WorldEx_Gall = QTide(index_WorldEx_Gall);
% % % QWave_WorldEx_Gall = QWave(index_WorldEx_Gall);
% % %     QWave_WorldEx_Gall(isnan(QWave_WorldEx_Gall)==1) = 0;
% % % 
% % % Sum_Q_WorldEx_Gall = QRiver_WorldEx_Gall...
% % %     + QWave_WorldEx_Gall...
% % %     + QTide_WorldEx_Gall;
% % % 
% % % % World deltas (Nienhuis et al.)
% % % QRiver_WorldEx_Nien = QRiver_prist(index_WorldEx_Nien);
% % % qRiver_WorldEx_Nien = Discharge_prist(index_WorldEx_Nien);
% % % QTide_WorldEx_Nien = QTide(index_WorldEx_Nien);
% % % QWave_WorldEx_Nien = QWave(index_WorldEx_Nien);
% % %     QWave_WorldEx_Nien(isnan(QWave_WorldEx_Nien)==1) = 0;
% % % 
% % % Sum_Q_WorldEx_Nien = QRiver_WorldEx_Nien...
% % %     + QWave_WorldEx_Nien...
% % %     + QTide_WorldEx_Nien;
% % % 
% % % % World deltas (Best)
% % % QRiver_WorldEx_Best = QRiver_prist(index_WorldEx_Best);
% % % qRiver_WorldEx_Best = Discharge_prist(index_WorldEx_Best);
% % % QTide_WorldEx_Best = QTide(index_WorldEx_Best);
% % % QWave_WorldEx_Best = QWave(index_WorldEx_Best);
% % %     QWave_WorldEx_Best(isnan(QWave_WorldEx_Best)==1) = 0;
% % % 
% % % Sum_Q_WorldEx_Best = QRiver_WorldEx_Best...
% % %     + QWave_WorldEx_Best...
% % %     + QTide_WorldEx_Best;
% % % 
% % % 
% % % 
% % % 
% % % % r_x values
% % % % World's examples
% % % rRiver_WorldEx_Gall = QRiver_WorldEx_Gall./Sum_Q_WorldEx_Gall;
% % % rWave_WorldEx_Gall = QWave_WorldEx_Gall./Sum_Q_WorldEx_Gall;
% % % rTide_WorldEx_Gall = QTide_WorldEx_Gall./Sum_Q_WorldEx_Gall;
% % % rRiver_WorldEx_Nien = QRiver_WorldEx_Nien./Sum_Q_WorldEx_Nien;
% % % rWave_WorldEx_Nien = QWave_WorldEx_Nien./Sum_Q_WorldEx_Nien;
% % % rTide_WorldEx_Nien = QTide_WorldEx_Nien./Sum_Q_WorldEx_Nien;
% % % rRiver_WorldEx_Best = QRiver_WorldEx_Best./Sum_Q_WorldEx_Best;
% % % rWave_WorldEx_Best = QWave_WorldEx_Best./Sum_Q_WorldEx_Best;
% % % rTide_WorldEx_Best = QTide_WorldEx_Best./Sum_Q_WorldEx_Best;
% % % rRiver_ColCar = QRiver_prist(indexM_CaribbeanCol)./SumQ_ColCar;
% % % rWave_ColCar = QWave(indexM_CaribbeanCol)./SumQ_ColCar;
% % % rTide_ColCar = QTide(indexM_CaribbeanCol)./SumQ_ColCar;
% % % rRiver_ColPac = QRiver_prist(indexM_PacificCol)./SumQ_ColPac;
% % % rWave_ColPac = QWave(indexM_PacificCol)./SumQ_ColPac;
% % % rTide_ColPac = QTide(indexM_PacificCol)./SumQ_ColPac;
% % % 
% % % [rRiver_WorldEx_Gall_log,rWave_WorldEx_Gall_log,rTide_WorldEx_Gall_log] = DeltaLogMaker(rRiver_WorldEx_Gall, rWave_WorldEx_Gall, rTide_WorldEx_Gall);
% % % [rRiver_WorldEx_Nien_log,rWave_WorldEx_Nien_log,rTide_WorldEx_Nien_log] = DeltaLogMaker(rRiver_WorldEx_Nien, rWave_WorldEx_Nien, rTide_WorldEx_Nien);
% % % [rRiver_WorldEx_Best_log,rWave_WorldEx_Best_log,rTide_WorldEx_Best_log] = DeltaLogMaker(rRiver_WorldEx_Best, rWave_WorldEx_Best, rTide_WorldEx_Best);
% % % [rRiver_ColCar_log,rWave_ColCar_log,rTide_ColCar_log] = DeltaLogMaker(rRiver_ColCar, rWave_ColCar, rTide_ColCar);
% % % [rRiver_ColPac_log,rWave_ColPac_log,rTide_ColPac_log] = DeltaLogMaker(rRiver_ColPac, rWave_ColPac, rTide_ColPac);
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % % Field-based ratios
% % % % Multiple distributaries (e.g., Patia) coalesce to a single delta, thus
% % % % the Qriver must be modified to account for this
% % % addpath(genpath('D:\Users\Jfpa\Dropbox\01_Research\03_Projects\08_Global-Deltas_NSF_2018-2020\2_Colombia\0_Exploracion\1_MorphologyPrediction'))
% % % GEdatapath = 'D:\Users\Jfpa\Dropbox\01_Research\03_Projects\08_Global-Deltas_NSF_2018-2020\2_Colombia\0_Exploracion\1_MorphologyPrediction\2_Morphology_GE\';
% % % WW3datapath = 'D:\Users\Jfpa\Dropbox\01_Research\03_Projects\08_Global-Deltas_NSF_2018-2020\2_Colombia\0_Exploracion\1_MorphologyPrediction\1_Waves_WW3\';
% % % Tidesdatapath = 'D:\Users\Jfpa\Dropbox\01_Research\03_Projects\08_Global-Deltas_NSF_2018-2020\2_Colombia\0_Exploracion\1_MorphologyPrediction\4_Tides_TPX\';
% % % 
% % % % % % dir([datapath '*.kml'])
% % % 
% % % % % % fileNameGE = {'Car_Atrato';...1
% % % % % %     'Car_Buritaca';...2
% % % % % %     'Car_Leon';...3
% % % % % %     'Car_MagdalenaAtBocasDeCeniza';...4
% % % % % %     'Car_MagdalenaAtCorrea';...5
% % % % % %     'Car_MagdalenaAtMatunilla';...6
% % % % % %     'Car_SinuAtTinajones';...7
% % % % % %     'Car_Turbo';...8
% % % % % %     'Pac_Baudo';...9
% % % % % %     'Pac_Caunapi';...10
% % % % % %     'Pac_Mira';...11
% % % % % %     'Pac_Patia';...12
% % % % % %     'Pac_SanJuan';...13
% % % % % %     'Pac_Tola';...14
% % % % % %     'Wor_Mississippi';...15
% % % % % %     'Wor_SaoFrancisco';...16
% % % % % %     'Wor_Copper';...17
% % % % % %     'Wor_Danube';...18
% % % % % %     'Wor_Mahakam';...19
% % % % % %     'Wor_Fly';...20
% % % % % %     'Wor_Amazon';...21
% % % % % %     'Wor_Eel';...22
% % % % % %     'Wor_Yangtze';...23
% % % % % %     'Wor_Columbia';...24
% % % % % %     'Wor_Elbe';...25
% % % % % %     'Wor_Grijalva';...26
% % % % % %     'Wor_Lena';...27
% % % % % %     'Wor_Nile';...28
% % % % % %     'Wor_GBM';...29
% % % % % %     'Wor_Irrawaddy';...30
% % % % % %     'Wor_Parana';...31
% % % % % %     'Wor_Orinoco';...32
% % % % % %     'Wor_Ob';...33
% % % % % %     'Wor_Yenisei';...34
% % % % % %     'Wor_MacKenzie';...35
% % % % % %     'Wor_Congo';...36
% % % % % %     'Wor_Niger';...37
% % % % % %     'Wor_Zambezi';...38
% % % % % %     'Wor_Limpopo'}; %39
% % % % % % 
% % % % % % fileNameWW3 = {'WaveWatch_timeseries_16085_Car_AtratoCO.nc';...1
% % % % % %     'WaveWatch_timeseries_16241_Car_BuritacaCO.nc';...2
% % % % % %     'WaveWatch_timeseries_16085_Car_LeonCO.nc';...3
% % % % % %     'WaveWatch_timeseries_16179_Car_MagdalenaAtBocasDeCenizaCO.nc';...4
% % % % % %     'WaveWatch_timeseries_16128_Car_MagdalenaAtCorreaCO.nc';...5
% % % % % %     'WaveWatch_timeseries_16128_Car_MagdalenaAtMatunillaCO.nc';...6
% % % % % %     'WaveWatch_timeseries_16130_Car_SinuAtTinajonesCO.nc';...7
% % % % % %     'WaveWatch_timeseries_16085_Car_TurboCO.nc';...8
% % % % % %     'WaveWatch_timeseries_16054_Pac_BaudoCO.nc';...9
% % % % % %     'WaveWatch_timeseries_15942_Pac_CaunapiCO.nc';...10
% % % % % %     'WaveWatch_timeseries_15943_Pac_MiraCO.nc';...11
% % % % % %     'WaveWatch_timeseries_15982_Pac_PatiaCO.nc';...12
% % % % % %     'WaveWatch_timeseries_16056_Pac_SanJuanCO.nc';...13
% % % % % %     'WaveWatch_timeseries_15942_Pac_TolaCA.nc';...14
% % % % % %     'WaveWatch_timeseries_213_Wor_Mississippi.nc';...15
% % % % % %     'WaveWatch_timeseries_17250_Wor_SaoFrancisco.nc';...16
% % % % % %     'WaveWatch_timeseries_14531_Wor_Copper.nc';...17
% % % % % %     'WaveWatch_timeseries_2882_Wor_Danube.nc';...18
% % % % % %     'WaveWatch_timeseries_11499_Wor_Mahakam.nc';...19
% % % % % %     'WaveWatch_timeseries_12886_Wor_Fly.nc';...20
% % % % % %     'WaveWatch_timeseries_17090_Wor_Amazon.nc';...21
% % % % % %     'WaveWatch_timeseries_9221_Wor_Eel.nc';...22
% % % % % %     'WaveWatch_timeseries_11794_Wor_Yangtze.nc';...23
% % % % % %     'WaveWatch_timeseries_9252_Wor_Columbia.nc';...24
% % % % % %     'WaveWatch_timeseries_5506_Wor_Elbe.nc';...25
% % % % % %     'WaveWatch_timeseries_15382_Wor_Grijalva.nc';...26
% % % % % %     'WaveWatch_timeseries_12247_Wor_Lena.nc';...27
% % % % % %     'WaveWatch_timeseries_2909_Wor_Nile.nc';...28
% % % % % %     'WaveWatch_timeseries_10753_Wor_GBM.nc';...29
% % % % % %     'WaveWatch_timeseries_10833_Wor_Irrawaddy.nc';...30
% % % % % %     'WaveWatch_timeseries_16925_Wor_Parana.nc';...31
% % % % % %     'WaveWatch_timeseries_16854_Wor_Orinoco.nc';...32
% % % % % % 	'WaveWatch_timeseries_10462_Wor_Ob.nc';...33
% % % % % % 	'WaveWatch_timeseries_10675_Wor_Yenisei.nc';...34
% % % % % % 	'WaveWatch_timeseries_14650_Wor_MacKenzie.nc';...35
% % % % % % 	'WaveWatch_timeseries_9668_Wor_Congo.nc';...36
% % % % % % 	'WaveWatch_timeseries_9597_Wor_Niger.nc';...37
% % % % % % 	'WaveWatch_timeseries_7662_Wor_Zambezi.nc';...38
% % % % % % 	'WaveWatch_timeseries_9972_Wor_Limpopo.nc'}; % 39
% % % % % % 
% % % % % % fileNameTides = {'Car_Atrato_tides.csv';...1
% % % % % %     'Car_Buritaca_tides.csv';...2
% % % % % %     'Car_Leon_tides.csv';...3
% % % % % %     'Car_MagdalenaAtBocasDeCeniza_tides.csv';...4
% % % % % %     'Car_MagdalenaAtCorrea_tides.csv';...5
% % % % % %     'Car_MagdalenaAtMatunilla_tides.csv';...6
% % % % % %     'Car_SinuAtTinajones_tides.csv';...7
% % % % % %     'Car_Turbo_tides.csv';...8
% % % % % %     'Pac_Baudo_tides.csv';...9
% % % % % %     'Pac_Caunapi_tides.csv';...10
% % % % % %     'Pac_Mira_tides.csv';...11
% % % % % %     'Pac_Patia_tides.csv';...12
% % % % % %     'Pac_SanJuan_tides.csv';...13
% % % % % %     'Pac_Tola_tides.csv';...14
% % % % % %     'Wor_Mississippi_tides.csv';...15
% % % % % %     'Wor_SaoFrancisco_tides.csv';...16
% % % % % %     'Wor_Copper_tides.csv';...17
% % % % % %     'Wor_Danube_tides.csv';...18
% % % % % %     'Wor_Mahakam_tides.csv';...19
% % % % % %     'Wor_Fly_tides.csv';...20
% % % % % %     'Wor_Amazon_tides.csv';...21
% % % % % %     'Wor_Eel_tides.csv';...22
% % % % % %     'Wor_Yangtze_tides.csv';...23
% % % % % %     'Wor_Columbia_tides.csv';...24
% % % % % %     'Wor_Elbe_tides.csv';...25
% % % % % %     'Wor_Grijalva_tides.csv';...26
% % % % % %     'Wor_Lena_tides.csv';...27
% % % % % %     'Wor_Nile_tides.csv';...28
% % % % % %     'Wor_GBM_tides.csv';...29
% % % % % %     'Wor_Irrawaddy_tides.csv';...30
% % % % % %     'Wor_Parana_tides.csv';...31
% % % % % %     'Wor_Orinoco_tides.csv';...32
% % % % % % 	'Wor_Ob_tides.csv';...33
% % % % % % 	'Wor_Yenisei_tides.csv';...34
% % % % % % 	'Wor_MacKenzie_tides.csv';...35
% % % % % % 	'Wor_Congo_tides.csv';...36
% % % % % % 	'Wor_Niger_tides.csv';...37
% % % % % % 	'Wor_Zambezi_tides.csv';...38
% % % % % % 	'Wor_Limpopo_tides.csv'}; % 39
% % % % % % 
% % % % % % 
% % % % % % % Profile data (in the order of 'kmlFile' struct')
% % % % % % % Difference in elevation (in m)
% % % % % % diff_elev = [1; ... Atrato 1
% % % % % %     7; ... Buritaca 2
% % % % % %     7; ... Leon 3
% % % % % %     2; ... Magdalena at Bocas de Ceniza 4
% % % % % %     3; ... Magdalena at Correa 5
% % % % % %     4; ... Magdalena at Matunilla 6
% % % % % %     2; ... Sinu at Tinajones 7
% % % % % %     7; ... Turbo 8
% % % % % %     6; ... Baudo 9
% % % % % %     2; ... Caunapi 10
% % % % % %     9; ... Mira 11
% % % % % %     41; ... Patia 12
% % % % % %     4; ... San Juan at San Juan 13
% % % % % %     5; ... Tola 14
% % % % % %     4;... Mississippi 15
% % % % % %     8;... Sao Francisco 16
% % % % % %     142;... Copper 17
% % % % % %     3;... Danube 18
% % % % % %     2;... Mahakam 19
% % % % % %     2;... Fly 20
% % % % % %     4;... Amazon 21
% % % % % %     10;... Eel 22
% % % % % %     1;... Yangtze 23
% % % % % %     3;... Columbia 24
% % % % % %     5;... Elbe 25
% % % % % %     2;... Grijalva 26
% % % % % %     11;... Lena 27
% % % % % %     3;... Nile 28
% % % % % %     8;... Ganges-Brahmaputra-Meghna 29
% % % % % %     5;... Irrawaddy 30
% % % % % %     4;... Parana 31
% % % % % %     12;... Orinoco 32
% % % % % % 	4;... Ob 33
% % % % % % 	1;... Yenisei 34
% % % % % % 	13;... MacKenzie 35
% % % % % % 	5;... Congo 36
% % % % % % 	15;... Niger 37
% % % % % % 	6;... Zambezi 38
% % % % % % 	2]; % Limpopo 39
% % % % % % 
% % % % % % % Distance from upstream location to mouth (in km)
% % % % % % diff_dist = [48;... Atrato 1
% % % % % %     (2.06-0); ... Buritaca 2
% % % % % %     (9.64-1.2); ... Leon 3
% % % % % %     59.5; ... Magdalena at Bocas de Ceniza 4
% % % % % %     (5.29-0.1); ... Magdalena at Correa 5
% % % % % %     (7.76-1.5); ... Magdalena at Matunilla 6
% % % % % %     (8.91-0); ... Sinu at Tinajones 7
% % % % % %     (4.03-0.25); ... Turbo 8
% % % % % %     (20-0); ... Baudo 9
% % % % % %     (15.6-0); ... Caunapi 0
% % % % % %     (49.1-0); ... Mira 11
% % % % % %     (89.5-0); ... Patia 12
% % % % % %     (58.6-3.0); ... San Juan at San Juan 13
% % % % % %     (22.8-1.5); ... Tola 14
% % % % % %     (382-40);... Mississippi 15
% % % % % %     146;... Sao Francisco 16
% % % % % %     157;... Copper 17
% % % % % %     198;... Danube 18
% % % % % %     210;... Mahakam 19
% % % % % %     106;... Fly 20
% % % % % %     1081;... Amazon 21
% % % % % %     21.1;... Eel 22
% % % % % %     470;... Yangtze 23
% % % % % %     160;... Columbia 24
% % % % % %     159;... Elbe 25
% % % % % %     41.6;...Grijalva 26
% % % % % %     359;... Lena 27
% % % % % %     211;... Nile 28
% % % % % %     318;... Ganges-Brahmaputra-Meghna 29
% % % % % %     245;... Irrawaddy 30
% % % % % %     296;... Parana 31
% % % % % %     375;... Orinoco 32
% % % % % % 	332;... Ob 33
% % % % % % 	307;... Yenisei 34
% % % % % % 	389;... MacKenzie 35
% % % % % % 	113;... Congo 36
% % % % % % 	221;... Niger 37
% % % % % % 	90.9;... Zambezi 38
% % % % % % 	40.4]; % Limpopo 39
% % % 
% % % 
% % % fileNameGE_WorldEx_Gall = {'Wor_Mississippi';...1
% % %     'Wor_SaoFrancisco';...2
% % %     'Wor_Copper';...3
% % %     'Wor_Danube';...4
% % %     'Wor_Mahakam';...5
% % %     'Wor_Fly'};% 6
% % % 
% % % 
% % % fileNameGE_WorldEx_Nien = {'Wor_Amazon';...7
% % %     'Wor_Eel';...8
% % %     'Wor_Yangtze';...9
% % %     'Wor_Columbia';...10
% % %     'Wor_Elbe';...11
% % %     'Wor_Grijalva';...12
% % %     'Wor_Lena';...13
% % %     'Wor_Nile'}; % 14
% % % 
% % % 
% % % fileNameGE_WorldEx_Best = {'Wor_GBM';...15
% % %     'Wor_Irrawaddy';...16
% % %     'Wor_Parana';...17
% % %     'Wor_Orinoco';...18
% % %     'Wor_Ob';...19
% % %     'Wor_Yenisei';...20
% % %     'Wor_MacKenzie';...21
% % %     'Wor_Congo';...22
% % %     'Wor_Niger';...23
% % %     'Wor_Zambezi';...24
% % %     'Wor_Limpopo'}; % 25
% % % 
% % % fileNameGE_ColCar_Rest = {'Car_MagdalenaAtBocasDeCeniza';...26
% % %     'Car_Atrato';...27
% % %     'Car_SinuAtTinajones'}; % 28
% % % 
% % % fileNameGE_ColPac_Rest = {'Pac_SanJuan';...29
% % %     'Pac_Patia';...30
% % %     'Pac_Mira'}; % 31
% % % 
% % % 
% % % 
% % % 
% % % fileNameWW3_WorldEx_Gall = {'WaveWatch_timeseries_213_Wor_Mississippi.nc';...15
% % %     'WaveWatch_timeseries_17250_Wor_SaoFrancisco.nc';...16
% % %     'WaveWatch_timeseries_14531_Wor_Copper.nc';...17
% % %     'WaveWatch_timeseries_2882_Wor_Danube.nc';...18
% % %     'WaveWatch_timeseries_11499_Wor_Mahakam.nc';...19
% % %     'WaveWatch_timeseries_12886_Wor_Fly.nc'}; % 20
% % % 
% % % fileNameWW3_WorldEx_Nien = {'WaveWatch_timeseries_17090_Wor_Amazon.nc';...21
% % %     'WaveWatch_timeseries_9221_Wor_Eel.nc';...22
% % %     'WaveWatch_timeseries_11794_Wor_Yangtze.nc';...23
% % %     'WaveWatch_timeseries_9252_Wor_Columbia.nc';...24
% % %     'WaveWatch_timeseries_5506_Wor_Elbe.nc';...25
% % %     'WaveWatch_timeseries_15382_Wor_Grijalva.nc';...26
% % %     'WaveWatch_timeseries_12247_Wor_Lena.nc';...27
% % %     'WaveWatch_timeseries_2909_Wor_Nile.nc'}; % 28
% % % 
% % % fileNameWW3_WorldEx_Best = {'WaveWatch_timeseries_10753_Wor_GBM.nc';...29
% % %     'WaveWatch_timeseries_10833_Wor_Irrawaddy.nc';...30
% % %     'WaveWatch_timeseries_16925_Wor_Parana.nc';...31
% % %     'WaveWatch_timeseries_16854_Wor_Orinoco.nc';...32
% % % 	'WaveWatch_timeseries_10462_Wor_Ob.nc';...33
% % % 	'WaveWatch_timeseries_10675_Wor_Yenisei.nc';...34
% % % 	'WaveWatch_timeseries_14650_Wor_MacKenzie.nc';...35
% % % 	'WaveWatch_timeseries_9668_Wor_Congo.nc';...36
% % % 	'WaveWatch_timeseries_9597_Wor_Niger.nc';...37
% % % 	'WaveWatch_timeseries_7662_Wor_Zambezi.nc';...38
% % % 	'WaveWatch_timeseries_9972_Wor_Limpopo.nc'}; % 39
% % % 
% % % fileNameWW3_ColCar_Rest = {'WaveWatch_timeseries_16179_Car_MagdalenaAtBocasDeCenizaCO.nc';...4
% % %     'WaveWatch_timeseries_16085_Car_AtratoCO.nc';...1
% % %     'WaveWatch_timeseries_16130_Car_SinuAtTinajonesCO.nc'}; % 7
% % % 
% % % fileNameWW3_ColPac_Rest = {'WaveWatch_timeseries_16056_Pac_SanJuanCO.nc';...
% % %     'WaveWatch_timeseries_15982_Pac_PatiaCO.nc';...12
% % %     'WaveWatch_timeseries_15943_Pac_MiraCO.nc'}; % 39
% % % 
% % % 
% % % 
% % % 
% % % fileNameTides_WorldEx_Gall = {'Wor_Mississippi_tides.csv';...15
% % %     'Wor_SaoFrancisco_tides.csv';...16
% % %     'Wor_Copper_tides.csv';...17
% % %     'Wor_Danube_tides.csv';...18
% % %     'Wor_Mahakam_tides.csv';...19
% % %     'Wor_Fly_tides.csv'}; % 20
% % % 
% % % fileNameTides_WorldEx_Nien = {'Wor_Amazon_tides.csv';...21
% % %     'Wor_Eel_tides.csv';...22
% % %     'Wor_Yangtze_tides.csv';...23
% % %     'Wor_Columbia_tides.csv';...24
% % %     'Wor_Elbe_tides.csv';...25
% % %     'Wor_Grijalva_tides.csv';...26
% % %     'Wor_Lena_tides.csv';...27
% % %     'Wor_Nile_tides.csv'}; % 28
% % % 
% % % fileNameTides_WorldEx_Best = {'Wor_GBM_tides.csv';...29
% % %     'Wor_Irrawaddy_tides.csv';...30
% % %     'Wor_Parana_tides.csv';...31
% % %     'Wor_Orinoco_tides.csv';...32
% % % 	'Wor_Ob_tides.csv';...33
% % % 	'Wor_Yenisei_tides.csv';...34
% % % 	'Wor_MacKenzie_tides.csv';...35
% % % 	'Wor_Congo_tides.csv';...36
% % % 	'Wor_Niger_tides.csv';...37
% % % 	'Wor_Zambezi_tides.csv';...38
% % % 	'Wor_Limpopo_tides.csv'}; % 39
% % % 
% % % fileNameTides_ColCar_Rest = {'Car_MagdalenaAtBocasDeCeniza_tides.csv';...4
% % %     'Car_Atrato_tides.csv';...1
% % %     'Car_SinuAtTinajones_tides.csv'}; % 7
% % % 
% % % fileNameTides_ColPac_Rest = {'Pac_SanJuan_tides.csv'
% % %     'Pac_Patia_tides.csv';...12
% % %     'Pac_Mira_tides.csv'};% 11
% % % 
% % % 
% % % 
% % % 
% % % % Profile data (in the order of 'kmlFile' struct')
% % % % Difference in elevation (in m)
% % % diff_elev_WorldEx_Gall = [4;... Mississippi 15
% % %     8;... Sao Francisco 16
% % %     142;... Copper 17
% % %     3;... Danube 18
% % %     2;... Mahakam 19
% % %     2]; % Fly 20
% % % 
% % % diff_elev_WorldEx_Nien =[4;... Amazon 21
% % %     10;... Eel 22
% % %     1;... Yangtze 23
% % %     3;... Columbia 24
% % %     5;... Elbe 25
% % %     2;... Grijalva 26
% % %     11;... Lena 27
% % %     3]; % Nile 28
% % % 
% % % diff_elev_WorldEx_Best = [8;... Ganges-Brahmaputra-Meghna 29
% % %     5;... Irrawaddy 30
% % %     4;... Parana 31
% % %     12;... Orinoco 32
% % % 	4;... Ob 33
% % % 	1;... Yenisei 34
% % % 	13;... MacKenzie 35
% % % 	5;... Congo 36
% % % 	15;... Niger 37
% % % 	6;... Zambezi 38
% % % 	2]; % Limpopo
% % % 
% % % diff_elev_ColCar_Rest = [2; ... Magdalena at Bocas de Ceniza 4
% % %     1; ... Atrato 1
% % %     2]; % Sinu at Tinajones 7
% % % 
% % % diff_elev_ColPac_Rest = [4; ... San Juan at San Juan 13
% % %     41; ... Patia 12
% % %     9]; %  Mira 11
% % % 
% % % 
% % % 
% % % 
% % % % Distance from upstream location to mouth (in km)
% % % diff_dist_WorldEx_Gall = [(382-40);... Mississippi 15
% % %     146;... Sao Francisco 16
% % %     157;... Copper 17
% % %     198;... Danube 18
% % %     210;... Mahakam 19
% % %     106]; % Fly 20
% % % 
% % % diff_dist_WorldEx_Nien = [1081;... Amazon 21
% % %     21.1;... Eel 22
% % %     470;... Yangtze 23
% % %     160;... Columbia 24
% % %     159;... Elbe 25
% % %     41.6;...Grijalva 26
% % %     359;... Lena 27
% % %     211]; % Nile 28
% % % 
% % % diff_dist_WorldEx_Best = [318;... Ganges-Brahmaputra-Meghna 29
% % %     245;... Irrawaddy 30
% % %     296;... Parana 31
% % %     375;... Orinoco 32
% % % 	332;... Ob 33
% % % 	307;... Yenisei 34
% % % 	389;... MacKenzie 35
% % % 	113;... Congo 36
% % % 	221;... Niger 37
% % % 	90.9;... Zambezi 38
% % % 	40.4]; % Limpopo
% % % 
% % % diff_dist_ColCar_Rest = [59.5; ... Magdalena at Bocas de Ceniza 4
% % %     48;... Atrato 1
% % %     (8.91-0)]; % Sinu at Tinajones 7
% % % 
% % % diff_dist_ColPac_Rest = [(58.6-3.0);... %  San Juan at San Juan 13
% % %     (89.5-0); ... Patia 12
% % %     (49.1-0)]; % Mira 11
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % S_WorldEx_Gall = diff_elev_WorldEx_Gall./(diff_dist_WorldEx_Gall*1000);
% % % S_WorldEx_Nien = diff_elev_WorldEx_Nien./(diff_dist_WorldEx_Nien*1000);
% % % S_WorldEx_Best = diff_elev_WorldEx_Best./(diff_dist_WorldEx_Best*1000);
% % % S_ColCar_Rest = diff_elev_ColCar_Rest./(diff_dist_ColCar_Rest*1000);
% % % S_ColPac_Rest = diff_elev_ColPac_Rest./(diff_dist_ColPac_Rest*1000);
% % % 
% % % 
% % % % World Examples
% % % % Galloway
% % % N_WorldEx_Gall = length(WorldEx_Deltanamefile_Gall);
% % % Delta_thetaM_WorldEx_Gall = nan(N_WorldEx_Gall,1);
% % % theta_rightM_Gall = nan(N_WorldEx_Gall,1);
% % % theta_leftM_Gall = nan(N_WorldEx_Gall,1);
% % % RM_WorldEx_Gall = nan(N_WorldEx_Gall,1);
% % % R_WorldEx_Gall = nan(N_WorldEx_Gall,1);
% % % theta_maxQw_Gall = nan(N_WorldEx_Gall,1);
% % % theta_minQw_Gall = nan(N_WorldEx_Gall,1);
% % % QWaveM_WorldEx_Gall = nan(N_WorldEx_Gall,1);
% % % QRiverM_WorldEx_Gall = nan(N_WorldEx_Gall,1);
% % % wu_WorldEx_Gall = nan(N_WorldEx_Gall,1);
% % % wm_WorldEx_Gall = nan(N_WorldEx_Gall,1);
% % % QAngTide_WorldEx_Gall = nan(N_WorldEx_Gall,1);
% % % TM_WorldEx_Gall = nan(N_WorldEx_Gall,1);
% % % for jj = 1:N_WorldEx_Gall
% % %     % r_Wave and r_River from morphology
% % %     [Delta_thetaM_WorldEx_Gall(jj,1), theta_rightM_Gall(jj,1), theta_leftM_Gall(jj,1), ...
% % %         RM_WorldEx_Gall(jj,1), QWaveM_WorldEx_Gall(jj,1), ~, ~, ...
% % %         QRiverM_WorldEx_Gall(jj,1), R_WorldEx_Gall(jj,1), theta_maxQw_Gall(jj,1), theta_minQw_Gall(jj,1),...
% % %         wu_WorldEx_Gall(jj,1), wm_WorldEx_Gall(jj,1), QAngTide_WorldEx_Gall(jj,1), TM_WorldEx_Gall(jj,1)] = ...
% % %         jfpa_MorphologyPred(QRiver_WorldEx_Gall(jj,1), qRiver_WorldEx_Gall(jj,1), ...
% % %         S_WorldEx_Gall(jj,1), fileNameGE_WorldEx_Gall{jj}, fileNameWW3_WorldEx_Gall{jj}, fileNameTides_WorldEx_Gall{jj}, ...
% % %         GEdatapath, WW3datapath, Tidesdatapath);
% % %         
% % % % % %     pause
% % % % % %     
% % % % % %     set(gcf,'PaperPositionMode','auto')
% % % % % %     % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % % % %     print(gcf,'-dpng','-r800',['z_' fileNameGE{jj} '_WaveRose'])
% % % % % %     close(gcf)
% % %     
% % % % % %     jj
% % % % % %     fileNameGE_WorldEx_Gall{jj}
% % % % % %     pause
% % % 
% % % end
% % % 
% % % 
% % % % Nienhuis
% % % N_WorldEx_Nien = length(WorldEx_Deltanamefile_Nien);
% % % Delta_thetaM_WorldEx_Nien = nan(N_WorldEx_Nien,1);
% % % theta_rightM_Nien = nan(N_WorldEx_Nien,1);
% % % theta_leftM_Nien = nan(N_WorldEx_Nien,1);
% % % RM_WorldEx_Nien = nan(N_WorldEx_Nien,1);
% % % R_WorldEx_Nien = nan(N_WorldEx_Nien,1);
% % % theta_maxQw_Nien = nan(N_WorldEx_Nien,1);
% % % theta_minQw_Nien = nan(N_WorldEx_Nien,1);
% % % QWaveM_WorldEx_Nien = nan(N_WorldEx_Nien,1);
% % % QRiverM_WorldEx_Nien = nan(N_WorldEx_Nien,1);
% % % wu_WorldEx_Nien = nan(N_WorldEx_Nien,1);
% % % wm_WorldEx_Nien = nan(N_WorldEx_Nien,1);
% % % QAngTide_WorldEx_Nien = nan(N_WorldEx_Nien,1);
% % % TM_WorldEx_Nien = nan(N_WorldEx_Nien,1);
% % % for jj = 1:N_WorldEx_Nien
% % %     % r_Wave and r_River from morphology
% % %     [Delta_thetaM_WorldEx_Nien(jj,1), theta_rightM_Nien(jj,1), theta_leftM_Nien(jj,1), ...
% % %         RM_WorldEx_Nien(jj,1), QWaveM_WorldEx_Nien(jj,1), ~, ~, ...
% % %         QRiverM_WorldEx_Nien(jj,1), R_WorldEx_Nien(jj,1), theta_maxQw_Nien(jj,1), theta_minQw_Nien(jj,1),...
% % %         wu_WorldEx_Nien(jj,1), wm_WorldEx_Nien(jj,1), QAngTide_WorldEx_Nien(jj,1), TM_WorldEx_Nien(jj,1)] = ...
% % %         jfpa_MorphologyPred(QRiver_WorldEx_Nien(jj,1), qRiver_WorldEx_Nien(jj,1), ...
% % %         S_WorldEx_Nien(jj,1), fileNameGE_WorldEx_Nien{jj}, fileNameWW3_WorldEx_Nien{jj}, fileNameTides_WorldEx_Nien{jj}, ...
% % %         GEdatapath, WW3datapath, Tidesdatapath);
% % %         
% % % % % %     pause
% % % % % %     
% % % % % %     set(gcf,'PaperPositionMode','auto')
% % % % % %     % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % % % %     print(gcf,'-dpng','-r800',['z_' fileNameGE_WorldEx_Nien{jj} '_WaveRose'])
% % % % % %     close(gcf)
% % %     
% % % % % %     jj
% % % % % %     fileNameGE_WorldEx_Nien{jj}
% % % % % %     pause
% % % 
% % % end
% % % 
% % % 
% % % % Best
% % % N_WorldEx_Best = length(WorldEx_Deltanamefile_Best);
% % % Delta_thetaM_Best = nan(N_WorldEx_Best,1);
% % % theta_rightM_Best = nan(N_WorldEx_Best,1);
% % % theta_leftM_Best = nan(N_WorldEx_Best,1);
% % % RM_WorldEx_Best = nan(N_WorldEx_Best,1);
% % % R_WorldEx_Best = nan(N_WorldEx_Best,1);
% % % theta_maxQw_Best = nan(N_WorldEx_Best,1);
% % % theta_minQw_Best = nan(N_WorldEx_Best,1);
% % % QWaveM_WorldEx_Best = nan(N_WorldEx_Best,1);
% % % QRiverM_WorldEx_Best = nan(N_WorldEx_Best,1);
% % % wu_WorldEx_Best = nan(N_WorldEx_Best,1);
% % % wm_WorldEx_Best = nan(N_WorldEx_Best,1);
% % % QAngTide_WorldEx_Best = nan(N_WorldEx_Best,1);
% % % TM_WorldEx_Best = nan(N_WorldEx_Best,1);
% % % for jj = 1:N_WorldEx_Best
% % %     % r_Wave and r_River from morphology
% % %     [Delta_thetaM_Best(jj,1), theta_rightM_Best(jj,1), theta_leftM_Best(jj,1), ...
% % %         RM_WorldEx_Best(jj,1), QWaveM_WorldEx_Best(jj,1), ~, ~, ...
% % %         QRiverM_WorldEx_Best(jj,1), R_WorldEx_Best(jj,1), theta_maxQw_Best(jj,1), theta_minQw_Best(jj,1),...
% % %         wu_WorldEx_Best(jj,1), wm_WorldEx_Best(jj,1), QAngTide_WorldEx_Best(jj,1), TM_WorldEx_Best(jj,1)] = ...
% % %         jfpa_MorphologyPred(QRiver_WorldEx_Best(jj,1), qRiver_WorldEx_Best(jj,1), ...
% % %         S_WorldEx_Best(jj,1), fileNameGE_WorldEx_Best{jj}, fileNameWW3_WorldEx_Best{jj}, fileNameTides_WorldEx_Best{jj}, ...
% % %         GEdatapath, WW3datapath, Tidesdatapath);
% % %         
% % % % % %     pause
% % % % % %     
% % % % % %     set(gcf,'PaperPositionMode','auto')
% % % % % %     % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % % % %     print(gcf,'-dpng','-r800',['z_' fileNameGE{jj} '_WaveRose'])
% % % % % %     close(gcf)
% % %     
% % % % % %     jj
% % % % % %     fileNameGE_WorldEx_Best{jj}
% % % % % %     pause
% % % 
% % % end
% % % 
% % % 
% % % % Caribbean
% % % N_ColCar = length(CaribbeanCol_Deltanamefile);
% % % Delta_thetaM_ColCar = nan(N_ColCar,1);
% % % theta_rightM_ColCar = nan(N_ColCar,1);
% % % theta_leftM_ColCar = nan(N_ColCar,1);
% % % RM_ColCar = nan(N_ColCar,1);
% % % R_ColCar = nan(N_ColCar,1);
% % % theta_maxQw_ColCar = nan(N_ColCar,1);
% % % theta_minQw_ColCar = nan(N_ColCar,1);
% % % QWaveM_ColCar = nan(N_ColCar,1);
% % % QRiverM_ColCar = nan(N_ColCar,1);
% % % w_u_ColCar = nan(N_ColCar,1);
% % % w_m_ColCar = nan(N_ColCar,1);
% % % QAngTide_ColCar = nan(N_ColCar,1);
% % % TM_ColCar = nan(N_ColCar,1);
% % % for jj = 1:N_ColCar
% % %     [Delta_thetaM_ColCar(jj,1), theta_rightM_ColCar(jj,1), theta_leftM_ColCar(jj,1), ...
% % %         RM_ColCar(jj,1), QWaveM_ColCar(jj,1), ~, ~, ...
% % %         QRiverM_ColCar(jj,1), R_ColCar(jj,1), theta_maxQw_ColCar(jj,1), theta_minQw_ColCar(jj,1),...
% % %         w_u_ColCar(jj,1), w_m_ColCar(jj,1), QAngTide_ColCar(jj,1), TM_ColCar(jj,1)] = ...
% % %         jfpa_MorphologyPred(QRiver_ColCar(jj,1), qRiver_ColCar(jj,1), S_ColCar_Rest(jj,1), ...
% % %         fileNameGE_ColCar_Rest{jj}, fileNameWW3_ColCar_Rest{jj}, fileNameTides_ColCar_Rest{jj},...
% % %         GEdatapath, WW3datapath, Tidesdatapath);
% % %     
% % % % % %     pause
% % % % % %     
% % % % % %     set(gcf,'PaperPositionMode','auto')
% % % % % %     % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % % % %     print(gcf,'-dpng','-r800',['z_' fileNameGE_ColCar_Rest{jj} '_WaveRose'])
% % % % % %     close(gcf)
% % % 
% % % % % %     jj
% % % % % %     fileNameGE_ColCar_Rest{jj}
% % %     
% % % end
% % % 
% % % 
% % % % Pacific
% % % N_ColPac = length(PacificCol_Deltanamefile);
% % % Delta_thetaM_ColPac = nan(N_ColPac,1);
% % % theta_rightM_ColPac = nan(N_ColPac,1);
% % % theta_leftM_ColPac = nan(N_ColPac,1);
% % % RM_ColPac = nan(N_ColPac,1);
% % % R_ColPac = nan(N_ColPac,1);
% % % theta_maxQw_ColPac = nan(N_ColPac,1);
% % % theta_minQw_ColPac = nan(N_ColPac,1);
% % % QWaveM_ColPac = nan(N_ColPac,1);
% % % QRiverM_ColPac = nan(N_ColPac,1);
% % % w_u_ColPac = nan(N_ColPac,1);
% % % w_m_ColPac = nan(N_ColPac,1);
% % % QAngTide_ColPac = nan(N_ColPac,1);
% % % TM_ColPac = nan(N_ColPac,1);
% % % for jj = 1:N_ColPac
% % %     % r_Wave and r_River from morphology
% % %     [Delta_thetaM_ColPac(jj,1), theta_rightM_ColPac(jj,1), theta_leftM_ColPac(jj,1), ...
% % %         RM_ColPac(jj,1), QWaveM_ColPac(jj,1), ~, ~, ...
% % %         QRiverM_ColPac(jj,1), R_ColPac(jj,1), theta_maxQw_ColPac(jj,1), theta_minQw_ColPac(jj,1),...
% % %         w_u_ColPac(jj,1), w_m_ColPac(jj,1), QAngTide_ColPac(jj,1), TM_ColPac(jj,1)] = ...
% % %         jfpa_MorphologyPred(QRiver_ColPac(jj,1), qRiver_ColPac(jj,1), S_ColPac_Rest(jj,1),...
% % %         fileNameGE_ColPac_Rest{jj}, fileNameWW3_ColPac_Rest{jj}, fileNameTides_ColPac_Rest{jj},...
% % %         GEdatapath, WW3datapath, Tidesdatapath);
% % %         
% % % % % %     pause
% % % % % %     
% % % % % %     set(gcf,'PaperPositionMode','auto')
% % % % % %     % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % % % %     print(gcf,'-dpng','-r800',['z_' fileNameGE_ColPac_Rest{jj} '_WaveRose'])
% % % % % %     close(gcf)
% % %     
% % % % % %     jj
% % % % % %     fileNameGE_ColPac_Rest{jj}
% % % % % %     pause
% % %     
% % % end
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % % Get a new Qriver for multiple distributaries
% % % % % % QRiverM_ColPac = QRiver_ColPac;
% % % 
% % % % Calculate r_x from Q_wave^M, Q_river, and T^M
% % % QTideM_ColPac = TM_ColPac.*QRiverM_ColPac;
% % % rRiverM_ColPac = QRiverM_ColPac./(QRiverM_ColPac + abs(QWaveM_ColPac) + QTideM_ColPac);
% % % rWaveM_ColPac = abs(QWaveM_ColPac)./(QRiverM_ColPac + abs(QWaveM_ColPac) + QTideM_ColPac);
% % % rTideM_ColPac = QTideM_ColPac./(QRiverM_ColPac+abs(QWaveM_ColPac) + QTideM_ColPac);
% % % 
% % % [rRiverM_ColPac_log, rWaveM_ColPac_log, rTideM_ColPac_log] = ...
% % %     DeltaLogMaker(rRiverM_ColPac, rWaveM_ColPac, rTideM_ColPac);
% % % 
% % % 
% % % 
% % % % % % QRiverM_ColCar = QRiver_ColCar;
% % % 
% % % QTideM_ColCar = TM_ColCar.*QRiverM_ColCar;
% % % rRiverM_ColCar = QRiverM_ColCar./(QRiverM_ColCar + abs(QWaveM_ColCar) + QTideM_ColCar);
% % % rWaveM_ColCar = abs(QWaveM_ColCar)./(QRiverM_ColCar + abs(QWaveM_ColCar) + QTideM_ColCar);
% % % rTideM_ColCar = QTideM_ColCar./(QRiverM_ColCar + abs(QWaveM_ColCar) + QTideM_ColCar);
% % % 
% % % [rRiverM_ColCar_log, rWaveM_ColCar_log, rTideM_ColCar_log] = ...
% % %     DeltaLogMaker(rRiverM_ColCar, rWaveM_ColCar, rTideM_ColCar);
% % % 
% % % 
% % % % Galloway
% % % % % % QRiverM_WorldEx_Gall = QRiver_WorldEx_Gall;
% % % 
% % % QTideM_WorldEx_Gall = TM_WorldEx_Gall.*QRiverM_WorldEx_Gall;
% % % rRiverM_WorldEx_Gall = QRiverM_WorldEx_Gall./(QRiverM_WorldEx_Gall + abs(QWaveM_WorldEx_Gall) + QTideM_WorldEx_Gall);
% % % rWaveM_WorldEx_Gall = abs(QWaveM_WorldEx_Gall)./(QRiverM_WorldEx_Gall + abs(QWaveM_WorldEx_Gall) + QTideM_WorldEx_Gall);
% % % rTideM_WorldEx_Gall = QTideM_WorldEx_Gall./(QRiverM_WorldEx_Gall + abs(QWaveM_WorldEx_Gall) + QTideM_WorldEx_Gall);
% % % 
% % % [rRiverM_WorldEx_Gall_log, rWaveM_WorldEx_Gall_log, rTideM_WorldEx_Gall_log] = ...
% % %     DeltaLogMaker(rRiverM_WorldEx_Gall, rWaveM_WorldEx_Gall, rTideM_WorldEx_Gall);
% % % 
% % % 
% % % % Nienhuis
% % % % % % QRiverM_WorldEx_Nien = QRiver_WorldEx_Nien;
% % % 
% % % QTideM_WorldEx_Nien = TM_WorldEx_Nien.*QRiverM_WorldEx_Nien;
% % % rRiverM_WorldEx_Nien = QRiverM_WorldEx_Nien./(QRiverM_WorldEx_Nien + abs(QWaveM_WorldEx_Nien) + QTideM_WorldEx_Nien);
% % % rWaveM_WorldEx_Nien = abs(QWaveM_WorldEx_Nien)./(QRiverM_WorldEx_Nien + abs(QWaveM_WorldEx_Nien) + QTideM_WorldEx_Nien);
% % % rTideM_WorldEx_Nien = QTideM_WorldEx_Nien./(QRiverM_WorldEx_Nien + abs(QWaveM_WorldEx_Nien) + QTideM_WorldEx_Nien);
% % % 
% % % [rRiverM_WorldEx_Nien_log, rWaveM_WorldEx_Nien_log, rTideM_WorldEx_Nien_log] = ...
% % %     DeltaLogMaker(rRiverM_WorldEx_Nien, rWaveM_WorldEx_Nien, rTideM_WorldEx_Nien);
% % % 
% % % 
% % % 
% % % % Best
% % % % % % QRiverM_WorldEx_Best = QRiver_WorldEx_Best;
% % % 
% % % QTideM_WorldEx_Best = TM_WorldEx_Best.*QRiverM_WorldEx_Best;
% % % rRiverM_WorldEx_Best = QRiverM_WorldEx_Best./(QRiverM_WorldEx_Best + abs(QWaveM_WorldEx_Best) + QTideM_WorldEx_Best);
% % % rWaveM_WorldEx_Best = abs(QWaveM_WorldEx_Best)./(QRiverM_WorldEx_Best + abs(QWaveM_WorldEx_Best) + QTideM_WorldEx_Best);
% % % rTideM_WorldEx_Best = QTideM_WorldEx_Best./(QRiverM_WorldEx_Best + abs(QWaveM_WorldEx_Best) + QTideM_WorldEx_Best);
% % % 
% % % [rRiverM_WorldEx_Best_log, rWaveM_WorldEx_Best_log, rTideM_WorldEx_Best_log] = ...
% % %     DeltaLogMaker(rRiverM_WorldEx_Best, rWaveM_WorldEx_Best, rTideM_WorldEx_Best);
% % % 
% % % 
% % % 
% % % 
% % % % % % figure
% % % % % % subplot(1,3,1)
% % % % % % histogram(r_river,round(sqrt(length(r_river))))
% % % % % % xlabel('$r_{river}$','Fontsize',16,'Interpreter','latex')
% % % % % % ylabel('Counts','Fontsize',16)
% % % % % % set(gca,'YLim',[0 100],'XLim',[0 1])
% % % % % % 
% % % % % % subplot(1,3,2)
% % % % % % histogram(r_wave,round(sqrt(length(r_wave))))
% % % % % % xlabel('$r_{wave}$','Fontsize',16,'Interpreter','latex')
% % % % % % set(gca,'YLim',[0 100],'XLim',[0 1])
% % % % % % 
% % % % % % subplot(1,3,3)
% % % % % % histogram(r_tide,round(sqrt(length(r_tide))))
% % % % % % xlabel('$r_{tide}$','Fontsize',16,'Interpreter','latex')
% % % % % % set(gca,'YLim',[0 100],'XLim',[0 1])
% % % % % % % % % 
% % % % % % % % % set(gcf,'PaperPositionMode','auto')
% % % % % % % % % print(gcf,'-dpng','-r400','rFractions_histograms')
% % % % % % % % % close(gcf)
% % % 
% % % 
% % % 
% % % 
% % % 
% % % % Galloway deltas
% % % % Fluxes-based vs. Morphology-based Predictions
% % % fontsiz = 16;
% % % fontnam = 'Arial';
% % % markersiz1 = 80;
% % % markersiz2 = 170;
% % % markercol1 = [0.6 0.6 0.6];
% % % markeredgecol1 = [0 0 0];
% % % markercol2 = [0 0.8 0];
% % % markeredgecol2 = [0 0 0];
% % % headsize = 0.00001;
% % % headangle = 22.5;
% % % linwidth = 1;
% % % lincolor = [0.6 0.6 0.6];
% % % tc_axis = [0.99 0.99 0.99];
% % % 
% % % figure
% % % [~, ~, ~, ~, ~, ~] = ...
% % %     jfpa_ternplotquiver(rTide_WorldEx_Gall_log,...
% % %     rRiver_WorldEx_Gall_log,...
% % %     rWave_WorldEx_Gall_log,...
% % %     rTideM_WorldEx_Gall_log,...
% % %     rRiverM_WorldEx_Gall_log,...
% % %     rWaveM_WorldEx_Gall_log,...
% % %     fontsiz, 0.10, 0.00,...
% % %     'scatter', markersiz1, markersiz2, markercol1, markeredgecol1, markercol2, markeredgecol2, ...
% % %     headsize, headangle, linwidth, lincolor, tc_axis); hold on
% % % % Legend
% % % xcoorleg = 0.03;
% % % ycoorleg = 0.82;
% % % dxcoorleg = 0.05;
% % % dycoorleg = 0.1;
% % % scatter(xcoorleg, ycoorleg, markersiz1, 'markerfacecolor', markercol1, 'markeredgecolor', markeredgecol1)
% % % quiver_tri(xcoorleg, ycoorleg, dxcoorleg, -dycoorleg, headsize, headangle, linwidth, lincolor)
% % % scatter(xcoorleg+dxcoorleg, ycoorleg-dycoorleg, markersiz2, 'markerfacecolor', markercol2, 'markeredgecolor', 'k')
% % % text(xcoorleg+0.02, ycoorleg,'Prediction','FontSize',fontsiz-5)
% % % text(xcoorleg+dxcoorleg+0.02, ycoorleg-dycoorleg,'Observation','FontSize',fontsiz-5)
% % % 
% % % title('World Deltas (Galloway, 1975) - Prediction vs. Observation')
% % % hlabel = ternlabel('Relative Q_t_i_d_e', 'Relative Q_r_i_v_e_r', 'Relative Q_w_a_v_e',...
% % %     0.12, 0.05);
% % % set(hlabel,'FontName',fontnam,'FontSize',fontsiz)
% % % 
% % % % % % set(gcf,'PaperPositionMode','auto')
% % % % % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % % % % print(gcf,'-dpng','-r600','y_Fig06_MorphologyDeltaPrediction_WorldGalloway')
% % % % % % close(gcf)
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % % Nienhuis deltas
% % % % Fluxes-based vs. Morphology-based Predictions
% % % fontsiz = 16;
% % % fontnam = 'Arial';
% % % markersiz1 = 80;
% % % markersiz2 = 170;
% % % markercol1 = [0.6 0.6 0.6];
% % % markeredgecol1 = [0 0 0];
% % % markercol2 = [0 0.8 0];
% % % markeredgecol2 = [0 0 0];
% % % headsize = 0.00001;
% % % headangle = 22.5;
% % % linwidth = 1;
% % % lincolor = [0.6 0.6 0.6];
% % % tc_axis = [0.99 0.99 0.99];
% % % 
% % % figure
% % % [~, ~, ~, ~, ~, ~] = ...
% % %     jfpa_ternplotquiver(rTide_WorldEx_Nien_log,...
% % %     rRiver_WorldEx_Nien_log,...
% % %     rWave_WorldEx_Nien_log,...
% % %     rTideM_WorldEx_Nien_log,...
% % %     rRiverM_WorldEx_Nien_log,...
% % %     rWaveM_WorldEx_Nien_log,...
% % %     fontsiz, 0.10, 0.00,...
% % %     'scatter', markersiz1, markersiz2, markercol1, markeredgecol1, markercol2, markeredgecol2, ...
% % %     headsize, headangle, linwidth, lincolor, tc_axis); hold on
% % % % Legend
% % % xcoorleg = 0.03;
% % % ycoorleg = 0.82;
% % % dxcoorleg = 0.05;
% % % dycoorleg = 0.1;
% % % scatter(xcoorleg, ycoorleg, markersiz1, 'markerfacecolor', markercol1, 'markeredgecolor', markeredgecol1)
% % % quiver_tri(xcoorleg, ycoorleg, dxcoorleg, -dycoorleg, headsize, headangle, linwidth, lincolor)
% % % scatter(xcoorleg+dxcoorleg, ycoorleg-dycoorleg, markersiz2, 'markerfacecolor', markercol2, 'markeredgecolor', 'k')
% % % text(xcoorleg+0.02, ycoorleg,'Prediction','FontSize',fontsiz-5)
% % % text(xcoorleg+dxcoorleg+0.02, ycoorleg-dycoorleg,'Observation','FontSize',fontsiz-5)
% % % 
% % % title('World Deltas (Nienhuis et al., 2020) - Prediction vs. Observation')
% % % hlabel = ternlabel('Relative Q_t_i_d_e', 'Relative Q_r_i_v_e_r', 'Relative Q_w_a_v_e',...
% % %     0.12, 0.05);
% % % set(hlabel,'FontName',fontnam,'FontSize',fontsiz)
% % % 
% % % % % % set(gcf,'PaperPositionMode','auto')
% % % % % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % % % % print(gcf,'-dpng','-r600','y_Fig06_MorphologyDeltaPrediction_WorldNienhuis')
% % % % % % close(gcf)
% % % 
% % % 
% % % 
% % % 
% % % 
% % % % Best deltas
% % % % Fluxes-based vs. Morphology-based Predictions
% % % fontsiz = 16;
% % % fontnam = 'Arial';
% % % markersiz1 = 80;
% % % markersiz2 = 170;
% % % markercol1 = [0.6 0.6 0.6];
% % % markeredgecol1 = [0 0 0];
% % % markercol2 = [0 0.8 0];
% % % markeredgecol2 = [0 0 0];
% % % headsize = 0.00001;
% % % headangle = 22.5;
% % % linwidth = 1;
% % % lincolor = [0.6 0.6 0.6];
% % % tc_axis = [0.99 0.99 0.99];
% % % 
% % % figure
% % % [~, ~, ~, ~, ~, ~] = ...
% % %     jfpa_ternplotquiver(rTide_WorldEx_Best_log,...
% % %     rRiver_WorldEx_Best_log,...
% % %     rWave_WorldEx_Best_log,...
% % %     rTideM_WorldEx_Best_log,...
% % %     rRiverM_WorldEx_Best_log,...
% % %     rWaveM_WorldEx_Best_log,...
% % %     fontsiz, 0.10, 0.00,...
% % %     'scatter', markersiz1, markersiz2, markercol1, markeredgecol1, markercol2, markeredgecol2, ...
% % %     headsize, headangle, linwidth, lincolor, tc_axis); hold on
% % % % Legend
% % % xcoorleg = 0.03;
% % % ycoorleg = 0.82;
% % % dxcoorleg = 0.05;
% % % dycoorleg = 0.1;
% % % scatter(xcoorleg, ycoorleg, markersiz1, 'markerfacecolor', markercol1, 'markeredgecolor', markeredgecol1)
% % % quiver_tri(xcoorleg, ycoorleg, dxcoorleg, -dycoorleg, headsize, headangle, linwidth, lincolor)
% % % scatter(xcoorleg+dxcoorleg, ycoorleg-dycoorleg, markersiz2, 'markerfacecolor', markercol2, 'markeredgecolor', 'k')
% % % text(xcoorleg+0.02, ycoorleg,'Prediction','FontSize',fontsiz-5)
% % % text(xcoorleg+dxcoorleg+0.02, ycoorleg-dycoorleg,'Observation','FontSize',fontsiz-5)
% % % 
% % % title('World Deltas (Best, 2019) - Prediction vs. Observation')
% % % hlabel = ternlabel('Relative Q_t_i_d_e', 'Relative Q_r_i_v_e_r', 'Relative Q_w_a_v_e',...
% % %     0.12, 0.05);
% % % set(hlabel,'FontName',fontnam,'FontSize',fontsiz)
% % % 
% % % % % % set(gcf,'PaperPositionMode','auto')
% % % % % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % % % % print(gcf,'-dpng','-r600','y_Fig06_MorphologyDeltaPrediction_WorldBest')
% % % % % % close(gcf)
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % % All Colombian Deltas
% % % % Fluxes-based vs. Morphology-based Predictions
% % % markersiz1 = 80;
% % % markersiz2 = 170;
% % % markercol1 = [0.6 0.6 0.6];
% % % markeredgecol1 = [0 0 0];
% % % markercol2 = [0 0.8 0];
% % % markeredgecol2 = [0 0 0];
% % % headsize = 0.00001;
% % % headangle = 22.5;
% % % linwidth = 1;
% % % lincolor = [0.6 0.6 0.6];
% % % tc_axis = [0.99 0.99 0.99];
% % % 
% % % figure
% % % [~, ~, ~, ~, ~, ~] = ...
% % %     jfpa_ternplotquiver([rTide_ColCar_log;rTide_ColPac_log],...
% % %     [rRiver_ColCar_log;rRiver_ColPac_log],...
% % %     [rWave_ColCar_log;rWave_ColPac_log],...
% % %     [rTideM_ColCar_log;rTideM_ColPac_log],...
% % %     [rRiverM_ColCar_log;rRiverM_ColPac_log],...
% % %     [rWaveM_ColCar_log;rWaveM_ColPac_log],...
% % %     fontsiz, 0.10, 0.00,...
% % %     'scatter', markersiz1, markersiz2, markercol1, markeredgecol1, markercol2, markeredgecol2, ...
% % %     headsize, headangle, linwidth, lincolor, tc_axis); hold on
% % % % Legend
% % % xcoorleg = 0.03;
% % % ycoorleg = 0.82;
% % % dxcoorleg = 0.05;
% % % dycoorleg = 0.1;
% % % scatter(xcoorleg, ycoorleg, markersiz1, 'markerfacecolor', markercol1, 'markeredgecolor', markeredgecol1)
% % % quiver_tri(xcoorleg, ycoorleg, dxcoorleg, -dycoorleg, headsize, headangle, linwidth, lincolor)
% % % scatter(xcoorleg+dxcoorleg, ycoorleg-dycoorleg, markersiz2, 'markerfacecolor', markercol2, 'markeredgecolor', 'k')
% % % text(xcoorleg+0.02, ycoorleg,'Prediction','FontSize',fontsiz-5)
% % % text(xcoorleg+dxcoorleg+0.02, ycoorleg-dycoorleg,'Observation','FontSize',fontsiz-5)
% % % 
% % % title('Colombian Deltas (Restrepo and Lopez, 2008) - Prediction vs. Observation')
% % % hlabel = ternlabel('Relative Q_t_i_d_e', 'Relative Q_r_i_v_e_r', 'Relative Q_w_a_v_e',...
% % %     0.12, 0.05);
% % % set(hlabel,'FontName',fontnam,'FontSize',fontsiz)
% % % 
% % % % % % set(gcf,'PaperPositionMode','auto')
% % % % % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % % % % print(gcf,'-dpng','-r600','y_Fig06_MorphologyDeltaPrediction_WorldRestrepo')
% % % % % % close(gcf)
%}
%% CONFUSION MATRIX

Conf_Matrix = nan(4,4);

% Cartesian coordinates of deltas
[x_NorthAmerica, y_NorthAmerica] = ...
    terncoords(rTide_NorthAmerica, rRiver_NorthAmerica, rWave_NorthAmerica);
[x_SouthAmerica, y_SouthAmerica] = ...
    terncoords(rTide_SouthAmerica, rRiver_SouthAmerica, rWave_SouthAmerica);
[x_AfricaEurope, y_AfricaEurope] = ...
    terncoords(rTide_AfricaEurope, rRiver_AfricaEurope, rWave_AfricaEurope);
[x_AsiaOceania, y_AsiaOceania] = ...
    terncoords(rTide_AsiaOceania, rRiver_AsiaOceania, rWave_AsiaOceania);

[x_obs_NorthAmerica, y_obs_NorthAmerica] = ...
    terncoords(rTide_obs_NorthAmerica, rRiver_obs_NorthAmerica, rWave_obs_NorthAmerica);
[x_obs_SouthAmerica, y_obs_SouthAmerica] = ...
    terncoords(rTide_obs_SouthAmerica, rRiver_obs_SouthAmerica, rWave_obs_SouthAmerica);
[x_obs_AfricaEurope, y_obs_AfricaEurope] = ...
    terncoords(rTide_obs_AfricaEurope, rRiver_obs_AfricaEurope, rWave_obs_AfricaEurope);
[x_obs_AsiaOceania, y_obs_AsiaOceania] = ...
    terncoords(rTide_obs_AsiaOceania, rRiver_obs_AsiaOceania, rWave_obs_AsiaOceania);



% Polygon coordinates -from ternplot- (tide, river, wave)
% Tide Dominance
% Ternary points: (1,0,0), (0.5,0.5,0), (1/3,1/3,1/3), (0.5,0,0.5)
[x_tidedom, y_tidedom] = ...
    terncoords([1 0.5 1/3 0.5 1], [0 0.5 1/3 0 0], [0 0 1/3 0.5 0]);
% River Dominance
% Ternary points: (0,1,0), (0.5,0,0.5), (1/3,1/3,1/3), (0,0.5,0.5)
[x_riverdom, y_riverdom] = ...
    terncoords([0 0 1/3 0.5 0], [1 0.5 1/3 0.5 1], [0 0.5 1/3 0 0]);
% Wave Dominance
% Ternary points: (0,0,1), (0.5,0,0.5), (1/3,1/3,1/3), (0,0.5,0.5)
[x_wavedom, y_wavedom] = ...
    terncoords([0 0.5 1/3 0 0], [0 0 1/3 0.5 0], [1 0.5 1/3 0.5 1]);

figure
plot(x_tidedom,y_tidedom,'-r')
hold on
plot(x_riverdom,y_riverdom,'-g')
plot(x_wavedom,y_wavedom,'-b')
plot(x_AsiaOceania,y_AsiaOceania,'.k')
plot(x_NorthAmerica,y_NorthAmerica,'ok')
plot(x_SouthAmerica,y_SouthAmerica,'ok')
plot(x_AfricaEurope,y_AfricaEurope,'ok')
plot(x_obs_AsiaOceania,y_obs_AsiaOceania,'.m')
plot(x_obs_NorthAmerica,y_obs_NorthAmerica,'om')
plot(x_obs_SouthAmerica,y_obs_SouthAmerica,'*m')
plot(x_obs_AfricaEurope,y_obs_AfricaEurope,'+m')
plot(x_obs_AsiaOceania(3),y_obs_AsiaOceania(3),'dy')


% Observations count
% Tide-dominated deltas
tidedom_obs_NorthAmerica = ...
    inpolygon(x_obs_NorthAmerica, y_obs_NorthAmerica, x_tidedom, y_tidedom);
tidedom_obs_SouthAmerica = ...
    inpolygon(x_obs_SouthAmerica, y_obs_SouthAmerica, x_tidedom, y_tidedom);
tidedom_obs_AfricaEurope = ...
    inpolygon(x_obs_AfricaEurope, y_obs_AfricaEurope, x_tidedom, y_tidedom);
tidedom_obs_AsiaOceania = ...
    inpolygon(x_obs_AsiaOceania, y_obs_AsiaOceania, x_tidedom, y_tidedom);

tidedom_obs = [tidedom_obs_NorthAmerica; tidedom_obs_SouthAmerica; tidedom_obs_AfricaEurope;...
    tidedom_obs_AsiaOceania];
index_tidedom_obs = find(tidedom_obs==1);


% River-dominated deltas
riverdom_obs_NorthAmerica = ...
    inpolygon(x_obs_NorthAmerica, y_obs_NorthAmerica, x_riverdom, y_riverdom);
riverdom_obs_SouthAmerica = ...
    inpolygon(x_obs_SouthAmerica, y_obs_SouthAmerica, x_riverdom, y_riverdom);
riverdom_obs_AfricaEurope = ...
    inpolygon(x_obs_AfricaEurope, y_obs_AfricaEurope, x_riverdom, y_riverdom);
riverdom_obs_AsiaOceania = ...
    inpolygon(x_obs_AsiaOceania, y_obs_AsiaOceania, x_riverdom, y_riverdom);

riverdom_obs = [riverdom_obs_NorthAmerica; riverdom_obs_SouthAmerica; riverdom_obs_AfricaEurope;...
    riverdom_obs_AsiaOceania];
index_riverdom_obs = find(riverdom_obs==1);


% Wave-dominated deltas
wavedom_obs_NorthAmerica = ...
    inpolygon(x_obs_NorthAmerica, y_obs_NorthAmerica, x_wavedom, y_wavedom);
wavedom_obs_SouthAmerica = ...
    inpolygon(x_obs_SouthAmerica, y_obs_SouthAmerica, x_wavedom, y_wavedom);
wavedom_obs_AfricaEurope = ...
    inpolygon(x_obs_AfricaEurope, y_obs_AfricaEurope, x_wavedom, y_wavedom);
wavedom_obs_AsiaOceania = ...
    inpolygon(x_obs_AsiaOceania, y_obs_AsiaOceania, x_wavedom, y_wavedom);

wavedom_obs = [wavedom_obs_NorthAmerica; wavedom_obs_SouthAmerica; wavedom_obs_AfricaEurope;...
    wavedom_obs_AsiaOceania];
index_wavedom_obs = find(wavedom_obs==1);


dominance_obs0 = [tidedom_obs riverdom_obs wavedom_obs];
dominance_obs = dominance_obs0;

% % %     rTideM_all = [rTideM_ColCar; rTideM_ColPac; rTideM_WorldEx];
% % %     rRiverM_all = [rRiverM_ColCar; rRiverM_ColPac; rRiverM_WorldEx];
% % %     rWaveM_all = [rWaveM_ColCar; rWaveM_ColPac; rWaveM_WorldEx];
% % %     dominanceM = dominanceM.*[rTideM_all rRiverM_all rWaveM_all];

indexrow_diffone_obs = find(sum(dominance_obs,2)~=1);

for ii = 1:length(indexrow_diffone_obs)
    index_one = find(dominance_obs(indexrow_diffone_obs(ii),:)==1);
    dominance_obs(indexrow_diffone_obs(ii),:) = 0;
    dominance_obs(indexrow_diffone_obs(ii),max(index_one)) = 1;
end

% Re-assign index values of tide, river, and wave dominance
tidedom_obs = dominance_obs(:,1);
    index_tidedom_obs = find(tidedom_obs==1);
riverdom_obs = dominance_obs(:,2);
    index_riverdom_obs = find(riverdom_obs==1);
wavedom_obs = dominance_obs(:,3);
    index_wavedom_obs = find(wavedom_obs==1);



% Predictions count
% Tide-dominated deltas
tidedom_NorthAmerica = ...
    inpolygon(x_NorthAmerica, y_NorthAmerica, x_tidedom, y_tidedom);
tidedom_SouthAmerica = ...
    inpolygon(x_SouthAmerica, y_SouthAmerica, x_tidedom, y_tidedom);
tidedom_AfricaEurope = ...
    inpolygon(x_AfricaEurope, y_AfricaEurope, x_tidedom, y_tidedom);
tidedom_AsiaOceania = ...
    inpolygon(x_AsiaOceania, y_AsiaOceania, x_tidedom, y_tidedom);

tidedom = [tidedom_NorthAmerica; tidedom_SouthAmerica; tidedom_AfricaEurope;...
    tidedom_AsiaOceania];
index_tidedom = find(tidedom==1);


% River-dominated deltas
riverdom_NorthAmerica = ...
    inpolygon(x_NorthAmerica, y_NorthAmerica, x_riverdom, y_riverdom);
riverdom_SouthAmerica = ...
    inpolygon(x_SouthAmerica, y_SouthAmerica, x_riverdom, y_riverdom);
riverdom_AfricaEurope = ...
    inpolygon(x_AfricaEurope, y_AfricaEurope, x_riverdom, y_riverdom);
riverdom_AsiaOceania = ...
    inpolygon(x_AsiaOceania, y_AsiaOceania, x_riverdom, y_riverdom);

riverdom = [riverdom_NorthAmerica; riverdom_SouthAmerica; riverdom_AfricaEurope;...
    riverdom_AsiaOceania];
index_riverdom = find(riverdom==1);


% Wave-dominated deltas
wavedom_NorthAmerica = ...
    inpolygon(x_NorthAmerica, y_NorthAmerica, x_wavedom, y_wavedom);
wavedom_SouthAmerica = ...
    inpolygon(x_SouthAmerica, y_SouthAmerica, x_wavedom, y_wavedom);
wavedom_AfricaEurope = ...
    inpolygon(x_AfricaEurope, y_AfricaEurope, x_wavedom, y_wavedom);
wavedom_AsiaOceania = ...
    inpolygon(x_AsiaOceania, y_AsiaOceania, x_wavedom, y_wavedom);

wavedom = [wavedom_NorthAmerica; wavedom_SouthAmerica; wavedom_AfricaEurope;...
    wavedom_AsiaOceania];
index_wavedom = find(wavedom==1);


dominance0 = [tidedom riverdom wavedom];
dominance = dominance0;

indexrow_diffone = find(sum(dominance,2)~=1);

for ii = 1:length(indexrow_diffone)
    index_one = find(dominance(indexrow_diffone(ii),:)==1);
    dominance(indexrow_diffone(ii),:) = 0;
    dominance(indexrow_diffone(ii),max(index_one)) = 1;
end

tidedom = dominance(:,1);
riverdom = dominance(:,2);
wavedom = dominance(:,3);





% Count for Confusion Matrix
% format: observation_prediction
tide_tide = sum(dominance_obs(index_tidedom,1));
tide_river = sum(dominance_obs(index_riverdom,1));
tide_wave = sum(dominance_obs(index_wavedom,1));

river_tide = sum(dominance_obs(index_tidedom,2));
river_river = sum(dominance_obs(index_riverdom,2));
river_wave = sum(dominance_obs(index_wavedom,2));

wave_tide = sum(dominance_obs(index_tidedom,3));
wave_river = sum(dominance_obs(index_riverdom,3));
wave_wave = sum(dominance_obs(index_wavedom,3));




Conf_Matrix = zeros(4,4);

Conf_Matrix(1,1) = tide_tide;
Conf_Matrix(1,2) = tide_river;
Conf_Matrix(1,3) = tide_wave;
Conf_Matrix(1,4) = sum(Conf_Matrix(1,1:3)); % 
Conf_Matrix(4,1) = sum(Conf_Matrix(1:3,1));

Conf_Matrix(2,1) = river_tide;
Conf_Matrix(2,2) = river_river;
Conf_Matrix(2,3) = river_wave;
Conf_Matrix(2,4) = sum(Conf_Matrix(2,1:3));
Conf_Matrix(4,2) = sum(Conf_Matrix(1:3,2));

Conf_Matrix(3,1) = wave_tide;
Conf_Matrix(3,2) = wave_river;
Conf_Matrix(3,3) = wave_wave;
Conf_Matrix(3,4) = sum(Conf_Matrix(3,1:3));
Conf_Matrix(4,3) = sum(Conf_Matrix(1:3,3));

accuracy = sum(diag(Conf_Matrix))/sum(Conf_Matrix(:,4))*100;








% R^M and T^M, R and T COMPARISON
QRiver_all = [QRiver_NorthAmerica; QRiver_SouthAmerica; QRiver_AfricaEurope; QRiver_AsiaOceania];
%T_AsiaOceania
QTide_all = [T_NorthAmerica; T_SouthAmerica; T_AfricaEurope; T_AsiaOceania].*QRiver_all;
%QTide_all = [QTide_NorthAmerica; QTide_SouthAmerica; QTide_AfricaEurope; QTide_AsiaOceania];
QWave_all = [QWave_NorthAmerica; QWave_SouthAmerica; QWave_AfricaEurope; QWave_AsiaOceania];

R_all = QRiver_all./QWave_all;
T_all = QTide_all./QRiver_all;


R_obs_all = [R_obs_NorthAmerica; R_obs_SouthAmerica; R_obs_AfricaEurope; R_obs_AsiaOceania];
T_obs_all = [T_obs_NorthAmerica; T_obs_SouthAmerica; T_obs_AfricaEurope; T_obs_AsiaOceania];

figure
subplot(1,2,1)
scatter(R_obs_all, R_all, 'ok','filled')
set(gca,'XScale','log','YScale','log')
xlabel('R observed'),ylabel('R predicted')
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1])
subplot(1,2,2)
scatter(T_obs_all, T_all, 'ok','filled')
set(gca,'XScale','log','YScale','log')
xlabel('T observed'),ylabel('T predicted')
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1])
disp('blub')
%}
%% Map of delta location
% % % close all
% % % addpath(genpath('D:\Users\Jfpa\m_map'))
% % % 
% % % 
% % % loadpathbasindata = 'D:\Users\Jfpa\Dropbox\01_Research\03_Projects\08_Global-Deltas_NSF_2018-2020\GlobalDeltaBasins_shp\';
% % % loadfilebasindata = 'GlobalDeltaBasins.shp';
% % % S = shaperead([loadpathbasindata loadfilebasindata]);
% % % [~, msg] = fopen([loadpathbasindata loadfilebasindata]);
% % % mapshow(S([index_PacificCol0; index_CaribbeanCol0]));
% % % S_Col = S([index_PacificCol0; index_CaribbeanCol0]);
% % % Discharge_prist_Col = Discharge_prist([index_PacificCol; index_AsiaOceania]);
% % % Qriver_prist_Col = QRiver_prist([index_PacificCol; index_AsiaOceania]);
% % % 
% % % 
% % % figure, hold on
% % % markersiz = 10;
% % % 
% % % m_coord('geographic');
% % % m_proj('UTM','longitudes',[-80 -70], ...
% % %     'latitudes',[0 13]);%,'direction','horizontal','aspect',.5);
% % % % % % [CS,CH] = m_etopo2('contour',[-70:10:-20 -10 -5 -2 -1],'edgecolor','b');
% % % % % % m_elev('contour',-70:10:-10,'edgecolor','b');
% % % m_gshhs_h('patch',[0.8 0.8 0.8],'edgecolor','k');
% % % % % % m_gshhs('fr','linewidth',2); % Rivers
% % % % Color scale for discharge
% % % % % % Dis_scale = [S_Col.Dis_dist]./max([S_Col.Dis_dist]);
% % % parula_int = interp1(1:256, parula, 256/99:256/99:256);
% % % [Dis_dist_sort, index_sort] = ...
% % %     sort(Discharge_prist_Col);
% % % for ii = 1:length(S_Col)
% % %     m_patch(S_Col(ii).X, S_Col(ii).Y, [0.6 0.6 0.6]);
% % % end
% % % m_grid('fontsize',fontsiz);
% % % % % % m_northarrow(-77.25,8.7,0.08,'type',2);
% % % % % % m_gshhs('fr','linewidth',2,'Color',[0 0.7 1]); % Rivers
% % % m_ruler([0.45 0.7],0.1,3,'fontsize',fontsiz-6)
% % % 
% % % m_plot(MouthLon_Col(index_AsiaOceania)-360,...
% % %     MouthLat_Col(index_AsiaOceania),'^','MarkerFaceColor',[0 0.7 0],...
% % %     'MarkerEdgeColor','k','MarkerSize', markersiz)
% % % m_plot(MouthLon_Col(index_PacificCol)-360,...
% % %     MouthLat_Col(index_PacificCol),'^','MarkerFaceColor',[0 0.7 0],...
% % %     'MarkerEdgeColor','k','MarkerSize', markersiz)


% % % set(gcf,'PaperPositionMode','auto')
% % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % print(gcf,'-dpng','-r600','y_Fig06_MorphBasedPrediction_Location')
% % % close(gcf)

%% FLUVIAL FACTOR (FF) TO INCREASE ACCURACY
QRiver_all = [QRiver_NorthAmerica; QRiver_SouthAmerica; QRiver_AfricaEurope; QRiver_AsiaOceania];
QTide_all = [QTide_NorthAmerica; QTide_SouthAmerica; QTide_AfricaEurope; QTide_AsiaOceania];
QWave_all = [QWave_NorthAmerica; QWave_SouthAmerica; QWave_AfricaEurope; QWave_AsiaOceania];

y_obs = [y_obs_NorthAmerica; y_obs_SouthAmerica; y_obs_AfricaEurope; y_obs_AsiaOceania];
x_obs = [x_obs_NorthAmerica; x_obs_SouthAmerica; x_obs_AfricaEurope; x_obs_AsiaOceania];


FF = nan(length(QRiver_all),1);

% ***Just fluvial (FF)***
for i = 1:length(FF)
fun_FF = @(FF_var) sqrt((x_obs(i) - (QTide_all(i)./(QRiver_all(i).*FF_var+QWave_all(i)+QTide_all(i))...
    + (QRiver_all(i).*FF_var./(QRiver_all(i).*FF_var+QWave_all(i)+QTide_all(i))).*sind(60).*cotd(60))).^2 ...
    + (y_obs(i) - ((QRiver_all(i).*FF_var./(QRiver_all(i)*FF_var+QWave_all(i)+QTide_all(i))).*sind(60))).^2);
options = optimset('TolFun', 5e-16, 'TolX', 5e-16);
% % % SR = fminsearch(fun_SR, SR0, options);
FF(i,1) = fminbnd(fun_FF, 0, 1000);
end

rRiver1 = QRiver_all.*FF./(QRiver_all.*FF+QWave_all+QTide_all);
rWave1 = QWave_all./(QRiver_all.*FF+QWave_all+QTide_all);
rTide1 = QTide_all./(QRiver_all.*FF+QWave_all+QTide_all);

y1 = rRiver1.*sind(60);
x1 = rTide1 + rRiver1.*sind(60).*cotd(60);

Error_rx1 = sqrt((x_obs-x1).^2+ (y_obs-y1).^2);

Error_rx1_NorthAmerica = sqrt((x_obs_NorthAmerica-x1(index_WorldEx_total_NA)).^2 ...
    + (y_obs_NorthAmerica-y1(index_WorldEx_total_NA)).^2);
Error_rx1_SouthAmerica = sqrt((x_obs_SouthAmerica-x1(index_WorldEx_total_SA)).^2 ...
    + (y_obs_SouthAmerica-y1(index_WorldEx_total_SA)).^2);
Error_rx1_AfricaEurope = sqrt((x_obs_AfricaEurope-x1(index_WorldEx_total_AE)).^2 ...
    + (y_obs_AfricaEurope-y1(index_WorldEx_total_AE)).^2);
Error_rx1_AsiaOceania = sqrt((x_obs_AsiaOceania-x1(index_WorldEx_total_AO)).^2 ...
    + (y_obs_AsiaOceania-y1(index_WorldEx_total_AO)).^2);




% Fluvial factor (FF)
figure
plot(0:32,(0:32).*0+1, '-k', 'Linewidth', 3), hold on
plot(index_riverdom_obs, FF(index_riverdom_obs), ...
    's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g',...
    'MarkerSize',15), hold on
plot(index_tidedom_obs, FF(index_tidedom_obs), ...
    'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r',...
    'MarkerSize',12)
plot(index_wavedom_obs, FF(index_wavedom_obs), ...
    'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b',...
    'MarkerSize',12)
set(gca,'XLim', [0 32], 'YLim', [1e-6 1e3])
set(gca,'XTick', 0:1:31, 'YTick', [1e-6 1e-3 1e0 1e3])
set(gca,'XTickLabel',{'' '1' '' '' '' '' '6' '' '' '' '' '11' '' '' '' '' '16' '' '' '' '' '21' '' '' '' '' '26' '' '' '' '' '31'})
set(gca,'FontSize',25)
set(gca,'YScale','log')
xtickangle(0)
xlabel('Delta','FontSize',25,'FontName','Arial')
ylabel('Fluvial Factor, $FF$','Interpreter','latex',...
    'FontSize',30,'FontName','Arial')
grid on, box on

% % % set(gcf,'PaperPositionMode','auto')
% % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % print(gcf,'-dpng','-r600','y_Fig06_MorphBasedRatios_FluvialFactor')
% % % close(gcf)

%% FACTORS (FF, TF, AND WF) TO INCREASE ACCURACY
%{
% % % QRiver_all = [QRiver_NorthAmerica; QRiver_SouthAmerica; QRiver_AfricaEurope; QRiver_AsiaOceania];
% % % QTide_all = [QTide_NorthAmerica; QTide_SouthAmerica; QTide_AfricaEurope; QTide_AsiaOceania];
% % % QWave_all = [QWave_NorthAmerica; QWave_SouthAmerica; QWave_AfricaEurope; QWave_AsiaOceania];
% % % 
% % % y_obs_all = [y_obs_NorthAmerica; y_obs_SouthAmerica; y_obs_AfricaEurope; y_obs_AsiaOceania];
% % % x_obs_all = [x_obs_NorthAmerica; x_obs_SouthAmerica; x_obs_AfricaEurope; x_obs_AsiaOceania];
% % % % % % delta = 4;
% % % % % % QRiver_i = QRiver_all(delta);
% % % % % % QTide_i = QTide_all(delta);
% % % % % % QWave_i = QWave_all(delta);
% % % QRiver_i = QRiver_all;
% % % QTide_i = QTide_all;
% % % QWave_i = QWave_all;
% % % 
% % % FF0 = zeros(length(QRiver_i),1);
% % % rRiver = QRiver_i.*FF0./(QRiver_i.*FF0+QWave_i+QTide_i);
% % % rWave = QWave_i./(QRiver_i.*FF0+QWave_i+QTide_i);
% % % rTide = QTide_i./(QRiver_i.*FF0+QWave_i+QTide_i);
% % % 
% % % y0 = rRiver.*sind(60);
% % % x0 = rTide + rRiver.*sind(60).*cotd(60);
% % % % % % yM = yM_WorldEx_Gall(delta);
% % % % % % xM = xM_WorldEx_Gall(delta);
% % % y_obs = y_obs_all;
% % % x_obs = x_obs_all;
% % % 
% % % FF = nan(length(QRiver_i),1);
% % % WF = nan(length(QWave_i),1);
% % % TF = nan(length(QTide_i),1);
% % % 
% % % % ***All factors (FF)***
% % % for i = 1:length(FF)
% % % fun_factors = @(FF_var) sqrt((x_obs(i) - ((QTide_i(i)*FF_var(3))./(QRiver_i(i).*FF_var(1)+QWave_i(i)*FF_var(2)+QTide_i(i)*FF_var(3))...
% % %     + (QRiver_i(i).*FF_var(1)./(QRiver_i(i).*FF_var(1)+QWave_i(i)*FF_var(2)+QTide_i(i)*FF_var(3))).*sind(60).*cotd(60))).^2 ...
% % %     + (y_obs(i) - ((QRiver_i(i).*FF_var(1)./(QRiver_i(i)*FF_var(1)+QWave_i(i)*FF_var(2)+QTide_i(i)*FF_var(3))).*sind(60))).^2);
% % % options = optimset('TolFun',1e-8);
% % % % % % SR = fminsearch(fun_SR, SR0, options);
% % % 
% % % FF_varsol = ...
% % %     fminsearch(fun_factors, [0,0,0], options);
% % % FF(i,1) = FF_varsol(1);
% % % WF(i,1) = FF_varsol(2);
% % % TF(i,1) = FF_varsol(3);
% % % end
% % % 
% % % rRiver2 = QRiver_i.*FF./(QRiver_i.*FF+QWave_i.*WF+QTide_i.*TF);
% % % rWave2 = QWave_i.*WF./(QRiver_i.*FF+QWave_i.*WF+QTide_i.*TF);
% % % rTide2 = QTide_i.*TF./(QRiver_i.*FF+QWave_i.*WF+QTide_i.*TF);
% % % 
% % % rRiver2 + rWave2 + rTide2
% % % 
% % % y2 = rRiver2.*sind(60);
% % % x2 = rTide2 + rRiver2.*sind(60).*cotd(60);
% % % 
% % % Error_rx2 = sqrt((x_obs-x2).^2+ (y_obs-y2).^2);
% % % 
% % % 
% % % 
% % % 
% % % figure
% % % plot(index_riverdom_obs, FF(index_riverdom_obs), ...
% % %     's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g',...
% % %     'MarkerSize',15), hold on
% % % plot(index_tidedom_obs, FF(index_tidedom_obs), ...
% % %     'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r',...
% % %     'MarkerSize',12)
% % % plot(index_wavedom_obs, FF(index_wavedom_obs), ...
% % %     'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b',...
% % %     'MarkerSize',12)
% % % set(gca,'XLim', [0 32], 'YLim', [1e-9 1e-2])
% % % set(gca,'XTick', 0:1:31, 'YTick', [1e-6 1e-3 1e0 1e3])
% % % set(gca,'XTickLabel',{'' '1' '' '' '' '' '6' '' '' '' '' '11' '' '' '' '' '16' '' '' '' '' '21' '' '' '' '' '26' '' '' '' '' '31'})
% % % set(gca,'FontSize',25)
% % % set(gca,'YScale','log')
% % % xtickangle(0)
% % % xlabel('Delta','FontSize',25,'FontName','Arial')
% % % ylabel('Fluvial Factor, $FF$','Interpreter','latex',...
% % %     'FontSize',30,'FontName','Arial')
% % % grid on, box on
% % % 
% % % % % % set(gcf,'PaperPositionMode','auto')
% % % % % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % % % % print(gcf,'-dpng','-r600','y_Fig06_MorphBasedRatios_FluvialFactor')
% % % % % % close(gcf)
% % % 
% % % 
% % % 
% % % % Wave Factor
% % % figure
% % % plot(index_riverdom_obs, WF(index_riverdom_obs), ...
% % %     's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g',...
% % %     'MarkerSize',15), hold on
% % % plot(index_tidedom_obs, WF(index_tidedom_obs), ...
% % %     'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r',...
% % %     'MarkerSize',12)
% % % plot(index_wavedom_obs, WF(index_wavedom_obs), ...
% % %     'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b',...
% % %     'MarkerSize',12)
% % % set(gca,'XLim', [0 32], 'YLim', [1e-9 1e-2])
% % % set(gca,'XTick', 0:1:31, 'YTick', [1e-6 1e-3 1e0 1e3])
% % % set(gca,'XTickLabel',{'' '1' '' '' '' '' '6' '' '' '' '' '11' '' '' '' '' '16' '' '' '' '' '21' '' '' '' '' '26' '' '' '' '' '31'})
% % % set(gca,'FontSize',25)
% % % set(gca,'YScale','log')
% % % xtickangle(0)
% % % xlabel('Delta','FontSize',25,'FontName','Arial')
% % % ylabel('Wave Factor, $WF$','Interpreter','latex',...
% % %     'FontSize',30,'FontName','Arial')
% % % grid on, box on
% % % 
% % % % % % set(gcf,'PaperPositionMode','auto')
% % % % % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % % % % print(gcf,'-dpng','-r600','y_Fig06_MorphBasedRatios_WaveFactor')
% % % % % % close(gcf)
% % % 
% % % 
% % % % Tide Factor
% % % figure
% % % plot(index_riverdom_obs, TF(index_riverdom_obs), ...
% % %     's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g',...
% % %     'MarkerSize',15), hold on
% % % plot(index_tidedom_obs, TF(index_tidedom_obs), ...
% % %     'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r',...
% % %     'MarkerSize',12)
% % % plot(index_wavedom_obs, TF(index_wavedom_obs), ...
% % %     'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b',...
% % %     'MarkerSize',12)
% % % set(gca,'XLim', [0 32], 'YLim', [1e-9 1e-2])
% % % set(gca,'XTick', 0:1:31, 'YTick', [1e-8 1e-4 1e0 1e4])
% % % set(gca,'XTickLabel',{'' '1' '' '' '' '' '6' '' '' '' '' '11' '' '' '' '' '16' '' '' '' '' '21' '' '' '' '' '26' '' '' '' '' '31'})
% % % set(gca,'FontSize',25)
% % % set(gca,'YScale','log')
% % % xtickangle(0)
% % % xlabel('Delta','FontSize',25,'FontName','Arial')
% % % ylabel('Tide Factor, $WF$','Interpreter','latex',...
% % %     'FontSize',30,'FontName','Arial')
% % % grid on, box on

% % % set(gcf,'PaperPositionMode','auto')
% % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % print(gcf,'-dpng','-r600','y_Fig06_MorphBasedRatios_TideFactor')
% % % close(gcf)

%% OCEAN FACTORS (TF AND WF) TO INCREASE ACCURACY
% % % QRiver_all = [QRiver_NorthAmerica; QRiver_SouthAmerica; QRiver_AfricaEurope; QRiver_AsiaOceania];
% % % QTide_all = [QTide_NorthAmerica; QTide_SouthAmerica; QTide_AfricaEurope; QTide_AsiaOceania];
% % % QWave_all = [QWave_NorthAmerica; QWave_SouthAmerica; QWave_AfricaEurope; QWave_AsiaOceania];
% % % 
% % % y_obs_all = [y_obs_NorthAmerica; y_obs_SouthAmerica; y_obs_AfricaEurope; y_obs_AsiaOceania];
% % % x_obs_all = [x_obs_NorthAmerica; x_obs_SouthAmerica; x_obs_AfricaEurope; x_obs_AsiaOceania];
% % % % % % delta = 4;
% % % % % % QRiver_i = QRiver_all(delta);
% % % % % % QTide_i = QTide_all(delta);
% % % % % % QWave_i = QWave_all(delta);
% % % QRiver_i = QRiver_all;
% % % QTide_i = QTide_all;
% % % QWave_i = QWave_all;
% % % 
% % % FF0 = zeros(length(QRiver_i),1);
% % % rRiver = QRiver_i.*FF0./(QRiver_i.*FF0+QWave_i+QTide_i);
% % % rWave = QWave_i./(QRiver_i.*FF0+QWave_i+QTide_i);
% % % rTide = QTide_i./(QRiver_i.*FF0+QWave_i+QTide_i);
% % % 
% % % y0 = rRiver.*sind(60);
% % % x0 = rTide + rRiver.*sind(60).*cotd(60);
% % % % % % yM = yM_WorldEx_Gall(delta);
% % % % % % xM = xM_WorldEx_Gall(delta);
% % % y_obs = y_obs_all;
% % % x_obs = x_obs_all;
% % % 
% % % WF = nan(length(QWave_i),1);
% % % TF = nan(length(QTide_i),1);
% % % 
% % % % ***Ocean Factors***
% % % for i = 1:length(WF)
% % % fun_factors = @(FF_var) sqrt((x_obs(i) - ((QTide_i(i)*FF_var(2))./(QRiver_i(i)+QWave_i(i)*FF_var(1)+QTide_i(i)*FF_var(2))...
% % %     + (QRiver_i(i)./(QRiver_i(i)+QWave_i(i)*FF_var(1)+QTide_i(i)*FF_var(2))).*sind(60).*cotd(60))).^2 ...
% % %     + (y_obs(i) - ((QRiver_i(i)./(QRiver_i(i)+QWave_i(i)*FF_var(1)+QTide_i(i)*FF_var(2))).*sind(60))).^2);
% % % options = optimset('TolFun',1e-8);
% % % % % % SR = fminsearch(fun_SR, SR0, options);
% % % 
% % % FF_varsol = ...
% % %     fminsearch(fun_factors, [0,0], options);
% % % WF(i,1) = FF_varsol(1);
% % % TF(i,1) = FF_varsol(2);
% % % end
% % % 
% % % rRiver2 = QRiver_i./(QRiver_i+QWave_i.*WF+QTide_i.*TF);
% % % rWave2 = QWave_i.*WF./(QRiver_i+QWave_i.*WF+QTide_i.*TF);
% % % rTide2 = QTide_i.*TF./(QRiver_i+QWave_i.*WF+QTide_i.*TF);
% % % 
% % % rRiver2 + rWave2 + rTide2
% % % 
% % % y2 = rRiver2.*sind(60);
% % % x2 = rTide2 + rRiver2.*sind(60).*cotd(60);
% % % 
% % % Error_rx2 = sqrt((x_obs-x2).^2+ (y_obs-y2).^2);
% % % 
% % % 
% % % 
% % % 
% % % % Wave Factor
% % % figure
% % % plot(index_riverdom_obs, WF(index_riverdom_obs), ...
% % %     's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g',...
% % %     'MarkerSize',15), hold on
% % % plot(index_tidedom_obs, WF(index_tidedom_obs), ...
% % %     'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r',...
% % %     'MarkerSize',12)
% % % plot(index_wavedom_obs, WF(index_wavedom_obs), ...
% % %     'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b',...
% % %     'MarkerSize',12)
% % % set(gca,'XLim', [0 32], 'YLim', [1e-2 1e4])
% % % set(gca,'XTick', 0:1:31, 'YTick', [1e-2 1e-0 1e2 1e4])
% % % set(gca,'XTickLabel',{'' '1' '' '' '' '' '6' '' '' '' '' '11' '' '' '' '' '16' '' '' '' '' '21' '' '' '' '' '26' '' '' '' '' '31'})
% % % set(gca,'FontSize',25)
% % % set(gca,'YScale','log')
% % % xtickangle(0)
% % % xlabel('Delta','FontSize',25,'FontName','Arial')
% % % ylabel('Wave Factor, $WF$','Interpreter','latex',...
% % %     'FontSize',30,'FontName','Arial')
% % % grid on, box on
% % % 
% % % % % % set(gcf,'PaperPositionMode','auto')
% % % % % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % % % % print(gcf,'-dpng','-r600','y_Fig06_MorphBasedRatios_WaveFactor')
% % % % % % close(gcf)
% % % 
% % % 
% % % % Tide Factor
% % % figure
% % % plot(index_riverdom_obs, TF(index_riverdom_obs), ...
% % %     's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g',...
% % %     'MarkerSize',15), hold on
% % % plot(index_tidedom_obs, TF(index_tidedom_obs), ...
% % %     'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r',...
% % %     'MarkerSize',12)
% % % plot(index_wavedom_obs, TF(index_wavedom_obs), ...
% % %     'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b',...
% % %     'MarkerSize',12)
% % % set(gca,'XLim', [0 32], 'YLim', [1e-2 1e4])
% % % set(gca,'XTick', 0:1:31, 'YTick', [1e-2 1e0 1e2 1e4])
% % % set(gca,'XTickLabel',{'' '1' '' '' '' '' '6' '' '' '' '' '11' '' '' '' '' '16' '' '' '' '' '21' '' '' '' '' '26' '' '' '' '' '31'})
% % % set(gca,'FontSize',25)
% % % set(gca,'YScale','log')
% % % xtickangle(0)
% % % xlabel('Delta','FontSize',25,'FontName','Arial')
% % % ylabel('Tidal Factor, $TF$','Interpreter','latex',...
% % %     'FontSize',30,'FontName','Arial')
% % % grid on, box on

% % % set(gcf,'PaperPositionMode','auto')
% % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % print(gcf,'-dpng','-r600','y_Fig06_MorphBasedRatios_TideFactor')
% % % close(gcf)
%}
%% ERROR in r_x (e_ter)

Error_rx_NorthAmerica = sqrt((x_obs_NorthAmerica-x_NorthAmerica).^2 ...
    + (y_obs_NorthAmerica-y_NorthAmerica).^2);

Error_rx_SouthAmerica = sqrt((x_obs_SouthAmerica-x_SouthAmerica).^2 ...
    + (y_obs_SouthAmerica-y_SouthAmerica).^2);

Error_rx_AfricaEurope = sqrt((x_obs_AfricaEurope-x_AfricaEurope).^2 ...
    + (y_obs_AfricaEurope-y_AfricaEurope).^2);

Error_rx_AsiaOceania = sqrt((x_obs_AsiaOceania-x_AsiaOceania).^2 ...
    + (y_obs_AsiaOceania-y_AsiaOceania).^2);



Error_rx = [Error_rx_NorthAmerica; Error_rx_SouthAmerica; Error_rx_AfricaEurope;...
    Error_rx_AsiaOceania];


figure
plot(index_riverdom_obs, Error_rx1(index_riverdom_obs), ...
    '*', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.5 0.5 0.5],...
    'MarkerSize',17), hold on
plot(index_tidedom_obs, Error_rx1(index_tidedom_obs), ...
    '*', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.5 0.5 0.5],...
    'MarkerSize',17)
plot(index_wavedom_obs, Error_rx1(index_wavedom_obs), ...
    '*', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.5 0.5 0.5],...
    'MarkerSize',17)
% % % plot(index_riverdomM, Error_rx2(index_riverdomM), ...
% % %     '^', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.5 0.5 0.5],...
% % %     'MarkerSize',17), hold on
% % % plot(index_tidedomM, Error_rx2(index_tidedomM), ...
% % %     '^', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.5 0.5 0.5],...
% % %     'MarkerSize',17)
% % % plot(index_wavedomM, Error_rx2(index_wavedomM), ...
% % %     '^', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.5 0.5 0.5],...
% % %     'MarkerSize',17)
plot(index_riverdom_obs, Error_rx(index_riverdom_obs), ...
    's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g',...
    'MarkerSize',15), hold on
plot(index_tidedom_obs, Error_rx(index_tidedom_obs), ...
    'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r',...
    'MarkerSize',12)
plot(index_wavedom_obs, Error_rx(index_wavedom_obs), ...
    'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b',...
    'MarkerSize',12)
set(gca,'XLim', [0 32], 'YLim', [0 1])
set(gca,'XTick', 0:1:31, 'YTick', 0:0.25:1)
set(gca,'XTickLabel',{'' '1' '' '' '' '' '6' '' '' '' '' '11' '' '' '' '' '16' '' '' '' '' '21' '' '' '' '' '26' '' '' '' '' '31'})
set(gca,'FontSize',25)
xtickangle(0)
xlabel('Delta','FontSize',25,'FontName','Arial')
ylabel('$e_{ter}$','Interpreter','latex',...
    'FontSize',30,'FontName','Arial')
grid on, box on

% % % set(gcf,'PaperPositionMode','auto')
% % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % print(gcf,'-dpng','-r600','y_Fig06_MorphBasedRatios_Error')
% % % close(gcf)

%% SEDIMENT RETENTION
QRiver_ret_all = [QRiver_ret_NorthAmerica; QRiver_ret_SouthAmerica; QRiver_ret_AfricaEurope; ...
    QRiver_ret_AsiaOceania];
% % % QWave_obs_all = abs([QWave_obs_NorthAmerica; QWave_obs_SouthAmerica; QWave_obs_AfricaEurope; ...
% % %     QWave_obs_AsiaOceania]);
% % % QTide_obs_all = [QTideM_NorthAmerica; QTide_obs_SouthAmerica; QTide_obs_AfricaEurope; ...
% % %     QTide_obs_AsiaOceania];


% % % SR = (QRiver_obs_all+QWave_obs_all+QTide_obs_all)...
% % %     ./(QRiver_all+QWave_all+QTide_all);

SR = QRiver_ret_all./QRiver_all



figure
plot(0:32,(0:32).*0+1, '-k', 'Linewidth', 3), hold on
plot(index_riverdom_obs, SR(index_riverdom_obs), ...
    's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g',...
    'MarkerSize',15)
plot(index_tidedom_obs, SR(index_tidedom_obs), ...
    'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r',...
    'MarkerSize',12)
plot(index_wavedom_obs, SR(index_wavedom_obs), ...
    'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b',...
    'MarkerSize',12)
set(gca,'XLim', [0 32], 'YLim', [1e-4 1e2])
set(gca,'XTick', 0:1:31, 'YTick', [1e-4 1e-2 1e0])
set(gca, 'YScale', 'log')
set(gca,'XTickLabel',{'' '1' '' '' '' '' '6' '' '' '' '' '11' '' '' '' '' '16' '' '' '' '' '21' '' '' '' '' '26' '' '' '' '' '31'})
set(gca,'FontSize',25)
xtickangle(0)
xlabel('Delta','FontSize',25,'FontName','Arial')
ylabel('Sediment Retention, $SR$','Interpreter','latex',...
    'FontSize',30,'FontName','Arial')
grid on, box on

% % % set(gcf,'PaperPositionMode','auto')
% % % % % % set(gcf,'PaperUnits','inches','PaperPosition',[0.5 0.5 figwidth figheight])
% % % print(gcf,'-dpng','-r600','y_Fig06_MorphBasedRatios_SedRetention')
% % % close(gcf)

%% Table of Morphology-Based Ratios
close all, clc

table_NorthAmerica = table(Deltaname_NorthAmerica, MouthLat(index_NorthAmerica), ...
    MouthLon(index_NorthAmerica), ...
    QRiver_NorthAmerica, QWave_NorthAmerica, QTide_NorthAmerica, ...
    rRiver_NorthAmerica, rWave_NorthAmerica, rTide_NorthAmerica, ...
    R_NorthAmerica, T_NorthAmerica, ...
    Ndist_pred_NorthAmerica, wm_wu_pred_NorthAmerica, theta_pred_NorthAmerica,...
    x_NorthAmerica, y_NorthAmerica, Error_rx_NorthAmerica, FF(index_WorldEx_total_NA), Error_rx1_NorthAmerica, ...
    QRiver_ret_NorthAmerica, QWave_NorthAmerica, abs(Qw_obs_NorthAmerica), QTide_obs_NorthAmerica, ...
    rRiver_obs_NorthAmerica, rWave_obs_NorthAmerica, rTide_obs_NorthAmerica, ...
    R_obs_NorthAmerica, T_obs_NorthAmerica, ...
    Ndist_obs_NorthAmerica, wm_NorthAmerica./wu_NorthAmerica, (theta_leftobs_NorthAmerica+theta_rightobs_NorthAmerica)/2, ...
    x_obs_NorthAmerica, y_obs_NorthAmerica, SR(index_WorldEx_total_NA));



table_SouthAmerica = table(Deltaname_SouthAmerica, MouthLat(index_SouthAmerica), ...
    MouthLon(index_SouthAmerica), ...
    QRiver_SouthAmerica, QWave_SouthAmerica, QTide_SouthAmerica, ...
    rRiver_SouthAmerica, rWave_SouthAmerica, rTide_SouthAmerica, ...
    R_SouthAmerica, T_SouthAmerica, ...
    Ndist_pred_SouthAmerica, wm_wu_pred_SouthAmerica, theta_pred_SouthAmerica,...
    x_SouthAmerica, y_SouthAmerica, Error_rx_SouthAmerica, FF(index_WorldEx_total_SA), Error_rx1_SouthAmerica, ...
    QRiver_ret_SouthAmerica, QWave_SouthAmerica, abs(Qw_obs_SouthAmerica), QTide_obs_SouthAmerica, ...
    rRiver_obs_SouthAmerica, rWave_obs_SouthAmerica, rTide_obs_SouthAmerica, ...
    R_obs_SouthAmerica, T_obs_SouthAmerica, ...
    Ndist_obs_SouthAmerica, wm_SouthAmerica./wu_SouthAmerica, (theta_leftobs_SouthAmerica+theta_rightobs_SouthAmerica)/2, ...
    x_obs_SouthAmerica, y_obs_SouthAmerica, SR(index_WorldEx_total_SA));



table_AfricaEurope = table(Deltaname_AfricaEurope, MouthLat(index_AfricaEurope), ...
    MouthLon(index_AfricaEurope), ...
    QRiver_AfricaEurope, QWave_AfricaEurope, QTide_AfricaEurope, ...
    rRiver_AfricaEurope, rWave_AfricaEurope, rTide_AfricaEurope, ...
    R_AfricaEurope, T_AfricaEurope, ...
    Ndist_pred_AfricaEurope, wm_wu_pred_AfricaEurope, theta_pred_AfricaEurope,...
    x_AfricaEurope, y_AfricaEurope, Error_rx_AfricaEurope, FF(index_WorldEx_total_AE), Error_rx1_AfricaEurope, ...
    QRiver_ret_AfricaEurope, QWave_AfricaEurope, abs(Qw_obs_AfricaEurope), QTide_obs_AfricaEurope, ...
    rRiver_obs_AfricaEurope, rWave_obs_AfricaEurope, rTide_obs_AfricaEurope, ...
    R_obs_AfricaEurope, T_obs_AfricaEurope, ...
    Ndist_obs_AfricaEurope, wm_AfricaEurope./wu_AfricaEurope, (theta_leftobs_AfricaEurope+theta_rightobs_AfricaEurope)/2, ...
    x_obs_AfricaEurope, y_obs_AfricaEurope, SR(index_WorldEx_total_AE));



table_AsiaOceania = table(Deltaname_AsiaOceania, MouthLat(index_AsiaOceania), ...
    MouthLon(index_AsiaOceania), ...
    QRiver_AsiaOceania, QWave_AsiaOceania, QTide_AsiaOceania, ...
    rRiver_AsiaOceania, rWave_AsiaOceania, rTide_AsiaOceania, ...
    R_AsiaOceania, T_AsiaOceania, ...
    Ndist_pred_AsiaOceania, wm_wu_pred_AsiaOceania, theta_pred_AsiaOceania,...
    x_AsiaOceania, y_AsiaOceania, Error_rx_AsiaOceania, FF(index_WorldEx_total_AO), Error_rx1_AsiaOceania, ...
    QRiver_ret_AsiaOceania, QWave_AsiaOceania, abs(Qw_obs_AsiaOceania), QTide_obs_AsiaOceania, ...
    rRiver_obs_AsiaOceania, rWave_obs_AsiaOceania, rTide_obs_AsiaOceania, ...
    R_obs_AsiaOceania, T_obs_AsiaOceania, ...
    Ndist_obs_AsiaOceania, wm_AsiaOceania./wu_AsiaOceania, (theta_leftobs_AsiaOceania+theta_rightobs_AsiaOceania)/2, ...
    x_obs_AsiaOceania, y_obs_AsiaOceania, SR(index_WorldEx_total_AO));


writetable(table_NorthAmerica,'table_NorthAmerica.xlsx') 
writetable(table_SouthAmerica,'table_SouthAmerica.xlsx') 
writetable(table_AfricaEurope,'table_AfricaEurope.xlsx') 
writetable(table_AsiaOceania,'table_AsiaOceania.xlsx') 