function [Ndist_obs, Delta_theta, theta_right, theta_left, R_obs, Qw_obs, Qw_obs_left, Qw_obs_right, ...
    R, Ndist_pred, theta_pred, theta_maxQw, theta_minQw, ...
    T, w_u, w_m, wm_wu_pred, T_obs, QTide_obs, ...
    QRiver_ret] = ...
    jfpa_MorphologyPred(QRiver, qRiver, QWave, QTide, S, fileNameGE, fileNameWW3, fileNameTides, GEdatapath, WW3datapath, Tidesdatapath)

    % [kml2struct function by James Slegers, 2012]
    % Order of coordinates in 'struct' files:
    % 1 - end point, 2 - starting point
    kmlStruct = kml2struct([GEdatapath fileNameGE '.kml']);
    shore_lon = kmlStruct(1).Lon;
    shore_lat = kmlStruct(1).Lat;
    wu_lon = kmlStruct(3).Lon;
    wu_lat = kmlStruct(3).Lat;
    LeftFlank_lon = kmlStruct(4).Lon;
    LeftFlank_lat = kmlStruct(4).Lat;
    RightFlank_lon = kmlStruct(5).Lon;
    RightFlank_lat = kmlStruct(5).Lat;
    wm_lon = nan(2*(length(kmlStruct)-5),1);
    wm_lat = nan(2*(length(kmlStruct)-5),1);
    Ndist_obs = length(wm_lon)/2;
    for ii = 1:Ndist_obs
        wm_lon(2*ii-1:2*ii,1) = kmlStruct(5+ii).Lon;
        wm_lat(2*ii-1:2*ii,1) = kmlStruct(5+ii).Lat;
    end
% % %     profile_lon = kmlStruct(6).Lon;
% % %     profile_lat = kmlStruct(6).Lat;

% % %     figure
% % %     plot(LeftFlank_lon,LeftFlank_lat,'-b'), hold on
% % %     plot(RightFlank_lon,RightFlank_lat,'-g')
% % %     plot(wm_lon,wm_lat,'-r')
% % %     plot(wu_lon,wu_lat,'-c')
% % %     axis equal

    % Calculation of flank angle, theta [Nienhuis et al. 2015 Geology 43]
    % Azimuthal angles (positive clockwise with 0° at N)
    % 
    % Left flank (distances in meters)
    LeftFlank_dx = (LeftFlank_lon(1) - LeftFlank_lon(2))*111300*cosd(LeftFlank_lat(1));
    LeftFlank_dy = (LeftFlank_lat(1) - LeftFlank_lat(2))*111300;
    LeftFlank_theta = atan2d(LeftFlank_dx,LeftFlank_dy);
    if LeftFlank_theta<0
        LeftFlank_theta = LeftFlank_theta+360;
    end

    % Right flank
    RightFlank_dx = (RightFlank_lon(1) - RightFlank_lon(2))*111300*cosd(RightFlank_lat(1));
    RightFlank_dy = (RightFlank_lat(1) - RightFlank_lat(2))*111300;
    RightFlank_theta = atan2d(RightFlank_dx,RightFlank_dy);
    if RightFlank_theta<0
        RightFlank_theta = RightFlank_theta+360;
    end
    
    % Difference in theta angles
    Delta_theta = RightFlank_theta - LeftFlank_theta;
    if Delta_theta<0
        Delta_theta = Delta_theta+360;
    end
    
    % Total fluvial dominance
    if Delta_theta>180
        Delta_theta = 180;
    end
    
    
    % General shoreline orientation, theta_shore
    % (in azimuthal degrees) - as in Nienhuis et al., 2015.
    % Reference shoreline orientation in degrees (measure in google earth the
    %angle along a straight -non deltaic- reference coastline, the angle from
    %right to left). a zero degree coastline runs from north to south with the
    %sea on the east side.
    shore_dx = (shore_lon(1) - shore_lon(2))*111300*cosd(shore_lat(1));
    shore_dy = (shore_lat(1) - shore_lat(2))*111300;
    theta_shore = round(atan2d(shore_dx, shore_dy));
    if theta_shore<0
        theta_shore = theta_shore + 360;
    end
    
    % WaveWatch3 data...
    dp = ncread([WW3datapath fileNameWW3],'dp'); % Peak direction, azimuthal degrees
    dp = double(dp)/100;
    dp = dp(~isnan(dp));
    hs = ncread([WW3datapath fileNameWW3],'hs'); % Significant wave height, mm 
    hs = double(hs)/1000;  % convert mm to m
    hs = hs(~isnan(hs));
    tp = ncread([WW3datapath fileNameWW3],'tp'); % Peak period, ms 
    tp = double(tp)/1000;  % convert ms to s
    tp = tp(~isnan(tp));
    
    % Q_wave values from wave climate and general shoreline orientation
    [~, Qw_net, Qw_uni, E, E_shore, theta] = ...
        jfpa_Qwave(theta_shore, dp, hs, tp);

% % %     % Wave rose
% % %     Options = {'anglenorth',0,'angleeast',90,'labels',...
% % %         {'','','',''},...
% % %         'ndirections',36,...
% % %         'nSpeeds',2,'vwinds',[0 100],...
% % %         'MaxFrequency',20,'nFreq',20,...
% % %         'LegendVariable','',...
% % %         'cMap','white',...
% % %         'scalefactor',0.5,...
% % %         'min_radius',0};
% % %     [~, ~, ~, ~, ~] = jfpa_WaveRose(dp, E, Options);
% % % 
% % %     % Energy distribution
% % %     figure
% % %     plot(theta, Qw_net)
% % %     
% % %     figure
% % %     bar(theta, E_shore)
% % %     
% % %     figure
% % %     plot(theta, Qw_uni)
    

    % Flank angles with respect to general shoreline
    % theta_left, theta_right
    % Check notes for positive/negative angles.
    % Delta flanks digitized towards river mouth.
    delta_thetaL = LeftFlank_theta - theta_shore;

    if delta_thetaL<=0
        delta_thetaL = delta_thetaL + 360;
    end
    
    % Calculate per quadrant...
    % Multi-channel
    if delta_thetaL>=0 && delta_thetaL<90
        theta_left = 0; % arbitrary
        thetaL_factor = Inf;
        
    % Cuspate
    elseif delta_thetaL>=90 && delta_thetaL<180
        theta_left = round(180 - delta_thetaL);
        thetaL_factor = 1;

    % Estuary
    elseif delta_thetaL>=180 && delta_thetaL<270
        theta_left = round(delta_thetaL - 180);
        thetaL_factor = -1;
        
% % %     % Spit (not used in Google Earth)
% % %     elseif delta_thetaL>=270 && delta_thetaL<360
% % %         theta_left = 0;
% % %         thetaL_factor = Inf;

    end
    
    
    
    delta_thetaR = RightFlank_theta - theta_shore;
    
    if delta_thetaR<=0
        delta_thetaR = delta_thetaR + 360;
    end
    
    % Calculate per quadrant...
    % Cuspate
    if delta_thetaR>=0 && delta_thetaR<90
        theta_right = round(delta_thetaR);
        thetaR_factor = -1;
        
    % Multi-channel
    elseif delta_thetaR>=90 && delta_thetaR<180
        theta_right = 0; % arbitrary
        thetaR_factor = Inf;
        
    % Spit (not used in Google Earth)
% % %     elseif delta_thetaR>=180 && delta_thetaR<270
% % %         theta_right = 0;
% % %         thetaR_factor = Inf;
        
    % Estuary
    elseif delta_thetaR>=270 && delta_thetaR<360
        theta_right = round(360 - delta_thetaR);
        thetaR_factor = 1;
        
    end

    % Net sediment flux at the mouth, as the difference between delta
    % flanks. Positive means net flux to the right of general shoreline.
    Qw_obs_left = thetaL_factor*Qw_net(theta==(theta_left*thetaL_factor));
    Qw_obs_right = thetaR_factor*Qw_net(theta==(theta_right*thetaR_factor));
    % QwaveM_right: negative add to mouth, positive take from mouth
    % QwaveM_left: negative take from mouth, positive add to mouth
    if isempty(Qw_obs_left)==1
        if isempty(Qw_obs_right)==1
            Qw_obs = NaN;
        else
            Qw_obs = Qw_obs_right;
        end
        
    elseif isempty(Qw_obs_right)==1
        Qw_obs = Qw_obs_left;
            
        
    elseif isinf(Qw_obs_left)==1
        if isinf(Qw_obs_right)==1
            Qw_obs = NaN;
        else
            Qw_obs = Qw_obs_right;
        end
        
    else
        if isinf(Qw_obs_right)==1
            Qw_obs = Qw_obs_left;
        else
            Qw_obs = Qw_obs_right - Qw_obs_left;
        end
    end
    
% % %     Qwave_M = QwaveM_right - QwaveM_left;
% % %     QWave = abs(max(Qw_net)) + abs(min(Qw_net));
    R = QRiver/QWave;
    theta_pred = 0.5*asind(R);
        theta_pred(R>1) = Inf;
    
    if R<=1
        Ndist_pred = 1;
    else
        Ndist_pred = round(R);
    end
    
    theta_maxQw = theta(Qw_net==max(Qw_net));
    theta_minQw = theta(Qw_net==min(Qw_net));
    
    
    % If river dominated (multiple mouths), then QwaveM=0 as calculated above. Calculate R^M
    % by counting the distributary mouths from the Google Earth kml file:
    if isnan(Qw_obs)==1 % N_dist_obs>1 (several distributaries)
        R_obs = length(kmlStruct) - 5; % distributaries start at position 6
        theta_right = Inf; % undefined flank angles
        theta_left = Inf; % undefined flank angles
        QRiver_ret = R_obs*QWave;
        QRiver_ret = QRiver;
        Qw_obs = QWave;

    else % N_dist_obs=1 (one distributary)
        R_obs = abs(Qw_obs)/QWave;
        QRiver_ret = abs(Qw_obs);
        QRiver_ret = QRiver;
        
    end
    
    
    % T^obs
    % Mouth and upstream channel widths
    wm_dx = (wm_lon(1:2:2*Ndist_obs) - wm_lon(2:2:2*Ndist_obs)).*111300.*cosd(wm_lat(1));
    wm_dy = (wm_lat(1:2:2*Ndist_obs) - wm_lat(2:2:2*Ndist_obs)).*111300;
    w_m = sum(sqrt(wm_dx.^2 + wm_dy.^2))/max(1,Ndist_obs^(1-0.5)); % Corrected mouth width
    
    
    wu_dx = (wu_lon(1) - wu_lon(2))*111300*cosd(wu_lat(1));
    wu_dy = (wu_lat(1) - wu_lat(2))*111300;
    w_u0 = sqrt(wu_dx^2 + wu_dy^2);
    
    if w_u0 > w_m
        w_u = w_m;
    else
        w_u = w_u0;
    end
        
        
    
    
    
    
    
    
    % TIDAL DATA
    tidal = readtable([Tidesdatapath fileNameTides], 'HeaderLines',1);
    
    a_tide = table2array(tidal(:,2)); % amplitude, in meters
        a_tide(a_tide==0) = NaN;
        a_tide = a_tide(~isnan(a_tide));
    Period = table2array(tidal(:,4)); % period, in hours
        Period = Period(~isnan(a_tide));
    omega = 2*pi./(Period*3600); % angular frequency (in 1/seconds)
    
% % %     % Preallocation: tidal velocities per constituent
% % %     u_tide0 = zeros(length(T),1);
% % %     d_tide0 = zeros(length(T),1);
% % %     
% % %     % Time series of tidal velocities per constituent
% % %     g = 9.81;
% % %     C_b = 0.0025;
% % %     n = 0.03; % Pasture floodplains
% % %     k_n = 1;
% % %     for ii = 1:length(a)
% % %         [u_tide0(ii,1), d_tide0(ii,1), ~, ~, ~] = ...
% % %             jfpa_tidalcurrents(g, C_b, n, k_n, S, w_m, a(ii,1), T(ii,1));
% % %     end
    
% % %     u_tide = sqrt(sum(u_tide0.^2));
% % %     d_tide = sqrt(sum(d_tide0.^2));
    
    
    % River velocity
% % %     R_h = w_u*d/(w_u+2*d);
% % %     d = -R_h*w_u/(2*R_h-w_u);
% % %     u_river0 = (k_n/n)*R_h^(2/3)*sqrt(S);
% % %     R_h = (u_river0*n/(k_n*sqrt(S)))^(3/2);
% % %     u_river0 = q_river/(w_u*d);
% % %     d = q_river/(w_u*u_river0);
% % %     R_h = w_u*(q_river/(w_u*u_river0))/(w_u+2*(q_river/(w_u*u_river0)));

% % %     fun_uriver = @(u_river0) u_river0 - (k_n/n)*(w_u*(q_river/(w_u*u_river0))/(w_u+2*(q_river/(w_u*u_river0))))^(2/3)*sqrt(S);
% % %     u_river = fsolve(fun_uriver, 0.1);
% % %     R_h = (u_river*n/(k_n*sqrt(S)))^(3/2);
% % %     d_river = -R_h*w_u/(2*R_h-w_u);
    d_u = 0.6*(qRiver)^(1/3); % Mikhailov, 1970, Eq. 21, Table 3
    u_river = qRiver/(d_u*w_u);
    
    omega = sum(omega.*a_tide)./sum(a_tide);
    a_tide = sum(a_tide);
    

    % T^M from tidal and fluvial properties
    u_tide = 0.5*(omega.*a_tide);
    U = u_tide./(S*u_river); % [Nienhuis et al. 2018 GRL Eq. 9]
    D_50 = 0.1; % Median grain size (mm)
    Theta_s = 0.2; % Critical Shields number for D_50=0.1 mm
    % Assume: R=1.65, C=55 (Chezy roughness coefficient)
    k = omega./(sqrt(Theta_s*D_50/1000)*1.65*55*pi);
    f = 1 + 2*S./(k.*a_tide);
    T_obs = ((w_m/w_u)-1).*U.*f;
    
    
    beta = w_u/d_u;
    L = d_u./S; %(w_m-w_u)./(beta.*k.*a_tide);
    QTide = 0.5*(L^2)*beta*omega.*k.*(a_tide.^2).*f;
    QTide = QTide.*QRiver./qRiver;
    T = QTide/QRiver;
    wm_pred = beta.*a_tide.*k.*L +w_u;
    wm_wu_pred = wm_pred./w_u;%(QTide/QRiver)/(U.*f) + 1;
    QTide_obs = T_obs*QRiver;
    
end
