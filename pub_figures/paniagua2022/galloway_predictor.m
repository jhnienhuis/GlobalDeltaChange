function [Qwave,Qtide,Qriver,shoreline_angle,w_upstream,w_downstream] = galloway_predictor(hs,a0,Qriver,discharge);


%% input
%hs = 1.5; %significant wave height m
tp = 5; %wave period s
wdir = [0]; %wave approach angles (deg)
rel_contrib = [1]; %relative contribution of waves

%discharge = 1000; %river discharge m3s-1

%a0 = 1.5; %tidal amplitude m
omega0 = 1.4e-4; %tidal frequency in rad/s (semidiurnal = 1.4e-4)
rel_dens = 1.65; %immersed density
d50 = 2e-4; %grain size in m
g = 9.81; %gravity
chezy = 55; %roughness

%% qriver
%the fluvial sediment flux. try to estimate the fraction that is retained
%nearshore. In the 2015 geology paper we simply used the bedload fraction.
%In the 2020 nature paper we used bedload and suspended load because we
%included tides and because bedload is not available at a global scale. 

%probably best to simply use sus+bed load flux, and to get it directly from
%the delft3d simulation. 
%Qriver = 100; %kg/s

%% qwave
energy = zeros(180,1); %array of wave energy approach angle divided over 180 approach angles (i.e. left to right looking offshore)
energy(wdir+90) = rel_contrib./sum(rel_contrib).*hs^2.4.*tp^0.2; 

%CERC like equation to convert wave energy to sediment flux (see Nienhuis 2015 Geology and Ashton & Murray 2006 JGR)
AngArray = linspace(-0.5*pi,0.5*pi,180);
k_wave = 5.3e-6*1050*(g^1.5)*(0.5^1.2)*(sqrt(g*0.78)/(2*pi))^0.2;
sedtrans = 2650 * (1-0.4) * k_wave .*  (cos(AngArray).^1.2) .* sin(AngArray); 

sedconv = conv(fliplr(energy),sedtrans,'same'); % sedconv in kg/s
Qwave = (max(sedconv)-min(sedconv)); %maximum potential flux away from the river mouth

Qwave_simple = 2*2650 * (1-0.4) * k_wave *hs^(12/5).*tp^0.2.*0.47;

%% qtide
%see also this table that is a supplementary table for the 2018 GRL paper 
%https://www.dropbox.com/s/uk94kkqcb1vnbbo/Nienhuis_TidalDeltas_S1_New.xlsx?dl=0

w_upstream = 6*discharge^0.5; %channel width upstream (m), example hydraulic geometry

qt = Qriver/1600/w_upstream; %m2/s
cf = g/(chezy^2); %non dimensional friction from chezy
rb = (cf*(discharge/w_upstream)^2/g)^(1/3);
slope = 5e-5; %((qt / (sqrt(rel_dens*g*d50)*d50)/(0.05/cf)).^(1/2.5)*rel_dens*d50/rb).^1.5; %slope estimate based on equilibrium normal flow profile from gary parker ebook ch 14.
d_upstream = (cf*(discharge/w_upstream)^2/g/slope)^(1/3); %water depth at bankful

beta = 100; %w_upstream/d_upstream; %channel_aspect_ratio at mouth

k_tide = omega0/(sqrt(d50/0.2)*rel_dens*chezy*pi); %cross-sectional area relationship k, see nienhuis GRL 2018
d_length = d_upstream / slope; %tidal intrusion length estimate (m)

Qtide_w = a0^2*k_tide*d_length^2*beta*omega0*0.5+omega0*a0*d_length*w_upstream; %qtide in water discharge m3/s
Qtide = Qtide_w*Qriver / discharge; %qtide in kg/s assuming sed concentration of tidal transport = sed concentration of river

%% combination:

R = Qriver / Qwave; %R>1 is river dominated, R<1 is wave-dominated
T = Qtide / Qriver; %T<1 is river dominated, T>1 is tide-dominated

%ternplot(Qtide,Qriver,Qwave,'scatter','filled')

%% morphologic prediction
%river mouth width in m (needs small correction if there are multiple river mouths, see GRL 2018
w_downstream = (a0*k_tide*d_length*beta)+w_upstream; 

%shoreline angle (deg) compared to reference shoreline assuming a balance
%between Qriver and actual (not potential!) wave-driven transport
if R<1, 
    sedconv = sedconv - sedconv(90); %correct for net alongshore transport from asymmetric wave climate
    if max(sedconv)<(0.5*Qriver), %right flank at maximum, river will migrate to the left, see Nienhuis 2016 EPSL
        [~,shoreline_angle(1)] = max(sedconv);
        shoreline_angle(1) = shoreline_angle(1)-90;
        shoreline_angle(2) = find(sedconv(90:-1:1)<(-(Qriver-max(sedconv))),1,'first'); %left flank of the wave-dominated delta
        
    elseif (-1 * min(sedconv)) < (0.5*Qriver), %left flank at maximum, river will migrate towards the right
        [~,shoreline_angle(2)] = min(sedconv);
        shoreline_angle(2) = 90-shoreline_angle(1);
        shoreline_angle(1) = find(sedconv(90:end)<(Qriver+min(sedconv)),1,'first'); %left flank of the wave-dominated delta
        
    else,
        shoreline_angle(1) = find(sedconv(90:end)>(0.5*Qriver),1,'first'); %right flank of the wave-dominated delta
        shoreline_angle(2) = find(sedconv(90:-1:1)<(-0.5*Qriver),1,'first'); %left flank of the wave-dominated delta
    end
else,
    shoreline_angle = nan; %river dominated deltas have no characteristic shoreline angle at the river mouth
end

shoreline_angle = mean(shoreline_angle);
