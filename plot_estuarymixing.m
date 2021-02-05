%https://www-annualreviews-org.proxy.lib.fsu.edu/doi/pdf/10.1146/annurev-fluid-010313-141302


load('D:\Dropbox\github\GlobalDeltaChange\GlobalDeltaData.mat','Discharge_dist','Discharge_tide','Discharge_prist','depth_mouth','width_mouth','delta_name');



H = depth_mouth; %mouth depth (m)

Ur = Discharge_dist./(2*width_mouth.*H); %river velocity (m/s) (at mouth?)
s_ocean = 34.5; %PSU is practical salinity unit??
beta = 7.7e-4; %another thing..
g = 9.81; %gravity

Ut = Discharge_tide./(width_mouth.*H); %tidal velocity (m/s) at mouth
Cd = 3e-3; %drag coefficient
omega=2*pi/(12.5*3600); %tidal frequency (s-1)
N0 = sqrt(beta.*g.*s_ocean./H); %buoyancy frequency (s-1)


%freshwater froude number
Frf = Ur./sqrt(beta.*g.*s_ocean.*H);

%estuary mixing coefficient
M = sqrt(Cd*Ut.^2./(omega.*N0.*H.^2));

scatter((M),(Frf))
set(gca,'Xscale','log','YScale','log'); set(gca,'XLim',[0.01 100],'YLim',[1e-4 1])
hold on
%phase_spaces

%condition for vertical mixing
Frf_mix = @(M) (M.^2./sqrt(3.4)).^(1/3);

fplot(Frf_mix,[0.01 1e2])


idx = ismember(delta_name,["Ebro","Mississippi","Amazon","Fraser"]);
table(delta_name(idx),Ut(idx),Ut(idx),Frf(idx),M(idx))
