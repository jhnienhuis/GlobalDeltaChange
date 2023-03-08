%get qs,max
% by JFPA, Tallahassee, 2018-06-06
% After Jaap Nienhuis, FSU, 2018-06-06
% Modified by JFPA, EAFIT, 2021-02-14 and beyond
function [QWave, Qw_net, Qw_uni, E, E_shore, theta] = ...
    jfpa_Qwave(theta_shore, Dp, Hs, Tp)

%reference shoreline orientation in degrees (measure in google earth the
%angle along a straight -non deltaic- reference coastline, the angle from
%right to left). a zero degree coastline runs from north to south with the
%sea on the east side.
% % % shoreline = 30; % For Ebro delta
shore = rem(theta_shore+(0:180),360)+1;

E = zeros(360,1);
%make array of the total potential qs for every wave approach angle (m^2.4 s^0.2)
for ii = 1:360
    if isempty(Hs(floor(Dp)==ii))==0
        E(ii) = nansum((Hs(floor(Dp)==ii).^2.4).*(Tp(floor(Dp)==ii).^0.2));
    end
end

%divide by total number of observations
E = E./length(Dp);

%the wave energy approaching the reference shoreline straight on is energy(90)
E_shore = E(shore);

%if more of the coast is shadowed by for instance an island, just set some
%remaining elements of the energy vector to zero. (e.g. energy(50:60) = 0
%if there is an island hiding wave approahc angles from 50 to 60 degrees


%potential approach angles of incoming waves
% Pani 2021-03-09: make theta to go from -90 to 90 every degree, so length
% must be 181 elements.
theta = round(linspace(-90, 90, 181));

%cerc formula in deep water, without energy component because that is
%already in the energy vector (kg/s)
g = 9.81;
p = 0.4;
K = 5.3e-6*1050*(g^1.5)*(0.5^1.2)*(sqrt(g*0.78)/(2*pi))^0.2;
Qw_uni = 2650 * (1-p) * K .* (cosd(theta).^1.2) .* sind(theta);

%do convolution, basically a running sum of the cerc formula with all the
%waves approaching a particular shoreline orientation. we flip the energy
%to get the direction right. (positive is flux to the right looking
%offshore)
Qw_net = conv(fliplr(E_shore),Qw_uni,'same');

%alongshore transport along the straight reference coastline (kg/s)
% % % Qw_along = Qw_net(180);

%maximum potential transport combined for the left and the right
QWave = max(Qw_net) - min(Qw_net);

%plot convolution function (shoreline angle is positive in the clockwise
%direction)
% % % plot(180*AngArray/pi,sedconv)

