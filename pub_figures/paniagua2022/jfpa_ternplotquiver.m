% TERNPLOT plot ternary phase diagram
%   TERNPLOT(A, B) plots ternary phase diagram for three components.  C is calculated
%      as 1 - A - B.
%
%   TERNPLOT(A, B, C) plots ternary phase data for three components A B and C.  If the values 
%       are not fractions, the values are normalised by dividing by the total.
%
%   TERNPLOT(A, B, C, LINETYPE) the same as the above, but with a user specified LINETYPE (see PLOT
%       for valid linetypes).
%   
%   NOTES
%   - An attempt is made to keep the plot close to the default plot type.  The code has been based largely on the
%     code for POLAR.       
%   - The regular TITLE and LEGEND commands work with the plot from this function, as well as incrimental plotting
%     using HOLD.  Labels can be placed on the axes using TERNLABEL
%
%   See also TERNLABEL PLOT POLAR

%       b
%      / \
%     /   \
%    c --- a 

% Author: Carl Sandrock 20020827

% To do

% Modifications

% Modifiers
% CS Carl Sandrock

function [handles1,x1,y1,handles2,x2,y2] = ...
    jfpa_ternplotquiver(A, B, C, D, E, F, fontsiz, xoffset, yoffset,...
    plottype, mksz1, mksz2, markercol1, markeredgecol1, markercol2, markeredgecol2, headsize, headangle, linwidth, lincolor, tc_axis, varargin)

majors = 4;
% % % fontsiz = 16;

if nargin < 3
    C = 1 - (A+B);
end;

if nargin>3 && ~any(strcmpi(plottype,{'scatter','plot'})),
   varargin = [{plottype},varargin];
end

[fA, fB, fC, fD, fE, fF] = fractions_quiver(A, B, C, D, E, F);

[x1, y1, x2, y2] = terncoords_quiver(fA, fB, fC, fD, fE, fF);

% % % % Sort data points in x order
% % % [x1, index_sort1] = sort(x01);
% % % y1 = y01(index_sort1);
% % % 
% % % [x2, index_sort2] = sort(x02);
% % % y2 = y02(index_sort2);
% % % 
% % % % Sort additional data as well
% % % siz = find(cellfun(@length,varargin)==length(x1));
% % % for ii=siz, varargin{ii} = varargin{ii}(index_sort1); end

% Make ternary axes
[hold_state, cax, next] = ternaxes(majors,fontsiz,xoffset,yoffset,tc_axis);
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');


% plot data
if nargin>3 && strcmpi(plottype,'scatter'),
    q1 = scatter(x1, y1, mksz1, markercol1, 'filled', 'markeredgecolor', markeredgecol1);
    %quiver_tri(x1, y1, (x2-x1), (y2-y1), headsize, headangle, linwidth, lincolor)
    q2 = scatter(x2, y2, mksz2, markercol2, 'filled', 'markeredgecolor', markeredgecol2);
else,
    q1 = line(x1, y1, varargin{:});
end
if nargout > 0
    handles1 = q1;
    handles2 = q2;
end
if ~hold_state
    set(gca,'dataaspectratio',[1 1 1]), axis off; set(cax,'NextPlot',next);
end

end

function [fA, fB, fC, fD, fE, fF] = fractions_quiver(A, B, C, D, E, F) % Pani: 2021-01-08, make coordinates in sigmoid scale
Total1 = (A+B+C);
fA = A./Total1;
fB = B./Total1;
fC = 1-(fA+fB);

Total2 = (D+E+F);
fD = D./Total2;
fE = E./Total2;
fF = 1-(fD+fE);
% TERNAXES create ternary axis
%   HOLD_STATE = TERNAXES(MAJORS) creates a ternary axis system using the system
%   defaults and with MAJORS major tickmarks.

% Author: Carl Sandrock 20050211

% To Do

% Modifications

% Modifiers
% (CS) Carl Sandrock
end

function [fA, fB, fC] = fractions(A, B, C) % Pani: 2021-01-08, make coordinates in sigmoid scale
Total1 = (A+B+C);
fA = A./Total1;
fB = B./Total1;
fC = 1-(fA+fB);

% TERNAXES create ternary axis
%   HOLD_STATE = TERNAXES(MAJORS) creates a ternary axis system using the system
%   defaults and with MAJORS major tickmarks.

% Author: Carl Sandrock 20050211

% To Do

% Modifications

% Modifiers
% (CS) Carl Sandrock
end


% Pani (2021-01-09): Include font size and x and y offsets in the inputs
function [hold_state, cax, next] = ternaxes(majors,fontsiz,xoffset,yoffset,tc_axis)

%TODO: Get a better way of offsetting the labels
% % % xoffset = 0.06;
% % % yoffset = 0.02;

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
% % % tc = get(cax,'xcolor'); % Pani: 2021-01-08, make grid lines gray
tc = get(cax,'xcolor'); % Pani: 2021-01-08
% % % tc_axis = [0.6 0.6 0.6]; % Pani: 2021-01-08, make axis color an input (2021-05-15)
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
% % % fSize   = get(cax, 'DefaultTextFontSize'); % Pani (2021-01-09):
% change font size
fSize   = fontsiz;
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');

set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   get(cax, 'FontSize'), ...
    'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
    'DefaultTextUnits','data')

% only do grids if hold is off
if ~hold_state
	%plot axis lines
	hold on;
	plot ([0 1 0.5 0],[0 0 sin(1/3*pi) 0], 'color', tc_axis, 'linewidth',1,...
                   'handlevisibility','off');
	set(gca, 'visible', 'off');

    % plot background if necessary
    if ~isstr(get(cax,'color')),
       patch('xdata', [0 1 0.5 0], 'ydata', [0 0 sin(1/3*pi) 0], ...
             'edgecolor',tc,'facecolor',get(gca,'color'),...
             'handlevisibility','off');
    end
    
	% Generate labels
	majorticks0 = linspace(0, 1, majors + 1);
	majorticks = majorticks0(1:end-1);
% % % 	labels = num2str(majorticks'*100); % Pani: 2021-01-08, make 0 to 1
% labels
    labels = num2str([0, 0.1, 0.5, 0.9]');
	
    zerocomp = zeros(size(majorticks)); % represents zero composition
    
	% Plot right labels (no c - only b a)
	[lxc, lyc] = terncoords(1-majorticks, majorticks, zerocomp);
% % % 	text(lxc, lyc, [repmat('  ', length(labels), 1) labels]);
	text(lxc, lyc, [repmat('  ', length(labels), 1) labels],'FontSize',fontsiz); % Pani (2021-01-09): Font size to tick labels
	
	% Plot bottom labels (no b - only a c)
	[lxb, lyb] = terncoords(majorticks, zerocomp, 1-majorticks); % fB = 1-fA
% % % 	text(lxb, lyb, labels, 'VerticalAlignment', 'Top');
	text(lxb, lyb, labels, 'VerticalAlignment', 'Top','FontSize',fontsiz);
	
	% Plot left labels (no a, only c b)
	[lxa, lya] = terncoords(zerocomp, 1-majorticks, majorticks);
% % % 	text(lxa-xoffset, lya, labels);
	text(lxa-xoffset, lya+yoffset, labels,'FontSize',fontsiz);
	
	nlabels = length(labels)-1;
    % Pani: 2021-01-08, make grid lines with the sigmoid scaling
% % % 	for i = 1:nlabels
% % %         plot([lxa(i+1) lxb(nlabels - i + 2)], [lya(i+1) lyb(nlabels - i + 2)], ls, 'color', 'r', 'linewidth',1,...
% % %            'handlevisibility','off');
% % %         plot([lxb(i+1) lxc(nlabels - i + 2)], [lyb(i+1) lyc(nlabels - i + 2)], ls, 'color', 'b', 'linewidth',1,...
% % %            'handlevisibility','off');
% % %         plot([lxc(i+1) lxa(nlabels - i + 2)], [lyc(i+1) lya(nlabels - i + 2)], ls, 'color', 'g', 'linewidth',1,...
% % %            'handlevisibility','off');
% % % 	end;
    
    labelsnum = [0.1, 0.5, 0.9];
	for i = 1:nlabels
        axistest = labelsnum(i);
        
        % Axis 1
        rAxis1_map0 = axistest+0*linspace(0,1,100)';
        rAxis2_map0 = linspace(0,(1-axistest),100)';
        rAxis3_map0 = 1 - rAxis1_map0 - rAxis2_map0;
        
        [rAxis1_map,rAxis2_map,rAxis3_map] = DeltaLogMaker(rAxis1_map0,rAxis2_map0,rAxis3_map0);
        
        [frAxis1_map, frAxis2_map, frAxis3_map] = fractions(rAxis1_map, rAxis2_map, rAxis3_map);
        
        [x_map1, y_map1] = terncoords(frAxis1_map, frAxis2_map, frAxis3_map);
        
        plot(x_map1, y_map1, ls, 'color', tc_axis, 'linewidth',1,...
           'handlevisibility','off');
       
        % Axis 2
        rAxis2_map0 = axistest+0*linspace(0,1,100)';
        rAxis3_map0 = linspace(0,(1-axistest),100)';
        rAxis1_map0 = 1 - rAxis2_map0 - rAxis3_map0;
        
        [rAxis1_map,rAxis2_map,rAxis3_map] = DeltaLogMaker(rAxis1_map0,rAxis2_map0,rAxis3_map0);
        
        [frAxis1_map, frAxis2_map, frAxis3_map] = fractions(rAxis1_map, rAxis2_map, rAxis3_map);
        
        [x_map2, y_map2] = terncoords(frAxis1_map, frAxis2_map, frAxis3_map);
        
        plot(x_map2, y_map2, ls, 'color', tc_axis, 'linewidth',1,...
           'handlevisibility','off');

        % Axis 3
        rAxis3_map0 = axistest+0*linspace(0,1,100)';
        rAxis1_map0 = linspace(0,(1-axistest),100)';
        rAxis2_map0 = 1 - rAxis3_map0 - rAxis1_map0;
        
        [rAxis1_map,rAxis2_map,rAxis3_map] = DeltaLogMaker(rAxis1_map0,rAxis2_map0,rAxis3_map0);
        
        [frAxis1_map, frAxis2_map, frAxis3_map] = fractions(rAxis1_map, rAxis2_map, rAxis3_map);
        
        [x_map3, y_map3] = terncoords(frAxis1_map, frAxis2_map, frAxis3_map);
        
        plot(x_map3, y_map3, ls, 'color', tc_axis, 'linewidth',1,...
           'handlevisibility','off');
       
	end;

end;

% Reset defaults
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits', fUnits );

end
% TERNCOORDS calculate rectangular coordinates of fractions on a ternary plot
%   [X, Y] = TERNCOORDS(FA, FB) returns the rectangular X and Y coordinates
%   for the point with a fraction defined by FA and FB.  It is assumed that
%   FA and FB are sensible fractions.
%
%   [X, Y] = TERNCOORDS(FA, FB, FC) returns the same.  FC is assumed to be
%   the remainder when subtracting FA and FB from 1.

%       b
%      / \
%     /   \
%    c --- a 

% Author: Carl Sandrock 20050211

% Modifications

% Modifiers
function [x, y] = terncoords(fA, fB, fC)
if nargin < 3
    fC = 1 - (fA + fB);
end

y = fB*sin(deg2rad(60));
x = fA + y*cot(deg2rad(60));
end

function [x1, y1, x2, y2] = terncoords_quiver(fA, fB, fC, fD, fE, fF)
if nargin < 3
    fC = 1 - (fA + fB);
end

y1 = fB*sin(deg2rad(60));
x1 = fA + y1*cot(deg2rad(60));

y2 = fE*sin(deg2rad(60));
x2 = fD + y2*cot(deg2rad(60));

end
