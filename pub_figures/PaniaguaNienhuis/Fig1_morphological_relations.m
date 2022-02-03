function Fig1_morphological_relations
%% first only R
subplot(2,2,1)
hs = 0:0.02:2; %m
a0=0; %m
Qriver = 10; %kg/s
discharge = 1000; %m3/s
for ii=1:length(hs),
    [Qwave(ii),Qtide(ii),~,shoreline_angle(ii),w_upstream(ii),w_downstream(ii)] = galloway_predictor(hs(ii),a0,Qriver,discharge);
end
shoreline_angle(isnan(shoreline_angle)) = 50;
plot(Qriver./Qwave,shoreline_angle,'-o'), set(gca,'XScale','log')
yticks([0 10 20 30 40 50])
yticklabels({'0','10','20','30','40','undefined'})
ylabel('Shoreline Angle')
xlabel('R (QRiver/QWave)')
xlim([1e-2 1e2])


%plot(1./(1+(Qriver./Qwave)),shoreline_angle,'-o'),

%plot(1./(1+(Qriver./Qwave)),shoreline_angle,'-o'),
%hold on
%plot((Qriver./Qwave),(Qriver./Qwave),'-o'), set(gca,'XScale','log')

%% and only T
subplot(2,2,2)
hs = 0; %m
a0=0:0.01:3; %m
Qriver = 10; %kg/s
discharge = 1000; %m3/s
for ii=1:length(a0),
    [Qwave,Qtide(ii),~,shoreline_angle(ii),w_upstream(ii),w_downstream(ii)] = galloway_predictor(hs,a0(ii),Qriver,discharge);
end
plot(Qtide./Qriver,w_downstream./w_upstream,'-o'), set(gca,'XScale','log')
ylabel('Channel widening')
xlabel('T (QTide/QRiver)')
xlim([1e-2 1e2])

%plot(1./(1+(Qtide./Qriver)),w_downstream./w_upstream,'-o')

%% now both in triangle
hs = logspace(-4,1,50); %m
a0=logspace(-4,1,50); %m
Qriver = 10; %kg/s
discharge = 1000; %m3/s

[hs_mesh,a0_mesh] = meshgrid(hs,a0);

for ii=1:length(hs_mesh(:)),
    [Qwave(ii),Qtide(ii),~,shoreline_angle(ii),w_upstream(ii),w_downstream(ii)] = galloway_predictor(hs_mesh(ii),a0_mesh(ii),Qriver,discharge);
end



%[Qriver_log,Qwave_log,Qtide_log] = DeltaLogMaker(Qriver.*ones(size(Qwave)),Qwave,Qtide);
s = Qriver+Qwave+Qtide;
[Qriver,Qwave,Qtide] = deal(Qriver./s,Qwave./s,Qtide./s);
y0 = Qriver*sin(deg2rad(60));
x0 = Qtide + y0*cot(deg2rad(60));

xn = 0:0.005:1;
yn = 0:0.005:1;

%plot shoreline angle
subplot(2,2,3)
[fn,xgrid,ygrid] = gridfit(x0,y0,shoreline_angle,xn,yn); %'interp','bilinear',
out = (yn>(xn.*2*sqrt(3/4))' | yn>(2*sqrt(3/4).*(1-xn)'))';
fn(out) = nan;
fn(fn>50) = nan;
contourf(xgrid,ygrid,fn,[0:5:50],'ShowText','on')
colormap(flipud(cbrewer('div', 'RdYlBu', 64)))
caxis([0,50])
set(gca,'Layer','top')
%set([findall(gcf,'String','  25'); findall(gcf,'String','25')],'String','0.1')
%set([findall(gcf,'String','  50'); findall(gcf,'String','50')],'String','0.5')
%set([findall(gcf,'String','  75'); findall(gcf,'String','75')],'String','0.9')
%hold on
%scatter(x0,y0,'.k')
h = colorbar('East');
axis equal
ylabel(h,'shoreline angle')

%plot river widening
subplot(2,2,4)
[fn,xgrid,ygrid] = gridfit(x0,y0,w_downstream./w_upstream,xn,yn); %'interp','bilinear',
out = (yn>(xn.*2*sqrt(3/4))' | yn>(2*sqrt(3/4).*(1-xn)'))';
fn(out) = nan;
contourf(xgrid,ygrid,fn,[1:0.1:2],'ShowText','on')
colormap(flipud(cbrewer('div', 'RdYlBu', 64)))
caxis([1,2])
set(gca,'Layer','top')
%set([findall(gcf,'String','  25'); findall(gcf,'String','25')],'String','0.1')
%set([findall(gcf,'String','  50'); findall(gcf,'String','50')],'String','0.5')
%set([findall(gcf,'String','  75'); findall(gcf,'String','75')],'String','0.9')
%hold on
%scatter(x0,y0,'.k')
h = colorbar('East');
axis equal
ylabel(h,'channel widening')

set(gca, 'FontSize', 8,'FontName','Helvetica')
saveas(gcf,'Fig1_morphological_relations.svg')