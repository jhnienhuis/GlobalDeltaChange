v1 = load('GlobalDeltaData_v1.mat','QRiver_dist','QRiver_prist','ee');
flux_ratio = v1.QRiver_dist./v1.QRiver_prist;
land = v1.ee.net_aqua;

bx = 0.5;

for ii=1:1000,
    idx = randsample(length(flux_ratio),9762);
    x = land(idx);
    b = flux_ratio(idx);
    
    ca(ii,1) = mean(x(b<bx));
    ca(ii,2) = mean(x(b>bx & b<(2-bx)));
    ca(ii,3) = mean(x(b>(2-bx)));

end
hold on,

violin(ca,'mc','','medc','');
legend off
scatter(1:3,ca,40,'xr')
xticklabels({'<-50%','-50% .. +50%','>+50%'})
xlabel('Human-induced Sediment Flux Change')
xticks([1 2 3]), xlim([0.5 3.5]), box on,
ylabel('Delta Area Change (km^2/yr)')
ylim([-0.02 0.02]), plot([0.5 3.5],[0 0],':k')
saveas(gcf,'Figure1.svg')
