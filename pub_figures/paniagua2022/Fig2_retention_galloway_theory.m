hs = 0.41; %m
a0=0.57; %m
ret = 0.64;
Qriver = 10.*ret; %kg/s
discharge = max(1,500); %m3/s

for ii=1:length(ret),
    [Qwave(ii),Qtide(ii),~,shoreline_angle(ii),w_upstream(ii),w_downstream(ii)] = galloway_predictor(hs,a0,Qriver(ii),discharge)
end


s = Qriver+Qwave+Qtide;
[Rriver,Rwave,Rtide] = deal(Qriver./s,Qwave./s,Qtide./s);
y0 = Rriver*sqrt(3)/2;
x0 = Rtide + y0*sqrt(3)/3;

plot(x0,y0,'o');
hold on,

plot([0,0.5,1,0],[0,sqrt(3)/2,0,0])


%example point:
%river, wave, tide
ret = 1;

obs = [10;13;12];
pred = [10;15;25];

pred = [10.*ret;15.*ones(size(ret));ret*25];
Robs = pred(1,:)/pred(2,:)
Tobs = pred(3,:)/pred(1,:)
obs = obs./sum(obs);
pred = pred./sum(pred)
scatter(obs(3)+0.5*obs(1),obs(1)*sqrt(3)/2);
hold on
scatter(pred(3,:)+0.5*pred(1,:),pred(1,:)*sqrt(3)/2);
plot([0,0.5,1,0],[0,sqrt(3)/2,0,0])

d = sqrt(sum((pred-obs).^2))


sqrt(sum((pred(1:2)-obs(1:2)).^2))



sqrt(((((obs(1)*sqrt(3)/2))-(pred(1)*sqrt(3)/2))).^2+((obs(3)+0.5*obs(1))-(pred(3)+0.5*pred(1))).^2)