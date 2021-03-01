function [delta_area,delta_width,delta_length] = get_delta_area(Discharge_prist,QRiver_prist,delta_name,shelf_depth,MouthLon,MouthLat,BasinID2)

[data,doug_id] = xlsread('Edmondsetal2020_NatCom_suppdata.xlsx','A3:K2177');
doug_id = str2double(doug_id(2:end,1));
doug_BasinID2 = zeros(size(doug_id));

for ii=1:length(doug_id)
doug_delta_area(ii) = areaint([data(ii,1) data(ii,5) data(ii,3) data(ii,7)] ,[data(ii,2) data(ii,6) data(ii,4) data(ii,8)],referenceSphere('earth'));
end

delta_area_proxy = 1.07.*Discharge_prist.^1.1.*QRiver_prist.^0.45./max(100,-shelf_depth).*1e6; %delta area from syvitski2009

%xx = sort(delta_area,'descend'); yy = sort(delta_area_proxy,'descend');
%plot(xx(1:2000),yy(1:2000)), grid on, axis equal
%set(gca,'XScale','log','YScale','log')

doug_delta_width = zeros(size(doug_id));
for ii=1:length(doug_id)
doug_delta_width(ii) = distance(data(ii,5),data(ii,6),data(ii,7),data(ii,8),referenceSphere('earth'));
end
doug_delta_length = zeros(size(doug_id));
for ii=1:length(doug_id)
doug_delta_length(ii) = distance(data(ii,1),data(ii,2),data(ii,3),data(ii,4),referenceSphere('earth'));
end

data(:,2) = mod(data(:,2)-1,360)+1;



%do the largest delta first, then work backwards
[~,delta_idx] = sort(doug_delta_area,'descend');
idx2 = zeros(size(doug_delta_area));


%define matching parameter
delta_coor = MouthLon+1i.*MouthLat;
dougdelta_coor = data(:,2)+1i.*data(:,1);
%do separate function for large and small delta
blub = abs(rot90(dougdelta_coor)-delta_coor) + abs(log10(doug_delta_area)-log10(delta_area_proxy)) + (abs(rot90(dougdelta_coor)-delta_coor)>2).*10;
%blub(:,delta_idx(1:200)) = abs(rot90(dougdelta_coor(delta_idx(1:200)))-delta_coor) + abs(log10(doug_delta_area(delta_idx(1:200)))-log10(delta_area_proxy)) + (abs(rot90(dougdelta_coor(delta_idx(1:200)))-delta_coor)>4).*10;
%blub(:,delta_idx(201:end)) = abs(rot90(dougdelta_coor(delta_idx(201:end)))-delta_coor) + abs(log10(doug_delta_area(delta_idx(201:end)))-log10(delta_area_proxy)) + (abs(rot90(dougdelta_coor(delta_idx(201:end)))-delta_coor)>1).*10;

%histogram(min(blub,[],2))
%


for ii=delta_idx
    [d,idx2(ii)] = min(blub(:,ii));
    %d
    %don't match if matching parameter exceeds 50
    if d>5, idx2(ii) = 0; continue, end
    
    %set blub to infinite to make sure a POC flux is not matched twice
    blub(idx2(ii),:) = inf;
    doug_BasinID2(ii) = BasinID2(idx2(ii));
    
    
end

%{
plot(delta_coor,'o'), hold on, plot(dougdelta_coor,'or')
hold on
plot([data(idx2>0,2) MouthLon(idx2(idx2>0))]',[data(idx2>0,1) MouthLat(idx2(idx2>0))]')

sum(QRiver_prist(idx2(idx2>0)))./sum(QRiver_prist)

histogram(doug_delta_area(idx2>0)-delta_area_proxy(idx2(idx2>0)))

scatter(doug_delta_area(idx2>0),delta_area_proxy(idx2(idx2>0)))
set(gca,'XScale','log','YScale','log')
%}

delta_area = nan(size(delta_name));
delta_width = nan(size(delta_name));
delta_length= nan(size(delta_name));

delta_area(idx2(idx2>0)) = doug_delta_area(idx2>0);
delta_width(idx2(idx2>0)) = doug_delta_width(idx2>0);
delta_length(idx2(idx2>0)) = doug_delta_length(idx2>0);

