function [] = WriteResults(output_raw,output_tracked,stats)
% output_raw.Y = Y;
% output_raw.Z = Z;
% output_raw.T = T;
% output_raw.W = W;
% output_raw.Z0 = Z0;
% output_raw.T0 = T0;
% output_raw.Z00 = Z00;
% output_raw.WM = WM;
% output_raw.WM0=WM0;
% output_tracked.track_Y=track_Y;
% output_tracked.track_Z=track_Z;
% output_tracked.track_T=track_T;
% output_tracked.track_W=track_W;
% output_tracked.track_Z0=track_Z0;
% output_tracked.track_T0=track_T0;
% output_tracked.track_Z00=track_Z00;
% output_tracked.track_WM=track_WM;
% output_tracked.track_WM0=track_WM0;
% 
% stats.salesOf = salesOf;
% stats.unitsSold = unitsSold;
% stats.salesAtLoc=salesAtLoc;
% stats.salesAtT=salesAtT;
% stats.TCosts=TCosts;
% stats.ZCosts=ZCosts;
% stats.Z0Costs=Z0Costs;
% stats.YCosts=YCosts;



%% Carrying over from SystemStats.m

viewOrderY = [3,4,2]; %[j,k,i,y(i,j,k)]
viewOrderZ = [3,4,2]; %[j,k,i,z(i,j,k)]
viewOrderT = [3,4,5,2]; %[j,k,L,i,t(i,j,k,L)]
viewOrderW = [3,4,2]; %[j,k,i,w(i,j,k)]
viewOrderZ0 = [3,4,2]; %[j,k,i,z0(i,j,k)]
viewOrderT0 = [3,4,5]; %[j,k,L,t0(j,k,L)]
viewOrderZ00 = [3,4]; %[j,k,z00(j,k)]
viewOrderWM = [3,4,2,5];%[j,k,i,m,wm(i,j,k)]
viewOrderWM0 = [3,4,2,5];%[j,k,i,m,wm0(i,j,k)]


salesOf = stats.salesOf;
unitsSold = stats.unitsSold;
salesAtLoc=stats.salesAtLoc;
salesAtT=stats.salesAtT;
TCosts=stats.TCosts;
ZCosts=stats.ZCosts;
Z0Costs=stats.Z0Costs;
YCosts=stats.YCosts;

%Writing Returned Values to .csv
filename = '/Users/Zach/Desktop/Theo_K./Central Intel Code/Tiered_Calculation/Supply_Output.csv';
c0 = {'Product (SKU)','Time-Step','Location (1)', 'Location (2)','Tier'};
% c = {{'Product (SKU)','Time-Step','Location','Amount Held-Over'};
%     {'Product (SKU)','Time-Step','Location','Amount Produced'};
%     {'Product (SKU)','Time-Step','Location 1','Location 2','Amount Shipped'};
%     {'Product (SKU)','Time-Step','Location','Amount Sold'};
%     {'Product (SKU)','Time-Step','Location','Start Production? (1 or 0)'};
%     {'Time-Step','Location 1','Location 2', 'Number Trucks Sent from 1 to 2'};
%     {'Time-Step','Location','"In-Production"? (1 or 0)'}};%Y,Z,T,W,Z0,T0,Z00
labels = {'Holding Orders';'Production Orders';'Transportation Orders';'Sales Benchmarks';'Production Site Un-Availability';'Truck Orders';'Production (availability) Orders';'Acquisitions'};
labels2 = {'Amount Held'; 'Amount Produced'; 'Amount Transported'; 'Amount Sold'; 'Start Production? (1 or 0)'; 'Number Trucks Sent'; '"In Production"? (1 or 0)';'Order Size'};
%labels2 = {'Amount Held - Y'; 'Amount Produced - Z'; 'Amount Transported - T'; 'Amount Sold - W'; 'Start Production? (1 or 0) - Z0'; 'Number Trucks Sent - T0'; '"In Production"? (1 or 0 - Z00'};
termHelp = {'%s,%s,%s,%s\n','%s,%s,%s,%s\n','%s,%s,%s,%s,%s\n','%s,%s,%s,%s\n','%s,%s,%s,%s\n','%s,%s,%s,%s\n','%s,%s,%s\n'}; %Y, Z, T, W, Z0, T0, Z00
fid = fopen(filename,'w');

if (size(output_tracked.track_Y,1)>0)
    fprintf(fid,'%s\n',labels{1,:});%Y
    %fprintf(fid,'%s,%s,%s,%s\n',c{1}{:});%Y
    c = {c0{viewOrderY-1},labels2{1}};%Y
    fprintf(fid,termHelp{1},c{:});%Y
    fprintf(fid,'%f,%f,%f,%f\n',output_tracked.track_Y');%Y
end
if (size(output_tracked.track_T,1)>0)
    fprintf(fid,'%s\n',labels{3,:});%T
    %fprintf(fid,'%s,%s,%s,%s,%s\n',c{3}{:});%T
    c = {c0{viewOrderT-1},labels2{3}};%T
    fprintf(fid,termHelp{3},c{:});%T
    fprintf(fid,'%f,%f,%f,%f,%f\n',output_tracked.track_T');%T
end
if (size(output_tracked.track_T0,1)>0)
    fprintf(fid,'%s\n',labels{6,:});%T0
    %fprintf(fid,'%s,%s,%s,%s\n',c{6}{:});%T0
    c = {c0{viewOrderT0-1},labels2{6}};%T0
    fprintf(fid,termHelp{6},c{:});%T0
    fprintf(fid,'%f,%f,%f,%f\n',output_tracked.track_T0');%T0
end
if (size(output_tracked.track_W,1)>0)
    fprintf(fid,'%s\n',labels{4,:});%W
    %fprintf(fid,'%s,%s,%s,%s\n',c{4}{:});%W
    c = {c0{viewOrderW-1},labels2{4}};%W
    fprintf(fid,termHelp{4},c{:});%W
    fprintf(fid,'%f,%f,%f,%f\n',output_tracked.track_W');%W
end
if (size(output_tracked.track_Z,1)>0)
    fprintf(fid,'%s\n',labels{2,:});%Z
    c = {c0{viewOrderZ-1},labels2{2}};%Z
    fprintf(fid,termHelp{2},c{:});%Z
    %fprintf(fid,'%s,%s,%s,%s\n',c{2}{:});%Z
    fprintf(fid,'%f,%f,%f,%f\n',output_tracked.track_Z');%Z
end

if (size(output_tracked.track_Z0,1)>0)
    fprintf(fid,'%s\n',labels{5,:});%Z0
    %fprintf(fid,'%s,%s,%s,%s\n',c{5}{:});%Z0
    c = {c0{viewOrderZ-1},labels2{5}};%Z0
    fprintf(fid,termHelp{5},c{:});%Z0
    fprintf(fid,'%f,%f,%f,%f\n',output_tracked.track_Z0');%Z0
end

if (size(output_tracked.track_Z00,1)>0)
    fprintf(fid,'%s\n',labels{7,:});%Z00
    %fprintf(fid,'%s,%s,%s\n',c{7}{:});%Z00
    c = {c0{viewOrderZ00-1},labels2{7}};%Z00
    fprintf(fid,termHelp{7},c{:});%Z00
    fprintf(fid,'%f,%f,%f\n',output_tracked.track_Z00');%Z00
end

if (size(output_tracked.track_WM,1)>0)
    fprintf(fid,'%s\n',labels{8,:});%WM
    c = {c0{[2,3,1,5]},labels2{8}};%WM
    fprintf(fid,'%s,%s,%s,%s,%s\n',c{:});%WM
    fprintf(fid,'%f,%f,%f,%f,%f\n',output_tracked.track_WM');%WM
end

fprintf(fid,'%s\n',['Final Holding Costs: ,', num2str(YCosts)]);
fprintf(fid,'%s\n',['Final Production Costs: ,', num2str(ZCosts + Z0Costs)]);
fprintf(fid,'%s\n',['Final Transportation Costs: ,', num2str(TCosts)]);
%fprintf(fid,'%s\n',['Final Profit: ,',num2str(-qOb*x)]);



%Now would be a good time to print all System Statistics

% salesOf = zeros(1,nz); %Revenues earned from each sku
% unitsSold = zeros(1,nz); %Number of units sold of each SKU
% salesAtLoc = zeros(1,nloc); %Revenues at each location
% salesAtT = zeros(1,nt); %Revenues from each time-step
% TCosts = 0;
% ZCosts = 0;
% Z0Costs = 0;
% YCosts = 0;


fclose(fid);

end