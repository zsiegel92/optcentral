%% Writing Returned Values to .csv
%Writing Returned Values to .csv
filename = '~/Desktop/Theo_K./Central Intel Code/Supply_Output.csv';
c0 = {'Product (SKU)','Time-Step','Location (1)', 'Location (2)'};
% c = {{'Product (SKU)','Time-Step','Location','Amount Held-Over'};
%     {'Product (SKU)','Time-Step','Location','Amount Produced'};
%     {'Product (SKU)','Time-Step','Location 1','Location 2','Amount Shipped'};
%     {'Product (SKU)','Time-Step','Location','Amount Sold'};
%     {'Product (SKU)','Time-Step','Location','Start Production? (1 or 0)'};
%     {'Time-Step','Location 1','Location 2', 'Number Trucks Sent from 1 to 2'};
%     {'Time-Step','Location','"In-Production"? (1 or 0)'}};%Y,Z,T,W,Z0,T0,Z00
labels = {'Holding Orders';'Production Orders';'Transportation Orders';'Sales Benchmarks';'Production (fixed) Orders';'Truck Orders';'Production (availability) Orders'};
labels2 = {'Amount Held'; 'Amount Produced'; 'Amount Transported'; 'Amount Sold'; 'Start Production? (1 or 0)'; 'Number Trucks Sent'; '"In Production"? (1 or 0)'};
%labels2 = {'Amount Held - Y'; 'Amount Produced - Z'; 'Amount Transported - T'; 'Amount Sold - W'; 'Start Production? (1 or 0) - Z0'; 'Number Trucks Sent - T0'; '"In Production"? (1 or 0 - Z00'};
termHelp = {'%s,%s,%s,%s\n','%s,%s,%s,%s\n','%s,%s,%s,%s,%s\n','%s,%s,%s,%s\n','%s,%s,%s,%s\n','%s,%s,%s,%s\n','%s,%s,%s\n'}; %Y, Z, T, W, Z0, T0, Z00
fid = fopen(filename,'w');

if (size(track_Y,1)>0)
    fprintf(fid,'%s\n',labels{1,:});%Y
    %fprintf(fid,'%s,%s,%s,%s\n',c{1}{:});%Y
    c = {c0{viewOrderY-1},labels2{1}};%Y
    fprintf(fid,termHelp{1},c{:});%Y
    fprintf(fid,'%f,%f,%f,%f\n',track_Y');%Y
end
if (size(track_T,1)>0)
    fprintf(fid,'%s\n',labels{3,:});%T
    %fprintf(fid,'%s,%s,%s,%s,%s\n',c{3}{:});%T
    c = {c0{viewOrderT-1},labels2{3}};%T
    fprintf(fid,termHelp{3},c{:});%T
    fprintf(fid,'%f,%f,%f,%f,%f\n',track_T');%T
end
if (size(track_T0,1)>0)
    fprintf(fid,'%s\n',labels{6,:});%T0
    %fprintf(fid,'%s,%s,%s,%s\n',c{6}{:});%T0
    c = {c0{viewOrderT0-1},labels2{6}};%T0
    fprintf(fid,termHelp{6},c{:});%T0
    fprintf(fid,'%f,%f,%f,%f\n',track_T0');%T0
end
if (size(track_W,1)>0)
    fprintf(fid,'%s\n',labels{4,:});%W
    %fprintf(fid,'%s,%s,%s,%s\n',c{4}{:});%W
    c = {c0{viewOrderW-1},labels2{4}};%W
    fprintf(fid,termHelp{4},c{:});%W
    fprintf(fid,'%f,%f,%f,%f\n',track_W');%W
end
if (size(track_Z,1)>0)
    fprintf(fid,'%s\n',labels{2,:});%Z
    c = {c0{viewOrderZ-1},labels2{2}};%Z
    fprintf(fid,termHelp{2},c{:});%Z
    %fprintf(fid,'%s,%s,%s,%s\n',c{2}{:});%Z
    fprintf(fid,'%f,%f,%f,%f\n',track_Z');%Z
end

if (size(track_Z0,1)>0)
    fprintf(fid,'%s\n',labels{5,:});%Z0
    %fprintf(fid,'%s,%s,%s,%s\n',c{5}{:});%Z0
    c = {c0{viewOrderZ-1},labels2{5}};%Z0
    fprintf(fid,termHelp{5},c{:});%Z0
    fprintf(fid,'%f,%f,%f,%f\n',track_Z0');%Z0
end

if (size(track_Z00,1)>0)
    fprintf(fid,'%s\n',labels{7,:});%Z00
    %fprintf(fid,'%s,%s,%s\n',c{7}{:});%Z00
    c = {c0{viewOrderZ00-1},labels2{7}};%Z00
    fprintf(fid,termHelp{7},c{:});%Z00
    fprintf(fid,'%f,%f,%f\n',track_Z00');%Z00
end
fclose(fid);