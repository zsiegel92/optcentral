%% Global Parameters
% nz = input('How many product skus are there? Include raw materials, components, and finished goods.\n')
% nt = 10; %Will be accounted for by forecast.m
% nloc = input('How many locations are there? Include factories, warehouses, and retail locations.\n')

autofill = [11,20,5];%Default [nz,nt,nloc]
%NOTE nt can be MAXIMUM size(extraTS,2)+nt0
c = {'nz (number skus)','nt (see forecast.m)','nloc (# Locations)'};
%AUTOFILLING
for i = 1:size(c,2)
    c{i} = [c{i},',',num2str(autofill(i))];
end


filename_global = '~/Desktop/Theo_K./UI/UI_csvs/global.csv';
filename = filename_global;
fid = fopen(filename,'w');
fprintf(fid, '%s\n',c{1,:});
fclose(fid);

unix(['open ', filename]);
input('Please fill out "constants.csv", save the file, and press return.');
constants = (csvread(filename,0,1));
nz=constants(1);
nt= constants(2);
nloc=constants(3);


 
%% Bom, FG (finished good), RM (raw material)
%(Probably just import)
if (nz ==11)
    [Bom,Gbom]=BOMCREATOR();%Creates a 11-product Bill of Materials Matrix
                            %Gbom is the BOM graph, p is the graphplot of GBom
p = BOMVisualizer(Gbom);
else
    Bom = triu(ones(nz,nz),4);
end
 %Bom(i,I) = number of units of item I required to produce 1 unit of item i DIRECTLY (aka "B").
%This dummy BOM means that the first four products aren't required to
%produce anything else. Also, the final four products are fully raw
%materials.
%% One-off Constants

autofill = [500,.1,10,15,.1];%Default [c_truck,perMileTruckCost,hourlyTruckCost,fixedTruckCost,mu_w]

c = {'c_truck','perMileTruckCost','hourlyTruckCost','fixedTruckCost','mu_w'};
%c = {'c_truck (capacity in "units"),100','perMileTruckCost ($/mile),.5','hourlyTruckCost ($/hr),25','fixedTruckCost (fixed $/truck),50','mu_w (demand backlog decay rate),.1'} %Autofilled
%AUTOFILLING
for i = 1:size(c,2)
    c{i} = [c{i},',',num2str(autofill(i))];
end





filename_constants = '~/Desktop/Theo_K./UI/UI_csvs/constants.csv';
filename = filename_constants;
fid = fopen(filename,'w');
fprintf(fid, '%s\n',c{1,:});
fclose(fid);


unix(['open ', filename]);
input('Please fill out "constants.csv", save the file, and press return.');
constants = (csvread(filename,0,1));

c_truck =constants(1); %Capacity of one truck
perMileTruckCost = constants(2);%Cost per mile of trucking
hourlyTruckCost = constants(3);%Cost per hour of trucking
fixedTruckCost = constants(4); %Cost to summon a truck
mu_w = constants(5)*ones(nz,nloc); %Decay rate of inventory
clear constants;

%% Is this a finished good? A raw material? A subassembly? (to make WC and ZC entry easier)

%% s_y (and f_w unless considering elasticity)

autofill = [.1*ones(nz,1),25*ones(nz,1)]; %Amount of storage space required per product, sales price of each product

c = {'Sku', 'Storage Requirement per Unit','Sales Price Point'};
format = '%f\n';
skus = (1:nz)';

%Autofilling
skus = [skus,autofill];%COMMENT THIS TO PREVENT AUTOFILL
for i = 1:size(skus,2)-1
    format = ['%f,',format];
end



filename_sku = '~/Desktop/Theo_K./UI/UI_csvs/sku.csv';
filename = filename_sku;
fid = fopen(filename,'w');
fprintf(fid, '%s,%s, %s\n',c{1,:});
fprintf(fid,format,skus');
fclose(fid);
disp('Please fill out "sku.csv", save the file, and press return.');
unix(['open ', filename]);
input('Please enter the storage requirements and sales price-points for each sku.');
constants = csvread(filename,1,1);

s_y = constants(:,1); %default 1
prices = constants(:,2); %default 25
f_w = repmat(prices,1,nt,nloc);%unused time-, location-varying pricing
clear prices; clear constants;

%% WC
c = {'Sku'};
for i = 1:nloc
    c{1,i+1} = ['Location ', int2str(i)];
end
for i = 1:nz
   c{i+1,1}=int2str(i); 
end
lineform = '%s\n';
for i = 1:nloc
    lineform = ['%s,',lineform];
end

filename_wc = '~/Desktop/Theo_K./UI/UI_csvs/wc.csv';
filename = filename_wc;
fid = fopen(filename,'w');
%fprintf(fid,'%s\n','In the Kth column and ith row: please enter a 1 if product i can be sold at location k.')
fprintf(fid, lineform,c{1,:});
fprintf(fid,'%s\n',c{2:end,:});
fclose(fid);

disp('Please fill out "wc.csv", save the file, and press return.');
unix(['open ', filename]);
input('In the Kth column and ith row: please enter a 1 if product i can be sold at location k, and a (-1) if it can be bought.');
wc = csvread(filename,1,1);

%% ZC
filename_zc = '~/Desktop/Theo_K./UI/UI_csvs/zc.csv';
filename = filename_zc;
fid = fopen(filename,'w');
%fprintf(fid,'%s\n','In the Kth column and ith row: please enter a 1 if product i can be sold at location k.')
fprintf(fid, lineform,c{1,:});
fprintf(fid,'%s\n',c{2:end,:});
fclose(fid);

disp('Please fill out "zc.csv", save the file, and press return.');
unix(['open ', filename]);
input('In the Kth column and ith row: please enter a 1 if product i can be produced at location k.');
zc = csvread(filename,1,1);


%% Transportation and location-based parameters R_t, f_t; c_y,f_y,coordinates,locNames

filename_trans = '~/Desktop/Theo_K./UI/UI_csvs/trans.csv';
filename = filename_trans;
fid = fopen(filename,'w');

autofillCities = {'Rochester','New York','Boston','Philadelphia','Pittsburg','Trenton','Durham','Long Beach','Torrance','Omaha'};
autofillparams = [10000,.1];%Capcity to store, cost to store per unit
autofill = {};
for i = 1:nloc
    autofill{i,1} = i;
    autofill{i,2}=autofillCities{i};
    autofill{i,3} =autofillparams(1);
    autofill{i,4}=autofillparams(2);
end

% Creating Address Template
c = {'Location ID', 'City','Storage Capacity in "units"','Cost to Store per "unit"'};
fprintf(fid, '%s,%s,%s,%s\n',c{1,:});
for i = 1:nloc
    fprintf(fid,'%f,%s,%f,%f\n',autofill{i,:}); %Print sheet with autofill
end
% fprintf(fid,'%f\n',1:nloc); %Print blank sheet with location IDs
fclose(fid);

% Reading Addresses, saving them in locString
disp('Please fill out "trans.csv", save the file, and press return.');
unix(['open ', filename]);
input('Please enter the city name and information for each location according to its "Location ID".');
% locations = csvread(filename,1,1);
fid = fopen(filename);
fgets(fid);
locations = textscan(fid,'%f %s %f %f','Delimiter',',');
fclose(fid);
locString = strjoin(locations{2},'|');
locString = strrep(locString,' ','+');
c_y=locations{3};%Capacity to hold inventory (in "units")
f_y=locations{4}; %Cost to hold inventory (per "unit")
f_y = repmat(f_y',nt,1);

locNames = locations{2};

% Querying Google Maps
%Google Maps API Credential: AIzaSyDDKHV64CvECHE8wf_KXpLiWr2fM0XAnrU
APIkey='AIzaSyDDKHV64CvECHE8wf_KXpLiWr2fM0XAnrU';
pre_url = 'https://maps.googleapis.com/maps/api/distancematrix/json?';
origins = ['origins=',locString];
dests = ['destinations=',locString];

url = [pre_url,origins,'&',dests,'&key=',APIkey];
%url = 'https://maps.googleapis.com/maps/api/distancematrix/json?origins=nottingham&destinations=london|manchester|liverpool&key=AIzaSyDDKHV64CvECHE8wf_KXpLiWr2fM0XAnrU';
result = webread(url);


% Copying API Result into shipping times R_t and costs f_t. NOTE: (k,L)'th entry denotes k->L route.

% perMileTruckCost = .2;
% hourlyTruckCost = .2;
% fixedTruckCost = 1.5;

R_t = zeros(nloc,nloc);
d_t = zeros(nloc,nloc);
for k = 1:nloc
    for L = 1:nloc
        d_t(k,L)=result.rows(k).elements(L).distance(1).value; %Distance of route in meters!
        R_t(k,L) = result.rows(k).elements(L).duration(1).value; %Duration of route in SECONDS!
    end
end
d_t = d_t*(1/1609.34); %Convert meters to miles.
R_t = R_t*(1/3600)*(1/24); %DAILY

f_t = fixedTruckCost*ones(nloc,nloc)+hourlyTruckCost*R_t + perMileTruckCost*d_t; %Set f_t based on pre-rounded trucking times.
R_t = ceil(R_t); %Hourly, rounded up to nearest hour - must convert to "time-steps"
R_t = R_t + eye(size(R_t)); %Takes one time-step to ship from location to self! Thus, trucks can be used for storage.

%Coordinates
%Note that 'locations{2}' contains the location-names as inputted by user
APIkey2='AIzaSyDVZPQrCVHtayFWu2LIWid6126hujMG5kQ';
pre_url2 = 'https://maps.googleapis.com/maps/api/geocode/json?';
coordinates = cell(1,size(locNames,1));
for i = 1:size(locNames,1)
% address = ['address=',locString];
address = ['address=',strrep(locNames{i},' ','+')];
url2 = [pre_url2,address,'&key=',APIkey2];
%url = 'https://maps.googleapis.com/maps/api/distancematrix/json?origins=nottingham&destinations=london|manchester|liverpool&key=AIzaSyDDKHV64CvECHE8wf_KXpLiWr2fM0XAnrU';
webreturn = webread(url2);
webreturn = webreturn.results.geometry;
webreturn = webreturn.location;
coordinates{i}= webreturn;
end
coordinates = cell2mat(cellfun(@(x) [x.lat;x.lng],coordinates,'UniformOutput',0)); %coordinates(1,k) is the latitude of the k'th location; coordinates(2,k) is the longitude.

%% All Production Parameters: f_z, f_z0, c_z, R_z
%%%TEMPLATE
c = {'', 'Production Cost', 'FIXED Prod. Cost', 'Prod. Capacity', 'Prod. Time';'sku','f_z','f_z0','c_z','r_z'};
autofill = [.1,5,1000,2];% = [f_z,f_z0,c_z,R_z]
filename_production = '~/Desktop/Theo_K./UI/UI_csvs/sku_vs_loc.csv';
filename = filename_production;
fid = fopen(filename,'w');
   fprintf(fid,'%s,%s,%s,%s,%s\n',c{1,:});
   for i = 1:nloc
       if(nnz(zc(:,i))>0)
           facnum = {'Location and ID:', locations{2}{i},i};
           fprintf(fid,'%s,%s,%f\n',facnum{1,:});
           fprintf(fid,'%s,%s,%s,%s,%s\n',c{2,:});
           %a = find(zc(:,i));
           %fprintf(fid,'%f\n',find(zc(:,i)));%USE THIS LINE IF NOT
           %AUTOFILLING
           fprintf(fid,'%f,%f,%f,%f,%f\n',[find(zc(:,i)),repmat(autofill,nnz(zc(:,i)),1)]');
           fprintf(fid,'\n');
       end
   end
   fclose(fid);
unix(['open ', filename]);
input('Fill in production parameters in "sku_vs_loc.csv" and press enter.');

%%%READ-IN
%f_z0 = ones(nz,nt,nloc); %f_z0(i,j,k)=fixed cost of starting production of product i at time-step j at location k
%f_z = ones(nz,nt,nloc); %f_z(i,j,k) = PRICE OF PRODUCTION per unit of product i, at location k, at time-step j

f_z0 = zeros(nz,nloc); %f_z0(i,j,k)=fixed cost of starting production of product i at time-step j at location k
f_z = zeros(nz,nloc); %f_z(i,j,k) = PRICE OF PRODUCTION per unit of product i, at location k, at time-step j
R_z = zeros(nz,nloc); %R_z(i,k) = NUMBER OF TIME-STEPS REQUIRED TO PRODUCE PRODUCT i AT LOCATION k
c_z = zeros(nz,nloc);%c_z(i,k) = number of units of product i produced in one production cycle at location k

linecount = 3;
for i = 1:nloc
    if (nnz(zc(:,i))>0)
        k = csvread(filename,linecount-2,2,[linecount-2,2,linecount-2,2]);%Location ID
        info = csvread(filename,linecount,0, [linecount,0,(nnz(zc(:,k))+linecount-1),4]);
        for j = 1:size(info,1)
            f_z(info(j,1),k)=info(j,2);
            f_z0(info(j,1),k)=info(j,3);
            c_z(info(j,1),k)=info(j,4);
            R_z(info(j,1),k)=info(j,5);
        end
        linecount = linecount + nnz(zc(:,i))+3;
    end
end


f_z = repmat(f_z,1,1,nt);%Add time-dimension to cost matrix
f_z = permute(f_z,[1,3,2]);%Arrance matrix to be (nz,nt,nloc);
f_z0= repmat(f_z0,1,1,nt);%Add time-dimension to cost matrix
f_z0=permute(f_z0,[1,3,2]);%Arrance matrix to be (nz,nt,nloc);



%% TC

filename_tc = '~/Desktop/Theo_K./UI/UI_csvs/tc.csv';
filename = filename_tc;
fid = fopen(filename,'w');

%fprintf(fid,'%s\n','In the Kth column and ith row: please enter a 1 if product i can be sold at location k.')
fprintf(fid, '%s\n',[',',strjoin(locations{2},',')]);
c = locations{2}';
fprintf(fid,'%s\n',c{1,:});
clear c;
fclose(fid);

disp('Please fill out "tc.csv", save the file, and press return.');
unix(['open ', filename]);
input('In the kth column and Lth row: please enter a 1 if the route from location k to location L is NEVER travelled.');
tc = (csvread(filename,1,1)==0);
for i = 1:size(tc,2)
    tc(i,i)=1;
end


