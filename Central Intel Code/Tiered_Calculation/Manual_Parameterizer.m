function params = Manual_Parameterizer()

%% Global Parameters

% nz = input('How many product skus are there? Include raw materials, components, and finished goods.\n')
% nt = 10; %Will be accounted for by forecast.m
% nloc = input('How many locations are there? Include factories, warehouses, and retail locations.\n')
nz = 11;
nt = 15;
nloc = 2;
time_step_hours = 4; %one time-step is half a day
nTiers = 3;
nt0 = 1;


%% One-off Constants

c_truck =500; %Capacity of one truck
perMileTruckCost = .1;%Cost per mile of trucking
hourlyTruckCost = 10;%Cost per hour of trucking
fixedTruckCost = 10; %Cost to summon a truck
mu_w = 0*ones(nz,nloc); %Decay rate of inventory



%% Bom, FG (finished good), RM (raw material)
%(Probably just import)
if (nz ==11)
    [Bom,Gbom,RBom,RGBom]=BOMCREATOR();%Creates a 11-product Bill of Materials Matrix
    %Gbom is the BOM graph, p is the graphplot of GBom
    %p = BOMVisualizer(Gbom);
else
    Bom = triu(ones(nz,nz),4);
end
%Bom(i,I) = number of units of item I required to produce 1 unit of item i DIRECTLY (aka "B").
%This dummy BOM means that the first four products aren't required to
%produce anything else. Also, the final four products are fully raw
%materials.
SKU_names= {'Wheat Berries (lb)';'Beef (lb)';'Egg (ea)';'Salt (g)';'Wheat Flour (lb)';'Ground Beef (lb)';'Roll/Bun (ea)';'Burger (ea) FG';'Egg Sandwich (ea) FG';'Mayonnaise (lb)';'Vegetable Oil (lb)'};
RM = [1;1;1;1;0;0;0;0;0;0;1]; %Raw Materials (purchase-able) %[1,2,3,4,11] %Wheat Berries, Beef, Egg, Salt,Vegetable Oil
FG = [0;0;0;0;0;0;0;1;1;0;0]; %[7,8,9] %Finished Goods (sell-able) %Roll/Bun, Egg Sandwich, Burger
SA = [0,0,0,0,1,1,1,0,0,1,0]; %Sub-Assemblies (called for in production) %[5,6,7,10]


%% f_w, c_w, f_ws, and s_y

f_w1 = zeros(1,nz);%Price to BUY product i raw (tier 1 - retail)
f_w2 = zeros(1,nz);%Price to BUY product i raw (tier 2 - wholesale)
f_w3 = zeros(1,nz);%Price to BUY product i raw (tier 3 - bulk)

f_ws = zeros(1,nz); %Price to SELL product i at retail.

f_ws([8,9])=[30,20]; %find(FG) = [8,9] <-> Roll/Bun, Burger, Egg Sandwich

f_w1 = [3,8,9/12,6.25/453,3.50,10,2,15,9,8,12]; %Premium retail prices (e.g. Belcampo, Grist and Toll, Frontier Co-Op, Mountain Rose Herbs
c_w1 = [1,1,12,453,1.5,1,1,1,1,1,1];

f_w2 = [2.5,7,12/30,44/(453*10),2.75,7.5,.5,-1,-1,30/8,10]; %NOTE: 9999 TERMS SHOULD BE ENFORCED-IRRELEVANT (true for all FG-ONLY goods?)
c_w2 = [25,25,30,453*10,5,2,24,-1,-1,8,6];

f_w3 =[1,2500/415,45/180,50/(453*25),2,6.75,.30,-1,-1,24/8,6.65];
c_w3 = [200,415, 180, 453*25,50,100,144,-1,-1,32,7];

f_w1 = repmat(f_w1',1,nt,nloc); %f_w can be varied over time and location
f_w2 = repmat(f_w2',1,nt,nloc); %f_w can be varied over time and location
f_w3 = repmat(f_w3',1,nt,nloc); %f_w can be varied over time and location

f_ws = repmat(f_ws',1,nt,nloc);

f_w(:,:,:,1) = f_w1;
f_w(:,:,:,2) = f_w2;
f_w(:,:,:,3) = f_w3;

f_w_easy = squeeze(f_w(:,1,1,:));


c_w = [c_w1;c_w2;c_w3];


s_y = [1/4,1/4,1/60,1/4000,1/5,1/5,1/30,1/4,1/4,1/10,1/10]; %in cubic feet (1/[amt that can be stored per cubic foot])


%% Materials Costs (automatically generated)

total_materials_cost = (RBom-diag(RM==0).*eye(size(RBom)))*f_w_easy; %total_materials_cost(i,m) is the total cost of the materials required to produce 1 unit of product i (at location 1 at time-step 1, by default), with goods purchased at tier m

%% WC
%Replaced with RM, SA, and FG
retail = 1:nloc;
wc = zeros(nz,nloc);
wc(find(FG),retail) = 1; %wc(i,k)==1 iff SKU i can be sold at location k
%% ZC
%For now, all locations are identical and can do every production process.
zc = zeros(nz,nloc);

%Identical-Facility Retail-Production (IFRP) Chain
for i = 1:nz
    if (RM(i) == 0)
        zc(i,:) = ones(1,nloc);
    end
end


%% Transportation and location-based parameters R_t, f_t; c_y,f_y,coordinates,locNames
%There are 12 cities

weeklyHoldingCost = .1; %per "unit" (cost to hold ONE unit for ONE week - one unit might be a pound of flour. It may be very cheap to hold one pound of flour for a week)
per_timeStep_holding_cost = weeklyHoldingCost*(time_step_hours/(7*24));

c_y=1200*ones(nloc,1);%[1200,1100,1000,900,800,700,600,500,400,300,200,100]';%Capacity to hold inventory (in "units")
f_y=per_timeStep_holding_cost*ones(1,nloc); %Cost to hold inventory (per "unit")
f_y = repmat(f_y,nt,1);

locNames = {'Rochester, NY','New York, NY','Boston, MA','Philadelphia, PA','Pittsburg, PA','Chicago, IL','Durham, NC','Long Beach, CA','Torrance, CA','Omaha, NB','Dallas, TX', 'Denver, CO'};
locNames = locNames(1:nloc);


%Querying Google API
% Querying Google Maps
%Google Maps API Credential: AIzaSyDDKHV64CvECHE8wf_KXpLiWr2fM0XAnrU

APIkey='AIzaSyAQ8hF_OJGNNSP9ekHrczqnUCsSPuznJ28';
pre_url = 'https://maps.googleapis.com/maps/api/distancematrix/json?';


d_t = zeros(nloc,nloc);
R_t = zeros(nloc,nloc);

inds = cell(nloc,nloc);
for k = 1:nloc
    for L = 1:nloc
        inds(k,L) ={[k,L]};
    end
end


for k = 1:ceil(length(locNames)/10)
    for L = 1:ceil(length(locNames)/10)
        mink = 10*(k-1)+1;
        minL = 10*(L-1)+1;
        limk = min(10*k,length(locNames));
        limL = min(10*L,length(locNames));
        
        locStringk = strrep(strjoin(locNames(mink:limk),'|'),' ','+');
        locStringL = strrep(strjoin(locNames(minL:limL),'|'),' ','+');

        
        url = [pre_url,'origins=',locStringk,'&','destinations=',locStringL,'&key=',APIkey];
        %url = 'https://maps.googleapis.com/maps/api/distancematrix/json?origins=nottingham&destinations=london|manchester|liverpool&key=AIzaSyDDKHV64CvECHE8wf_KXpLiWr2fM0XAnrU';
        result = webread(url);
        
        %To account for API overuse...
        if (size(result.rows,1)<1)
            d_t = 100000*ones(nloc,nloc);
            R_t = 100000*ones(nloc,nloc);
        else
            d_t(mink:limk,minL:limL)= cellfun(@(x) result.rows(x(1)).elements(x(2)).distance(1).value,inds(1:limk-mink+1,1:limL-minL+1));
            R_t(mink:limk,minL:limL)= cellfun(@(x) result.rows(x(1)).elements(x(2)).duration(1).value,inds(1:limk-mink+1,1:limL-minL+1));
        end
    end
end
clear mink; clear minL; clear limk; clear limL;


d_t = d_t*(1/1609.34); %Convert meters to miles.
R_t = R_t*(1/3600)*(1/time_step_hours);


R_t = ceil(R_t); %Integer number of time-steps
R_t(logical(eye(size(R_t)))) = 1; %Takes one time-step to ship from location to self! Thus, trucks can be used for storage.
f_t = fixedTruckCost*ones(nloc,nloc)+hourlyTruckCost*R_t*time_step_hours + perMileTruckCost*d_t; %Set f_t based on pre-rounded trucking times.

% Copying API Result into shipping times R_t and costs f_t. NOTE: (k,L)'th entry denotes k->L route.

% perMileTruckCost = .2;
% hourlyTruckCost = .2;
% fixedTruckCost = 1.5;

% 
% R_t = zeros(nloc,nloc);
% d_t = zeros(nloc,nloc);
% for k = 1:nloc
%     for L = 1:nloc
%         d_t(k,L)=result.rows(k).elements(L).distance(1).value; %Distance of route in meters!
%         R_t(k,L) = result.rows(k).elements(L).duration(1).value; %Duration of route in SECONDS!
%     end
% end

%Coordinates
%Note that 'locations{2}' contains the location-names as inputted by user
APIkey2='AIzaSyDVZPQrCVHtayFWu2LIWid6126hujMG5kQ';
pre_url2 = 'https://maps.googleapis.com/maps/api/geocode/json?';
coordinates = cell(1,size(locNames,1));
for i = 1:length(locNames)
    address = ['address=',strrep(locNames{i},' ','+')];
    url2 = [pre_url2,address,'&key=',APIkey2];
    %url = 'https://maps.googleapis.com/maps/api/distancematrix/json?origins=nottingham&destinations=london|manchester|liverpool&key=AIzaSyDDKHV64CvECHE8wf_KXpLiWr2fM0XAnrU';
    webreturn = webread(url2);
    webreturn = webreturn.results.geometry;
    webreturn = webreturn.location;
    coordinates{i}= webreturn;
end
coordinates = cell2mat(cellfun(@(x) [x.lat;x.lng],coordinates,'UniformOutput',0)); %coordinates(1,k) is the latitude of the k'th location; coordinates(2,k) is the longitude.


%% TC


% % Visualization Helper
% locTable = num2cell(R_t);
% locTable=[locNames;locTable];
% locTable = [['x';locNames'],locTable];
% locTable = [num2cell(0:nloc);locTable];
% locTable{1,1}='Location ID';

%TO AUTOMATE HUB DETERMINATION: 1) define number of hubs, 2) implement 
%K-Means (K= # hubs) clustering on 'coordinates' array, 3) define hubs to
%be locations within cluster with minimal sum of squares distances to other
%points in cluster

%A location's hub is the hub closest to it (by drive-cost from location to hub in f_t)
%Consider changing this to R_t
if (nloc > 10)
hubs = [4,9,10]; %Philadelphia, Torrance, Omaha
else
    hubs = 1:nloc;
end
myhub=hubs(arrayfun(@(x) find(f_t(x,hubs)==min(f_t(x,hubs)),1,'first'),(1:nloc)));


tc = zeros(nloc,nloc); %tc(k,l) = 0 iff route from k to l is never travelled
% Locations with the same hub can ship amongst themselves, and hubs 
%can ship between each other. Other routes are not allowed.
for k = 1:size(tc,2)
    tc(k,k)=1;
    tc(k,find(myhub==myhub(k)))=1;
end
%Hubs can ship between each other
for k = 1:length(hubs)
    for L = 1:length(hubs)
        tc(hubs(k),hubs(L))=1;
    end
end

%% All Production Parameters: f_z, f_z0, c_z, R_z
%All sites can produce identically
%rows [5,6,7,8,9,10] consist of NON-raw materials, aka MG (manufactured goods)

%Production Capacity
c_z = zeros(nz,1);%c_z(i,k) = number of units of product i produced in one production cycle at location k
c_z = [0,0,0,0,...%Wheat berries,Beef,Egg,Salt
    275,...%Flour (lb)
    200,...%Ground Beef (lb)
    500,...%Roll/Bun (ea)
    75,...%Burger (ea)
    60,...%Egg Sandwich (ea)
    250,...%Mayonnaise (lb) %~31 gallons! ~7.8 gallons per hour, or ~8 minutes per gallon..
    0]';
c_z = repmat(c_z,1,nloc);%Identical locations

%Production Time
time_z = zeros(nz,1);%Time to produce (one batch of) good i at location k, measured in HOURS.
time_z=repmat(time_z,1,nloc);
time_z = [0,0,0,0,...%Wheat berries,Beef,Egg,Salt
    4,...%Flour (lb)
   4,...%Ground Beef (lb)
    4,...%Roll/Bun (ea)
   4,...%Burger (ea)
    4,...%Egg Sandwich (ea
   4,...%Mayonnaise (lb)
    0]'; %Currently, "one batch" is defined as the amount that can be made in a 4-hour time-step.
R_z = zeros(nz,1); %R_z(i,k) = NUMBER OF TIME-STEPS REQUIRED TO PRODUCE a batch of PRODUCT i AT LOCATION k.
R_z = time_z./time_step_hours;
R_z = ceil(R_z);
R_z = repmat(R_z,1,nloc);
R_z = R_z + (R_z==0); %Ensure no production process takes zero time-steps.


%Production Costs
hourly_wage = 10; %Hourly wage per employee
numEmployees = ones(nz,1); %number of employees required in production of each sku.

%FIXED COSTS
f_z0 = zeros(nz,nloc); %f_z0(i,j,k)=fixed cost of starting production of product i at time-step j at location k
f_z0 = hourly_wage*time_step_hours*R_z;
f_z0= repmat(f_z0,1,1,nt);%Add time-dimension to cost matrix
f_z0=permute(f_z0,[1,3,2]);%Arrange matrix to be (nz,nt,nloc);

% f_z0=[zeros(1,12);zeros(1,12);zeros(1,12);zeros(1,12);...
%     7*ones(1,12);...%SKU 5: Wheat Flour
%     5*ones(1,12);...%SKU 6: Ground Beef (lb)
%     10*ones(1,12);... %SKU 7: Roll/Bun (ea)
%     0*ones(1,12);... % SKU 8: Burger (ea)
%     0*ones(1,12);... %SKU 9: Egg Sandwich (ea)
%     7*ones(1,12);
%     zeros(1,12)]; %SKU 10: Mayonnaise (lb)

%VARIABLE COSTS
% f_z = zeros(nz,nloc); %f_z(i,j,k) = PRICE OF PRODUCTION per unit of product i, at location k, at time-step j
% f_z=[zeros(1,12);zeros(1,12);zeros(1,12);zeros(1,12);...
%     .2*ones(1,12);...%SKU 5: Wheat Flour
%     1*ones(1,12);...%SKU 6: Ground Beef (lb)
%     .1*ones(1,12);... %SKU 7: Roll/Bun (ea)
%     1.5*ones(1,12);... % SKU 8: Burger (ea)
%     1.5*ones(1,12);... %SKU 9: Egg Sandwich (ea)
%     .5*ones(1,12);...
%     zeros(1,12)]; %SKU 10: Mayonnaise (lb)
% f_z = repmat(f_z,1,1,nt);%Add time-dimension to cost matrix
% f_z = permute(f_z,[1,3,2]);%Arrance matrix to be (nz,nt,nloc);

f_z = zeros(nz,nt,nloc);




%% Initialization: Storing inital values
init = zeros(9,nz,nt,nloc,max(nTiers,nloc)); %init(J(n,i,j,k,L)) = initial state of system, and MUST be defined for all j<=nt0
%Making a dummy initialization matrix, with nt0 = 1, here:
%Nothing in production, nothing in transit
%Alternatively, inventory OR sales OR production quantities at various
%times can be hard-set to meet orders
for j = 1:nt0
    for i = 1:nz
        for k=1:nloc
            init(1,i,j,k,1)=0; %10*(i<(floor(nz/3)));%Only allow initial storage of floor(nz/3) products until real parameters used.
            init(2,i,j,k,1)=0;%z(i,j,k)
            for L = 1:nloc
                init(3,i,j,k,L)=0;%t(i,j,k,L)
            end
            init(4,i,j,k,1)=0;%w(i,j,k)
            init(5,i,j,k,1)=0;%z0(i,j,k)
        end
    end
    for k = 1:nloc
        init(7,1,j,k,1)=0; %z00(j,k)
        for L = 1:nloc
            init(6,1,j,k,L)=0; %t0(j,k,L)
        end
    end
end


%% Storing Output in 'params' struct
params = struct();

params.nz = nz;
params.nt = nt;
params.nloc = nloc;
params.time_step_hours = time_step_hours;
params.nTiers = nTiers;
params.nt0 = nt0;

params.c_truck = c_truck;
params.perMileTruckCost = perMileTruckCost;
params.hourlyTruckCost = hourlyTruckCost;
params.fixedTruckCost = fixedTruckCost;
params.mu_w = mu_w;

params.Bom = Bom;
params.Gbom = Gbom;
params.SKU_names = SKU_names;
params.RM= RM;
params.FG = FG;
params.SA = SA;
params.RBom = RBom;
params.RGBom = RGBom;

params.s_y = s_y;
params.f_ws = f_ws;
params.f_w = f_w;
params.c_w = c_w;
params.f_w1 = f_w1;
params.f_w2 = f_w2;
params.f_w3 = f_w3;
params.c_w1 = c_w1;
params.c_w2 = c_w2;
params.c_w3=c_w3;
params.f_w_easy = f_w_easy;
params.total_materials_cost = total_materials_cost;

params.retail = retail;
params.wc = wc;

params.zc = zc;

params.c_y = c_y;
params.f_y = f_y;
params.locNames = locNames;
params.R_t = R_t;
params.d_t = d_t;
params.f_t = f_t;
params.coordinates = coordinates;

params.f_z0 = f_z0;
params.f_z = f_z;
params.R_z = R_z;
params.c_z = c_z;

params.tc = tc;

params.init = init;

save('prev_parameters.mat');
end
