close all
clear all

%% Load Demand Forecast
load('forecast_workspace.mat')
clearvars -except net extraTS extraTimes xlTimes realTS xlSkus trainForesight


%% Generating Parameter Templates (call 'CSV_Parameterizer.m')
 CSV_Parameterizer;

 
%% Cutting Data per model size
%load('/Users/Zach/Desktop/Theo K./forecast_workspace.mat')
% nz = 7; %Choose number of products to optimize for! As of now, the maximum calculable on my computer is approximately 30 with 4 locations; 100 with 2 locations
% nt = 15; %Choose a number of time-steps to optimize - take these two lines out eventually
%nt = size(extraTS,2);%Use all forecasted time-steps in optimization

nt0 = max([max(max(R_t)),max(max(R_z))]); %Number of time-steps included as initial conditions.

extraTS = extraTS(:,1:(nt-nt0));
shiftFactor = abs(min(min(min(extraTS)),min(min(realTS))));
extraTS = extraTS + shiftFactor;%Imperfect solution to negative demand predictions
realTS = realTS + shiftFactor;%Imperfect solution to negative demand predictions
extraTS = extraTS(1:nz,:);%Optimize based on nz products.
realTS = realTS(1:nz,:);
xlSkus = xlSkus(1:nz);

%% Randomization Parameters %Uncomment if not using CSV_Parameterizer.m
% %Note - these parameters are simply to model a plausible supply chain,
% %and will be irrelevant once real supply chain parameters are included
% numProdSites = 3;
% numSalesSites=2;
% nloc = numProdSites + numSalesSites;
%% Setting nt and nt0
%UNCOMMENT WHEN REAL SYSTEM BEING MODELED
%nloc = 5; %Number of locations (including production, warehousing, retail)
% nt0 = 5; %Number of time-steps included as initial conditions.
% nt = size(extraTS,2)+nt0; %Number of time-steps included

ntOld = size(realTS,2);

dimTotal = (4+nloc)*nz*nloc*nt + (nloc^2)*nt + nt*nloc;
dimReal = (nloc+3)*nz*nloc*nt;
dimInt = nz*nloc*nt + (nloc^2)*nt + nt*nloc;
%dimFinal = 4*nt*nloc*nz + 6*nt*nloc + 2*nt*ceil(1.5*nz);
%% Index Matrix J
%{
%Index Description
Populate Index Matrix J such that:
J(1,i,j,k,1) = io y(i,j,k) = quantity of i held at location k at time j
J(2,i,j,k,1) = io z(i,j,k) = quantity of i produced at location k at time j
J(3,i,j,k,L) = io t(i,j,k,L) = quantity of i transported from location k  to L at time j
J(4,i,j,k,1) = io w(i,j,k) = quantity of i sold at location k at time j
J(5,i,j,k,1) = io z0(i,j,k) = (1 or 0) whether i is produced at location k
at time j
J(6,1,j,k,L) = t0(j,k,L) = number of trucks from location k to L at time
j
J(7,i,j,k,1) = io z00(i,j,k) = (1 or 0) whether i is "in production"
(at least one time-step after start)at location k and time j

(io means "index of")
%}
J = zeros(7, nz, nt, nloc, nloc); %Index matrix
Jinv = zeros(dimTotal,5); %Inverse index matrix s.t. Jinv(J(n,i,j,k,L),:) = (n,i,j,k,L)
intcon = zeros(dimInt,1);
count2=1;%Counts number of integers
count = 1;
for nmat = 1:2
    for i = 1:nz
        for j = 1:nt
            for k = 1:nloc
                J(nmat,i,j,k,1)= count; %Index all y(i,j,k) [A1], z(i,j,k) [A2]
                Jinv(count,:) = [nmat,i,j,k,1];
                count = count + 1;
                
            end
        end
    end
end
for nmat = 3
    for i = 1:nz
        for j = 1:nt
            for k = 1:nloc
                for L = 1:nloc
                    J(nmat,i,j,k,L)= count; %Index all t(i,j,k,L) [A3]
                    Jinv(count,:) = [nmat,i,j,k,L];
                    count = count + 1;
                end
            end
        end
    end
end
for nmat = 4
    for i = 1:nz
        for j = 1:nt
            for k = 1:nloc
                J(nmat,i,j,k,1)= count; %Index all w(i,j,k) [A4], z0(i,j,k) [A5]
                Jinv(count,:) = [nmat,i,j,k,1];
                count = count + 1;
            end
        end
    end
end
for nmat = 5
    for i = 1:nz
        for j = 1:nt
            for k = 1:nloc
                J(nmat,i,j,k,1)= count; %Index all z0(i,j,k) [A5]
                Jinv(count,:) = [nmat,i,j,k,1];
                intcon(count2)=count;
                count2 = count2+1;
                count = count + 1;
            end
        end
    end
end
for nmat = 6
    for j = 1:nt
        for k = 1:nloc
            for L = 1:nloc
                J(nmat,1,j,k,L) = count; %Index all t0(j,k,L) [A6]
                Jinv(count,:) = [nmat,1,j,k,L];
                intcon(count2)=count;
                count2 = count2+1;
                count = count+1;
            end
        end
    end
end
for nmat = 7
    for j = 1:nt
        for k = 1:nloc
            J(nmat,1,j,k,1)=count; %Index all z00(j,k) [A7]
            Jinv(count,:)=[nmat,1,j,k,1];
            intcon(count2)=count;
            count2=count2+1;
            count=count+1;
        end
    end
end
clear count; clear nmat; clear i; clear j; clear k; clear L;
%% Parameter Arrays %Uncomment if not using CSV_Parameterizer.m
% 
% 
% %dem = 47*ones(nz,nt,nloc); %dem(i,j,k) = demand for product i, at time-step j, at location k
% Bom = triu(ones(nz,nz),1).*randi([0,1],nz,nz); %Bom(i,I) = number of units of item I required to produce 1 unit of item i DIRECTLY (aka "B").
% %Note - Direct vs Indirect Requirements:
% %{
% Let R(i,I) = number of units of item I required to produce 1 unit of item i
% DIRECTLY AND INDIRECTLY. By convention, B(i,i) = 0, R(i,i) = 1 for all i.
% Then R = (I - B)^(-1), where I is the nz x nz identity matrix.
% %}
% %To incorporate R vs B: Use R to calculate "total cost per produced unit" of each producable good
% %Then, by convention set f_w = 2*(total cost), as a "good" dummy parameter.
% 
% 
% R_z = 3*ones(nz,nloc); %R_z(i,k) = NUMBER OF TIME-STEPS REQUIRED TO PRODUCE PRODUCT i AT LOCATION k
% R_t = 4*ones(nloc, nloc); %R_t(k,L) = NUMBER OF TIME-STEPS REQUIRED TO TRANSPORT FROM LOCATION k TO LOCATION L
% 
% s_y = 1*ones(nz,1); %s_y(i) denotes the number of "units" (of space, load, etc.) required to hold one unit of product/component i
% 
% c_y = 250*ones(nloc,1); %c_y(k) = the storage capacity (in "units") at location k
% c_truck = 75; %c_truck is the number of "units" that a truck can fit
% %c_z = 10*ones(nz,nloc); %c_z(i,k) = capacity to produce product i at location k in one production cycle
% c_z = [250*randi([0,1],nz,numProdSites),zeros(nz,nloc-numProdSites)]; %nz X nloc
% %c_z = 10*ones(nz,nloc); %Extra easy dummy
% f_z0 = 5*ones(nz,nt,nloc); %f_z0(i,j,k)=fixed cost of starting production of product i at time-step j at location k
% f_z = .1*ones(nz,nt,nloc); %f_z(i,j,k) = PRICE OF PRODUCTION per unit of product i, at location k, at time-step j
% f_y = .1*ones(nt,nloc); %f_y(j,k) = PRICE OF HOLDING one "unit" at a given location, at a given time (possibly not so variable over time)
% f_t = 10*ones(nloc,nloc); %f_t(k,L) = price of sending one truck from location k to location L. If f_t(k,L) <0, then no truck will make that route. By convention, set to -1.
% dummy1= randi([0,1],nloc,nloc); %Experiment
% dummy2= randi([0,1],nloc,nloc); %Experiment
% f_t = f_t - 2*(dummy1 + dummy2 - dummy1.*dummy2); %Experiment - want lots of non-travellable routes.
% f_t = 5*f_t;
% clear dummy1; clear dummy2;
% 
% f_w = 25*ones(nz,nt,nloc); %f_w(i,j,k) = REVENUE per sale of product type i at
% f_wf = f_w*(1/nt);
% %time-step j @ location k
% 
% %Note - Pricing Variability:
% %{
% Note: Can consider different price points to be separate
% "products," specify demand at different price points, and
% constrain such that only one price point can be used using
% binary constraints (as in discrete production
% variables), or even one price point PER LOCATION if
% desired
% %}
% 
% mu_w = .01*ones(nz,nloc); %Sales backlog decay constants.
% perMileTruckCost = .2;
% hourlyTruckCost = .2;
% fixedTruckCost = 1.5;
% 


%% Classification parameter array %Uncomment if not using CSV_Parameterizer.m
% %Use this when I start removing unused variables and create upper and lower
% %bounds.
% wc = zeros(nz,nloc); %%wc(i,k)= 1 <=> product i can be sold at location k; wc(i,k) = -1 <=> i can be bought at k; wc(i,k)=0 otherwise
% % wc = wc + randi([0,1],nz,nloc); %Experiment
% % wc= wc - randi([0,1],nz,nloc); %Experiment
% % wc(1:floor(nz/2),[1,2])=1;
% % wc(floor(nz/2):end,3)=1;
% % wc= wc.*(c_z==0);%Arbitrary BUT this ensures that locations
% % %that sell or buy something aren't the ones that can produce it. And not
% % %all of them can produce it (only those that intersect with the random 1's
% % %in this randi matrix).
% % %wc = ones(nz,nloc);%Extra dumb dummy
% % wc = wc - randi([0,1],size(wc));
% wc(:,(numProdSites+1:end)) = randi([0,1],size(wc,1),nloc-numProdSites);
% 
% zc = (c_z >0); %zc(i,k) ==1 iff capacity to produce c_z(i,k) is > 0
% 
% tc = ones(nloc,nloc); %tc(k,L) = 1 if trucking can happen between locations k and L, 0 otherwise.
% tc = tc - (f_t<0);


zNonProducers = ~any(c_z,1); %Locations that have NO production can be stripped of their z00 variable at all time-steps


zPrune = nnz(zc<1); %Number of pairs of product and location at which production will never happen.
tPrune = nnz(tc<1); %Number of pairs of locations between which no trucks will every be sent.
wPrune = nnz(wc==0); %Number of pairs of product and location at which sales and purchases will never occur.
nPrune = 2*zPrune*nt + (tPrune*nt + tPrune*nt*nz) + wPrune*nt +nnz(zNonProducers)*nt; %number of y,z,t,w,z0,t0,z00 variables
intPrune = zPrune*nt + tPrune*nt + nnz(zNonProducers)*nt; %number of pruned z0, z00, and t0 variables


%% Imported Demand Forecast
dem = zeros(nz,nt,nloc);
for i = 1:nz
    %divFactor = 1./nnz(wc(i,:)==1); %Distribute demand equally over every location that can sell the product. Change when demand is resolved by location.
    divFactor = 1/nloc;
    for j = 1:nt0
        for k = 1:nloc
            dem(i,j,k) = realTS(i,ntOld-nt0+j)*divFactor*(wc(i,k)>0);
        end
    end
    for j = (nt0+1):nt
        for k = 1:nloc
            dem(i,j,k) =  extraTS(i,j-nt0)*divFactor*(wc(i,k)>0);
        end
    end
end

%% Bounds:
%Lower: >=0 for all (except <=0 for acquisitions - this requires NO constraint!)
%Upper Bounds: variable classification (and demand)

%If you can buy, you can buy infinitely much
%Can't sell beyond demand! (best method?)
LB = zeros(dimTotal,1);
UB = Inf*ones(dimTotal,1);
for i = 1:nz
    for k = 1:nloc
        %NOTE: Cannot BUY and SELL at same location! WC==1, or WC==-1,
        %or WC==0. Not configurable for it to be otherwise.
        if (wc(i,k) == -1) %Can BUY product i at location k
            for j = 1:nt
                %LB(J(4,i,j,k,1)) = -Inf; %Can buy infinitely much product
                LB(J(4,i,j,k,1)) = -c_y(k)/s_y(i); %Can buy only what can be held at location.
                UB(J(4,i,j,k,1)) = 0; %Can't sell any product.
            end
        elseif (wc(i,k)==1) %Can SELL product i at location k
            %Sales default maximum is Inf given demand backlog constraints.
            for j = 1:nt
                UB(J(4,i,j,k,1)) = 0;
                for T = 1:j
                    UB(J(4,i,j,k,1)) = UB(J(4,i,j,k,1))+ mu_w(i,k)^(j-T)*dem(i,T,k); %Demand Backlog constraint in SALES equation.
                    %UB(J(4,i,j,k,1)) = dem(i,j,k); %w(i,j,k) <= dem(i,j,k)
                end
                UB(J(4,i,j,k,1)) = min(UB(J(4,i,j,k,1)),c_y(k)/s_y(i));
            end
        else %wc(i,k)==0 %Can NEITHER SELL NOR BUY product i at location k
            for j = 1:nt
                LB(J(4,i,j,k,1))=0;
                UB(J(4,i,j,k,1))=0; %Can't sell at all if classified as non-retail
            end
        end
    end
end





%Upper Bounds

%Can't ship to self!
for j = 1:nt
    for k = 1:nloc
        UB(J(6,1,j,k,k))=0; %t_0(j,k,k) <= 0 forall j,k
        for i = 1:nz
            UB(J(3,i,j,k,k)) = 0; %t(i,j,k,k) <=0 forall i,j,k
        end
    end
end


%Can't produce more than capacity
%NOTE: This is accounted for in PRODUCTION constraint
%z(i,j,k)<=z0(i,j,k)c_z(i,k),but may increase efficiency by reducing
%feasible region.
for i=1:nz
    for j = 1:nt
        for k = 1:nloc
            UB(J(2,i,j,k,1))=c_z(i,k); %z(i,j,k)<=c_z(i,k)
        end
    end
end

%Holding: Can't hold more than capacity
%NOTE: This is accounted for in the INVENTORY/HOLDING constraint
%Sum_i y(i,j,k) <= c_y(k) forall j,k
%since y(i,j,k)>=0 forall i,j,k
for i = 1:nz
    for j = 1:nt
        for k = 1:nloc
            UB(J(1,i,j,k,1))=c_y(k); %y(i,j,k)<=c_y(k)
        end
    end
end

%Binary variables ARE binary
for i = 1:nz
    for j = 1:nt
        for k = 1:nloc
            UB(J(5,i,j,k,1)) = 1; %z_0(i,j,k)<=1 NOTE:
            %this is NOT necessary (it will be 0 or 1 due to optimality,
            %but this might help the algorithm proceed more quickly/efficiently.
            %Same goes for z00(i,j,k)<=1
        end
    end
end
for n = 7
    for j = 1:nt
        for k = 1:nloc
            UB(J(n,1,j,k,1)) = 1; %z_00(i,j,k)<=1
        end
    end
end

%% Initialization: Storing inital values
init = zeros(size(J)); %init(J(n,i,j,k,L)) = initial state of system, and MUST be defined for all j<=nt0
%Making a dummy initialization matrix, with nt0 = 1, here:
%Nothing in production, nothing in transit
%Alternatively, inventory OR sales OR production quantities at various
%times can be hard-set to meet orders
for j = 1:nt0
    for i = 1:nz
        for k=1:nloc
            init(1,i,j,k,1)=10*(i<(floor(nz/3)));%Only allow initial storage of floor(nz/3) products until real parameters used.
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
%% Initialization: hard-constraining values
% Let init be a matrix the same size as J (i.e., accommodates all decision
% variables) consisting of all zeros EXCEPT for entries init(n,i,j,k,L)
% where j<=nt0, which means the time-step is included in initialization. In
% that case, init(n,i,j,k,L) contains the initialized value of the decision
% variable being referenced.
for j = 1:nt0
    for k = 1:nloc
        for n = [1,2,4,5] %Initialize y,z,t,w,z0,t0,z00 (i,j,k)
            for i = 1:nz
                UB(J(n,i,j,k,1))= init(n,i,j,k,1);
                LB(J(n,i,j,k,1))= init(n,i,j,k,1);
            end
        end
        for n = 3 %Initialize t(i,j,k,L)
            for i = 1:nz
                for L = 1:nloc
                    UB(J(n,i,j,k,L))= init(n,i,j,k,L);
                    LB(J(n,i,j,k,L))= init(n,i,j,k,L);
                end
            end
        end
        for n = 6 %Initialize t0(j,k,L)
            for L = 1:nloc
                UB(J(n,1,j,k,L))= init(n,1,j,k,L);
                LB(J(n,1,j,k,L))= init(n,1,j,k,L);
            end
        end
        for n = 7 %Initialize z00(j,k)
            UB(J(n,1,j,k,1))=init(n,1,j,k,1);
            LB(J(n,1,j,k,1))=init(n,1,j,k,1);
        end
    end
end
%% Objective Function: qOb*x
%Note that this objective function, which maximizes profit (sales revenue-costs)
%DOES NOT maximize the value of the inventory that is IN TRANSIT or BEING
%HELD. This is an extreme simplification that can be accounted for easily.

qOb = zeros(1,dimTotal);
for i = 1:nz
    for j = 1: nt
        for k = 1:nloc
            qOb(J(1,i,j,k,1)) = -f_y(j,k)*s_y(i); %holding cost
            qOb(J(2,i,j,k,1)) = -f_z(i,j,k); %proportional price of production
            qOb(J(4,i,j,k,1)) = f_w(i,j,k); %sales revenue/purchasing cost
            qOb(J(5,i,j,k,1)) = -f_z0(i,j,k);%fixed cost of production
            for L = 1:nloc
                qOb(J(6,1,j,k,L)) = -f_t(k,L); %price of sending trucks
            end
        end
    end
end
qOb = -qOb; %intlinprog MINIMIZES objective; qOb was formulated for maximization.

%% Product Equity (Added to objective function)
qEquity = zeros(1,dimTotal);

qOb = qOb - qEquity;

%% Evolution Equation: qEv*x == 0 forall i, j, k
qEv = zeros(nz*(nt-nt0)*nloc,dimTotal);
bEv = zeros(nz*(nt-nt0)*nloc,1); %Equality constraint with no constant term
count = 1;
for i = 1:nz
    for j = (nt0+1):nt %INITIALIZE AT j=nt0 !!!!
        for k = 1:nloc
            q = zeros(1,dimTotal);
            q(J(1,i,(j-1),k,1)) = 1; %holdover
            q(J(4,i,(j-1),k,1)) = -1; %sales
            if (j>R_z(i)) %if product can come out of production at this time-step
                q(J(2,i,(j-R_z(i,k)),k,1)) = 1; %product i @ location k coming out of production at time j, which started at time-step j - R_z(i,k) by definition
            end
            for I = 1:nz
                q(J(2,I,j-1,k,1)) = -Bom(I,i); %product i @ location k being used for production of product I at time-step (j-1)
            end
            for I = 1:nloc
                if (j>R_t(I,k))
                    q(J(3,i,j-R_t(I,k),I,k)) = 1; %transported inventory: incoming
                end
                q(J(3,i,j-1,k,I)) = -1; %transported inventory: outgoing
            end
            q(J(1,i,j,k,1)) = -1; %Set equality for current (i,j,k) in equation
            qEv(count,:) = q;
            count = count +1;
        end
    end
end
clear count; clear I; clear q;

%% Evolution FINAL-TIMESTEP Equation q*x <= 0
%w(i,nt,k) + Sum_I z(I,nt,k)Bom(I,i) + Sum_L t(i,nt,k,L) + w(i,nt,k) -
%y(i,nt,k) <= 0
qEvF = zeros(nz*nloc,dimTotal);
bEvF = zeros(nz*nloc,1);
count = 1;
for i = 1:nz
    for k = 1:nloc
        q = zeros(1,dimTotal);
        q(J(1,i,nt,k,1)) = -1; %y(i,nt,k)
        q(J(4,i,nt,k,1))=1; %w(i,nt,k)
        for I = 1:nz
            q(J(2,I,nt,k,1))=Bom(I,i);  %Bom(I,i)*z(I,nt,k)
        end
        for I = 1:nloc
            q(J(3,i,nt,k,I))=1; %Sum_L t(i,nt,k,I)
        end
        qEvF(count,:) = q;
        count = count +1;
    end
end
clear count; clear i; clear k; clear q; clear I;
%% Inventory/Holding: Sum_i y(i,j,k) <= c_y(k) forall j,k
qInv = zeros(nt*nloc,dimTotal);
bInv = zeros(nt*nloc,1); %Less than constraint
count = 1;
for j = 1:nt
    for k = 1:nloc
        q = zeros(1,dimTotal);
        for i = 1:nz
            q(J(1,i,j,k,1))=s_y(i); %Each item held at a location at a time step y(i,j,k) contributes to "total units held", which much be less than holding capacity at that location c_y(k)
        end
        qInv(count,:) = q;
        bInv(count) = c_y(k);
        count = count+1;
    end
end

clear count; clear i; clear j; clear k; clear q;
%% Production: z(i,j,k) <- c_z(i,k)z_0(i,j,k) forall i,j,k
qProd = zeros(nz*nt*nloc,dimTotal);
bProd = zeros(nz*nt*nloc,1); %Less than constraint

count = 1;
for i = 1:nz
    for j = 1:nt
        for k = 1:nloc
            q = zeros(1,dimTotal);
            q(J(2,i,j,k,1))=1; %z(i,j,k)
            q(J(5,i,j,k,1)) = -c_z(i,k); %-c_z(i,k)z_0(i,j,k)
            qProd(count,:) = q;
            %bProd(count) = 0;
            count = count+1;
        end
    end
end

clear i; clear j; clear k; clear L; clear q; clear count;
%% Transportation: Sum_i t(i,j,k,L)*s_y(i) <= c_truck*t_0(j,k,L) forall j,k,L
qTrans = zeros(nt*nloc*nloc,dimTotal);
bTrans = zeros(nt*nloc*nloc,1); %Less than constraint
count = 1;
for j = 1:nt
    for k = 1:nloc
        for L=1:nloc
            q = zeros(1,dimTotal);
            for i =1:nz
                q(J(3,i,j,k,L)) = s_y(i);
            end
            q(J(6,1,j,k,L)) = -c_truck;
            qTrans(count,:) = q;
            count = count +1;
        end
    end
end


clear i; clear j; clear k; clear L; clear q; clear count;
%% Sales (Backlog-Capable) qSales*x <= bSales
%Sum_{T=1}^j mu_w(i,k)^(j-T)w(i,T,k)  <= Sum_{T=1}^j mu_w(i,k)^(j-T)dem(i,T,k) forall i,j,k


qSales = zeros(nz*nt*nloc,dimTotal);
bSales = zeros(nz*nt*nloc,1); %Less than constraint

count = 1;
for i = 1:nz
    for j = 1:nt
        for k = 1:nloc
            q = zeros(1,dimTotal);
            for T =1:j
                q(J(4,i,T,k,1))=mu_w(i,k)^(j-T); %w(i,j,k)
                bSales(count) = bSales(count) + mu_w(i,k)^(j-T)*dem(i,T,k); %mu_w(i,k)*dem(i,T,k)
            end
            qSales(count,:) = q;
            %bProd(count) = 0;
            count = count+1;
        end
    end
end

clear i; clear j; clear k; clear L; clear q; clear count;

%% Production DISCRETE z_00(j,k)+ Sum_i z_0(i,j,k)  <= 1 forall j,k
qProd0 = zeros(nt*nloc,dimTotal);
bProd0 = ones(nt*nloc,1); %Less than constraint
count = 1;
for j = 1:nt
    for k = 1:nloc
        q = zeros(1,dimTotal);
        q(J(7,1,j,k,1)) = 1;
        for i = 1:nz
            q(J(5,i,j,k,1)) = 1;
        end
        qProd0(count,:) = q;
        count = count+1;
    end
end

clear i; clear j; clear k; clear L; clear q; clear count;
%% Production Discrete2 -z00(j,k) + Sum_i Sum_(1<=TT<=R_z(i,k)) z0(i,j-TT,k) == 0 forall j,k
% This constraint ensures that only one product can be "in production" at a
% time
qProd00 = zeros(nt*nloc,dimTotal);
bProd00 = zeros(nt*nloc,1); %Less than constraint
count = 1;
for j = 1:nt
    for k = 1:nloc
        q = zeros(1,dimTotal);
        q(J(7,1,j,k,1))=-1; %z00(j,k)
        for i = 1:nz
            for TT = 1:(R_z(i,k)-1) %%NOTE: sum runs to R_z(i,k)-1 instead of R_z(i,k) because
                %j-R_z(i,k) is the time-step of the production order, not an "in-production" time-step
                if (j>TT)
                    q(J(5,i,j-TT,k,1)) = 1;
                end
            end
        end
        qProd00(count,:) = q;
        %         bProd0(count) = 1;
        count = count+1;
    end
end



clear i; clear j; clear k; clear L; clear q; clear T; clear count;
%% Transportation DISCRETE t(i,j,k,k)=0 forall i,j,k && t_0(j,k,k)=0 forall j,k
%See upper bound
%% Concatenating constraints AKA final constraint matrix
qFinal = [qEvF; qInv; qProd; qTrans; qSales; qProd0]; %Complete constraint matrix
bFinal = [bEvF; bInv; bProd; bTrans; bSales; bProd0];
qEqFinal = [qEv;qProd00];
bEqFinal = [bEv;bProd00];
%% Pruning: Recording Indices
%TRANSPORT (cts & integer):

%{
If f_t(k,L) < 0 (set to -1 as convention), then remove columns
corresponding to each t(i,j,k,L) and t_0(j,k,L) forall i,j, which are
columns J(3,i,j,k,L) and J(6,1,j,k,L), respectively, forall i,j.
This corresponds to trucks NEVER travelling that route.
NOTE: Use tc for transportation classification
%}
pruner = zeros(nPrune,1);
count =1;
for k = 1:nloc
    for L = 1:nloc
        if (tc(k,L) ==0) %if transportation cannot happen between locations k and L
            for j =1:nt
                pruner(count) = J(6,1,j,k,L); %prune all unused t0(j,k,L) variables
                count = count+1;
                for i = 1:nz
                    pruner(count) = J(3,i,j,k,L); %prune all unused t(i,j,k,L) variables
                    count = count+1;
                end
            end
        end
    end
end


%SALES:
for i = 1:nz
    for k = 1:nloc
        if (wc(i,k) ==0) %if product i cannot be sold or bought at location k
            for j = 1:nt
                pruner(count) = J(4,i,j,k,1); %prune all unused w(i,j,k)
                count = count+1;
            end
        end
    end
end


%PRODUCTION (cts & binary):

for k = 1:nloc
    if (zNonProducers(k)==1)
        for j = 1:nt
            pruner(count) = J(7,1,j,k,1); %prune all unused z00(j,k) variables
            count = count+1;
        end
    end
    for i =1:nz
        if (zc(i,k) ==0)
            for j = 1:nt
                pruner(count) = J(2,i,j,k,1); %prune all unused z(i,j,k)
                count = count+1;
                pruner(count) = J(5,i,j,k,1); %prune all unused z0(i,j,k)
                count = count+1;
            end
        end
    end
end

%INVENTORY:
%Eventually, implement a "check" procedure to determine if a product will
%NEVER be at a given location (no incoming transport variables, no purchase
%variables, no production variables), in which case the inventory variables
%are automatically eliminated.
%Also implement usage of tc, the "route not travelled" parameter.

pruner = sort(pruner);
%% Reparameterization

J2 = zeros(7,nz,nt,nloc,nloc); %Index matrix
Jinv2 = Jinv;
Jinv2(pruner,:) = []; %Delete pruned indices
intcon2=zeros(dimInt-intPrune,1);
count2=1; %Index of
for count = 1:(dimTotal-nPrune)
    J2(Jinv2(count,1),Jinv2(count,2),Jinv2(count,3),Jinv2(count,4),Jinv2(count,5)) = count;
    if (Jinv2(count,1)>=5)
        intcon2(count2)=count;
        count2=count2+1;
    end
end
%J2(Jinv(intcon))
%intcon = intersect(setxor(intcon,pruner),intcon);
%% Pruning: Deleting Columns and rows
qFinal(:,pruner) = [];
qEqFinal(:,pruner) = [];
qOb(pruner) = [];
LB(pruner) = [];
UB(pruner) = [];


bFinal(~any(qFinal,2))=[]; %Delete unused rows (unused constraints)
qFinal(~any(qFinal,2),:) = []; %Delete unused rows (unused constraints)

% disp('Number of zero constraint rows with nonzero constant')
% nnz(~any(qFinal,2).*any(bFinal,2)) %Should be zero

%% Specifying which decisions are integer
%SEE REPARAMETERIZATION SECTION
%% MILP

%The objective function (which is minimized) is total COST. Therefore, if
%the optimal value is negative, there is a PROFIT. If it is positive, there
%is a deficit.

%--'OutputFcn',@savemilpsolutions "puts the feasible points in a matrix named
%xIntSol in your base workspace, where each column is one integer feasible 
%point. It saves the objective function values in a vector named fIntSol,"
%--'RootLPAlgorithm','primal-simplex' alternative is 'dual-simplex'

options = optimoptions('intlinprog','OutputFcn',@savemilpsolutions,'RootLPAlgorithm','primal-simplex');
%options.MaxTime = 7200; %Default
options.MaxTime = 15000; %4.16 hours
options.CutGeneration = 'advanced'; %Default is 'basic'
options.CutGenMaxIter = 50; %Default is 10;
x = intlinprog(qOb,intcon2,qFinal,bFinal,qEqFinal,bEqFinal,LB,UB,options);
%% Storing Returned Values & System Statistics
Y = zeros(nz,nt,nloc); %(i,j,k)'th entry is quantity of i held over at loc. k at time-step j
Z = zeros(nz,nt,nloc);%(i,j,k)'th entry is quantity of i produced at loc. k at time-step j
T = zeros(nz,nt,nloc,nloc);%(i,j,k,L)'th entry is quantity of i transported from loc. k to L at time-step j
W = zeros(nz,nt,nloc);%(i,j,k)'th entry is quantity of i sold (bought if neg.) at loc. k at time-step j
Z0 = zeros(nz,nt,nloc);%(i,j,k)'th entry is WHETHER OR NOT product i went into production at location k at time-step j (0=no, 1=yes)
T0 = zeros(nt,nloc,nloc);%(i,k,L)'th entry is number of trucks was sent from loc. k to loc. L at time-step j
Z00 = zeros(nt,nloc);%(j,k)'th entry is WHETHER OR NOT location k is "in production" at time-step j
track_Y = zeros(nnz(Y),4);%Row = [product, time-step,location,holdover]
track_Z = zeros(nnz(Z),4);%Row = [product,time-step,location,production order]
track_T = zeros(nnz(T),5);%Row = [product,time-step,location1,location2,shipment order]
track_W = zeros(nnz(W),4);%Row = [product,time-step,location,sales order]
track_Z0 = zeros(nnz(Z0),4);%Row = [product,time-step,location,"start production?" (1 or 0)]
track_T0 = zeros(nnz(T0),4);%Row = [time-step,location1,location2,number of trucks sent from 1 to 2]
track_Z00 = zeros(nnz(Z00),3);%Row = [time-step,location,"in production?" (1 or 0)]

viewOrderY = [3,4,2]; %[j,k,i,y(i,j,k)]
viewOrderZ = [3,4,2]; %[j,k,i,z(i,j,k)]
viewOrderT = [3,4,5,2]; %[j,k,L,i,t(i,j,k,L)]
viewOrderW = [3,4,2]; %[j,k,i,w(i,j,k)]
viewOrderZ0 = [3,4,2]; %[j,k,i,z0(i,j,k)]
viewOrderT0 = [3,4,5]; %[j,k,L,t0(j,k,L)]
viewOrderZ00 = [3,4]; %[j,k,z00(j,k)]

%System Statistics
salesOf = zeros(1,nz); %Revenues earned from each sku
unitsSold = zeros(1,nz); %Number of units sold of each SKU
salesAtLoc = zeros(1,nloc); %Revenues at each location
salesAtT = zeros(1,nt); %Revenues from each time-step
TCosts = 0;
ZCosts = 0;
Z0Costs = 0;
YCosts = 0;

YRow = (@(i) [arrayfun(@(n) Jinv2(i,n),viewOrderY),x(i)]);
ZRow = (@(i) [arrayfun(@(n) Jinv2(i,n),viewOrderZ),x(i)]);
TRow = (@(i) [arrayfun(@(n) Jinv2(i,n),viewOrderT),x(i)]);
WRow = (@(i) [arrayfun(@(n) Jinv2(i,n),viewOrderW),x(i)]);
Z0Row = (@(i) [arrayfun(@(n) Jinv2(i,n),viewOrderZ0),x(i)]);
T0Row = (@(i) [arrayfun(@(n) Jinv2(i,n),viewOrderT0),x(i)]);
Z00Row = (@(i) [arrayfun(@(n) Jinv2(i,n),viewOrderZ00),x(i)]);



counts = ones(7,1); %Index of how many of each type of decision has been dealt with. [Y,Z,T,W,Z0,T0,Z00]
%NOTE: Jinv2(i,:) = [n,i,j,k,L], n maps [1,2,3,4,5]->[Y,Z,T,W,Z0,T0,Z00]
for i = 1:size(Jinv2,1)
    if (x(i)~=0)
        %track_Y(i,:)= [sku,time-step,location,holding]
        if (Jinv2(i,1)==1)
            Y(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4)) = x(i);
            %track_Y(counts(Jinv2(i,1)),:) = [Jinv2(i,2),Jinv2(i,3),Jinv2(i,4),x(i)];
            track_Y(counts(Jinv2(i,1)),:) = YRow(i);
            YCosts = YCosts + x(i)*s_y(Jinv2(i,2))*f_y(Jinv2(i,3),Jinv2(i,4)); %+= y(i,j,k)*s_y(i)*f_y(j,k)
            counts(Jinv2(i,1)) = counts(Jinv2(i,1))+1;
        end
        %track_Z(i,:)= [sku,time-step,location,production]
        if (Jinv2(i,1)==2)
            Z(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4)) = x(i);
            %track_Z(counts(Jinv2(i,1)),:) = [Jinv2(i,2),Jinv2(i,3),Jinv2(i,4),x(i)];
            track_Z(counts(Jinv2(i,1)),:) = ZRow(i);
            ZCosts = ZCosts + x(i)*f_z(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4)); %+= z(i,j,k)*f_z(i,j,k)
            counts(Jinv2(i,1)) = counts(Jinv2(i,1))+1;
        end
        %track_T(i,:)= [sku,time-step,location1,location2,holding]
        if (Jinv2(i,1)==3)
            T(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4),Jinv2(i,5)) = x(i);
            %track_T(counts(Jinv2(i,1)),:) = [Jinv2(i,2),Jinv2(i,3),Jinv2(i,4),Jinv2(i,5),x(i)];
            track_T(counts(Jinv2(i,1)),:) = TRow(i);
            counts(Jinv2(i,1)) = counts(Jinv2(i,1))+1;
        end
        if (Jinv2(i,1)==4)
            W(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4)) = x(i);
            %track_W(counts(Jinv2(i,1)),:) = [Jinv2(i,2),Jinv2(i,3),Jinv2(i,4),x(i)];
            track_W(counts(Jinv2(i,1)),:) = WRow(i);
            unitsSold(Jinv2(i,2))= unitsSold(Jinv2(i,2)) + x(i); %+= w(i,j,k)
            salesOf(Jinv2(i,2)) = salesOf(Jinv2(i,2))+f_w(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4))*x(i); %+= f_w(i,j,k)*w(i,j,k)
            salesAtT(Jinv2(i,3)) = salesAtT(Jinv2(i,3)) + f_w(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4))*x(i);%+= f_w(i,j,k)*w(i,j,k)
            salesAtLoc(Jinv2(i,4)) = salesAtLoc(Jinv2(i,4))+ f_w(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4))*x(i);%+= f_w(i,j,k)*w(i,j,k)
            counts(Jinv2(i,1)) = counts(Jinv2(i,1))+1;
        end
        if (Jinv2(i,1)==5)
            Z0(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4)) = x(i);
            %track_Z0(counts(Jinv2(i,1)),:) = [Jinv2(i,2),Jinv2(i,3),Jinv2(i,4),x(i)];
            track_Z0(counts(Jinv2(i,1)),:) = Z0Row(i);
            Z0Costs = Z0Costs + f_z0(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4))*x(i); %+=f_z0(i,j,k)*z0(i,j,k)=f_z0(i,j,k)
            counts(Jinv2(i,1)) = counts(Jinv2(i,1))+1;
        end
        if (Jinv2(i,1)==6)
            T0(Jinv2(i,3),Jinv2(i,4),Jinv2(i,5)) = x(i);
            %track_T0(counts(Jinv2(i,1)),:) = [Jinv2(i,3),Jinv2(i,4),Jinv2(i,5),x(i)];
            track_T0(counts(Jinv2(i,1)),:) = T0Row(i);
            TCosts = TCosts + x(i)*f_t(Jinv2(i,4),Jinv2(i,5)); %+= T0(j,k,L)*f_t(k,L)
            counts(Jinv2(i,1)) = counts(Jinv2(i,1))+1;
        end
        if (Jinv2(i,1)==7)
            Z00(Jinv2(i,3),Jinv2(i,4)) = x(i);
            %track_Z00(counts(Jinv2(i,1)),:) = [Jinv2(i,3),Jinv2(i,4),x(i)];%[i,k,Z00(j,k)]
            track_Z00(counts(Jinv2(i,1)),:) = Z00Row(i);
            counts(Jinv2(i,1)) = counts(Jinv2(i,1))+1;
        end
    end
end

track_Y = sortrows(track_Y);
track_Z = sortrows(track_Z);
track_T = sortrows(track_T);
track_W = sortrows(track_W);
track_Z0 = sortrows(track_Z0);
track_T0 = sortrows(track_T0);
track_Z00 = sortrows(track_Z00);



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
labels = {'Holding Orders';'Production Orders';'Transportation Orders';'Sales Benchmarks (negative value = material purchase)';'Production Site Un-Availability';'Truck Orders';'Production (availability) Orders'};
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


fprintf(fid,'%s\n',['Final Holding Costs: ,', num2str(YCosts)]);
fprintf(fid,'%s\n',['Final Production Costs: ,', num2str(ZCosts + Z0Costs)]);
fprintf(fid,'%s\n',['Final Transportation Costs: ,', num2str(TCosts)]);
fprintf(fid,'%s\n',['Final Profit: ,',num2str(-qOb*x)]);

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

%% Saving solution for horizon shift/turnover
oldSol.MILP.x = x;
oldSol.MILP.qOb = qOb;
oldSol.MILP.intcon2 = intcon2;
oldSol.MILP.qFinal = qFinal;
oldSol.MILP.bFinal = bFinal;
oldSol.MILP.qEqFinal = qEqFinal;
oldSol.MILP.bEqFinal = bEqFinal;
oldSol.MILP.LB = LB;
oldSol.MILP.UB = UB;
oldSol.MILP.J2 = J2;
oldSol.MILP.Jinv2 = Jinv2;

oldSol.SOLUTION.W = W;
oldSol.SOLUTION.Z = Z;
oldSol.SOLUTION.T = T;
oldSol.SOLUTION.T0 = T0;
oldSol.SOLUTION.Z0 = Z0;
oldSol.SOLUTION.Z00 = Z00;
oldSol.SOLUTION.Y = Y;
oldSol.SOLUTION.track_T = track_T;
oldSol.SOLUTION.track_T0 = track_T0;
oldSol.SOLUTION.track_W = track_W;
oldSol.SOLUTION.track_Y = track_Y;
oldSol.SOLUTION.track_Z = track_Z;
oldSol.SOLUTION.track_Z0 = track_Z0;
oldSol.SOLUTION.track_Z00 = track_Z00;

oldSol.oldParam.Jinv = Jinv;
oldSol.oldParam.J = J;
oldSol.oldParam.intcon = intcon;
oldSol.oldParam.qEvF = qEvF;
oldSol.oldParam.qInv = qInv;
oldSol.oldParam.qProd = qProd;
oldSol.oldParam.qTrans = qTrans;
oldSol.oldParam.qSales = qSales;
oldSol.oldParam.qProd0 = qProd0;
oldSol.oldParam.bEvF = bEvF;
oldSol.oldParam.bInv = bInv;
oldSol.oldParam.bProd = bProd;
oldSol.oldParam.bTrans = bTrans;
oldSol.oldParam.bSales = bSales;
oldSol.oldParam.bProd0 = bProd0;
oldSol.oldParam.qEv = qEv;
oldSol.oldParam.qProd00 = qProd00;
oldSol.oldParam.bEv = bEv;
oldSol.oldParam.bProd00 = bProd00;
oldSol.oldParam.pruner = pruner;
oldSol.oldParam.Gbom = Gbom;

oldSol.oldTS.extraTimes = extraTimes;
oldSol.oldTS.extraTS = extraTS;
oldSol.oldTS.xlSkus = xlSkus;
oldSol.oldTS.xlTimes= xlTimes;
oldSol.oldTS.dem = dem;
oldSol.oldTS.realTS = realTS;
oldSol.oldTS.shiftFactor = shiftFactor;
oldSol.oldTS.trainForesight = trainForesight;

%Save all parameters in a struct, for organization.
%oldSol.params = sysParams;
oldSol.params.Bom = Bom;
oldSol.params.c_truck = c_truck;
oldSol.params.c_y = c_y;
oldSol.params.c_z = c_z;
oldSol.params.d_t = d_t;
oldSol.params.dem = dem;
oldSol.params.locations = locations;
oldSol.params.locNames = locNames;
oldSol.params.f_t = f_t;
oldSol.params.f_w = f_w;
oldSol.params.f_y = f_y;
oldSol.params.f_z = f_z;
oldSol.params.f_z0 = f_z0;
oldSol.params.fixedTruckCost = fixedTruckCost;
oldSol.params.hourlyTruckCost = hourlyTruckCost;
oldSol.params.mu_w = mu_w;
oldSol.params.perMileTruckCost = perMileTruckCost;
oldSol.params.R_t = R_t;
oldSol.params.R_z = R_z;
oldSol.params.s_y = s_y;


oldSol.params.nt = nt;
oldSol.params.nt0 = nt0;
oldSol.params.ntOld = ntOld;
oldSol.params.nz = nz;
oldSol.params.nloc = nloc;


filename = '~/Desktop/Theo_K./Central Intel Code/xold_12_18.mat';

save('prevSol_12_18.mat','oldSol');
