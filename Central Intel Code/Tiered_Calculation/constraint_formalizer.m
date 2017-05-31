function [qFinal,bFinal,qEqFinal,bEqFinal,qOb,LB,UB] = constraint_formalizer(param,J)
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
J(8,i,j,k,m)= wm(i,j,k)
J(9,i,j,k,m)=wm0(i,j,k)
(io means "index of")
%}


%% Dimension determination
nz = param.nz;
nt = param.nt;
nloc = param.nloc;
nTiers = param.nTiers;
nt0 = param.nt0;

dimTotal = (4+nloc+2*nTiers)*nz*nloc*nt + (nloc^2)*nt + nt*nloc;
%dimReal = (nloc+3+2*nTiers)*nz*nloc*nt;
dimInt = nz*nloc*nt*(1+nTiers) + (nloc^2)*nt + nt*nloc;
%% Initializing constraint matrices

qFinal = sparse(0,dimTotal);
bFinal = sparse(0,1);
qEqFinal = sparse(0,dimTotal);
bEqFinal = sparse(0,1);
%% Objective Function: qOb*x
%Note that this objective function, which maximizes profit (sales revenue-costs)
%DOES NOT maximize the value of the inventory that is IN TRANSIT or BEING
%HELD. This is an extreme simplification that can be accounted for easily.

qOb = sparse(1,dimTotal);
for i = 1:nz
    for j = 1: nt
        for k = 1:nloc
            qOb(J(1,i,j,k,1)) = -param.f_y(j,k)*param.s_y(i); %holding cost
            qOb(J(2,i,j,k,1)) = -param.f_z(i,j,k); %proportional price of production
            qOb(J(4,i,j,k,1)) = param.f_ws(i,j,k); %sales revenue
            qOb(J(5,i,j,k,1)) = -param.f_z0(i,j,k);%fixed cost of production
            for m = 1:nTiers
                qOb(J(8,i,j,k,m)) = -param.f_w(i,j,k,m);
            end
            for L = 1:nloc
                qOb(J(6,1,j,k,L)) = -param.f_t(k,L); %price of sending trucks
            end
        end
    end
end
qOb = -qOb; %intlinprog MINIMIZES objective; qOb was formulated for maximization.
%% Product Equity (Added to objective function)
qEquity = sparse(1,dimTotal);

qOb = qOb - qEquity;
%% Evolution Equation: qEv*x == 0 forall i, j, k
qEv = sparse(nz*(nt-nt0)*nloc,dimTotal);
bEv = sparse(nz*(nt-nt0)*nloc,1); %Equality constraint with no constant term
count = 1;
for i = 1:nz
    for j = (nt0+1):nt %INITIALIZE AT j=nt0 !!!!
        for k = 1:nloc
            q = zeros(1,dimTotal);
            q(J(1,i,(j-1),k,1)) = 1; %holdover
            q(J(4,i,(j-1),k,1)) = -1; %sales
            for m = 1:nTiers
                q(J(8,i,(j-1),k,m))=1;%acquisitions
                
            end
            if (j>param.R_z(i,k)) %if product can come out of production at this time-step
                q(J(2,i,(j-param.R_z(i,k)),k,1)) = 1; %product i @ location k coming out of production at time j, which started at time-step j - param.R_z(i,k) by defparam.inition
            end
            for I = [1:(i-1),(i+1):nz];%Production of an item never subtracts from its own inventory
                q(J(2,I,j-1,k,1)) = -param.Bom(I,i); %product i @ location k being used for production of product I at time-step (j-1)
            end
            for I = 1:nloc
                if (j>param.R_t(I,k))
                    q(J(3,i,j-param.R_t(I,k),I,k)) = 1; %transported inventory: incoming
                end
                q(J(3,i,j-1,k,I)) = -1; %transported inventory: outgoing
            end
            %Shipping to self is forbidden, so it's not necessary to
            %carefully constrain self-shipping-driven inventory changes.
            q(J(3,i,j-1,k,k))=0;%Very strange fix. Not sure if this is correct.
            
            q(J(1,i,j,k,1)) = -1; %Set equality for current (i,j,k) in equation
            qEv(count,:) = q;
            count = count +1;
        end
    end
end
clear count; clear I; clear q;

qEqFinal = cat(1,qEqFinal,qEv);
bEqFinal = cat(1,bEqFinal,bEv);
%% Evolution FINAL-TIMESTEP Equation q*x <= 0
%Note: No acquisitions in final time-step

%w(i,nt,k) + Sum_I z(I,nt,k)param.Bom(I,i) + Sum_L t(i,nt,k,L) + w(i,nt,k) -
%y(i,nt,k) <= 0
qEvF = sparse(nz*nloc,dimTotal);
bEvF = sparse(nz*nloc,1);
count = 1;
for i = 1:nz
    for k = 1:nloc
        q = zeros(1,dimTotal);
        q(J(1,i,nt,k,1)) = -1; %y(i,nt,k)
        q(J(4,i,nt,k,1))=1; %w(i,nt,k)
                            %Note: not necessary to control final time-step
                            %acquisitions in this way.
        for I = 1:nz
            q(J(2,I,nt,k,1))=param.Bom(I,i);  %param.Bom(I,i)*z(I,nt,k)
        end
        for I = 1:nloc
            q(J(3,i,nt,k,I))=1; %Sum_L t(i,nt,k,I)
        end
        qEvF(count,:) = q;
        count = count +1;
    end
end
clear count; clear i; clear k; clear q; clear I;

qFinal = cat(1,qFinal,qEvF);
bFinal = cat(1,bFinal,bEvF);
%% Inventory/Holding: Sum_i s_y(i)*y(i,j,k) <= param.c_y(k) forall j,k
qInv = sparse(nt*nloc,dimTotal);
bInv = sparse(nt*nloc,1); %Less than constraint
count = 1;
for j = 1:nt
    for k = 1:nloc
        q = zeros(1,dimTotal);
        for i = 1:nz
            q(J(1,i,j,k,1))=param.s_y(i); %Each item held at a location at a time step y(i,j,k) contributes to "total units held", which much be less than holding capacity at that location param.c_y(k)
        end
        qInv(count,:) = q;
        bInv(count) = param.c_y(k);
        count = count+1;
    end
end

clear count; clear i; clear j; clear k; clear q;

qFinal = cat(1,qFinal,qInv);
bFinal = cat(1,bFinal,bInv);
%% Production: z(i,j,k) <- param.c_z(i,k)z_0(i,j,k) forall i,j,k
qProd = sparse(nz*nt*nloc,dimTotal);
bProd = sparse(nz*nt*nloc,1); %Less than constraint

count = 1;
for i = 1:nz
    for j = 1:nt
        for k = 1:nloc
            q = zeros(1,dimTotal);
            q(J(2,i,j,k,1))=1; %z(i,j,k)
            q(J(5,i,j,k,1)) = -param.c_z(i,k); %-param.c_z(i,k)z_0(i,j,k)
            qProd(count,:) = q;
            %bProd(count) = 0;
            count = count+1;
        end
    end
end

clear i; clear j; clear k; clear L; clear q; clear count;

qFinal = cat(1,qFinal,qProd);
bFinal = cat(1,bFinal,bProd);
%% Transportation: Sum_i t(i,j,k,L)*param.s_y(i) <= param.c_truck*t_0(j,k,L) forall j,k,L
qTrans = sparse(nt*nloc*nloc,dimTotal);
bTrans = sparse(nt*nloc*nloc,1); %Less than constraint
count = 1;
for j = 1:nt
    for k = 1:nloc
        for L=1:nloc
            q = zeros(1,dimTotal);
            for i =1:nz
                q(J(3,i,j,k,L)) = param.s_y(i);
            end
            q(J(6,1,j,k,L)) = -param.c_truck;
            qTrans(count,:) = q;
            count = count +1;
        end
    end
end


clear i; clear j; clear k; clear L; clear q; clear count;

qFinal = cat(1,qFinal,qTrans);
bFinal = cat(1,bFinal,bTrans);
%% Sales qSales*x <= bSales
%Sum_{T=1}^j param.mu_w(i,k)^(j-T)w(i,T,k)  <= Sum_{T=1}^j param.mu_w(i,k)^(j-T)dem(i,T,k) forall i,j,k

qSales = sparse(nz*nt*nloc,dimTotal);
bSales = sparse(nz*nt*nloc,1); %Less than constraint

count = 1;
for i = 1:nz
    for j = 1:nt
        for k = 1:nloc
            q = zeros(1,dimTotal);
            for T =1:j %This accounts for BACKLOG
                q(J(4,i,T,k,1))=param.mu_w(i,k)^(j-T); %w(i,j,k) %NOTE that 0^0 = 1 in MATLAB.
                bSales(count) = bSales(count) + param.mu_w(i,k)^(j-T)*param.dem(i,T,k); %param.mu_w(i,k)*dem(i,T,k)
            end
            qSales(count,:) = q;
            %bProd(count) = 0;
            count = count+1;
        end
    end
end

clear i; clear j; clear k; clear L; clear q; clear count;

qFinal = cat(1,qFinal,qSales);
bFinal = cat(1,bFinal,bSales);
%% Tiered Acquisitions qAcquisitions*x <= bAcquisitions
%c_wm(i)*w_m0(i,j,k) <= w_m(i,j,k)
%wm(i,j,k) < = w_m0(i,j,k)*c_y(k)/s_y(i) %Can't purchase anything if
%un-tiered, and if tiered, can only purchase maximum holding capacity.
qAcquisitions = sparse(2*nTiers*nz*nt*nloc,dimTotal);
bAcquisitions = sparse(2*nTiers*nz*nt*nloc,1); %Less than constraint

count = 1;
for i = 1:nz
    for j = 1:nt
        for k = 1:nloc
            for m = 1:nTiers
                %c_wm(i)*w_m0(i,j,k) <= w_m(i,j,k)
                q = zeros(1,dimTotal);
                q(J(8,i,j,k,m)) = -1; %-w_m(i,j,k)
                q(J(9,i,j,k,m))=param.c_w(m,i); %c_w(m,i)w_m0(i,j,k)
                %bAcquisitions(count) = 0; %q*x<= 0
                qAcquisitions(count,:) = q;
                %bProd(count) = 0;
                count = count+1;
                
                %wm(i,j,k) < = w_m0(i,j,k)*c_y(k)/s_y(i)
                q = zeros(1,dimTotal);
                q(J(8,i,j,k,m)) = 1; %w_m(i,j,k)
                q(J(9,i,j,k,m))=-param.c_y(k)/param.s_y(i); %-w_m0(i,j,k)
                %bAcquisitions(count) = 0; %q*x<= 0
                qAcquisitions(count,:) = q;
                %bProd(count) = 0;
                count = count+1;
            end
        end
    end
end

clear i; clear j; clear k; clear L; clear q; clear count;

qFinal = cat(1,qFinal,qAcquisitions);
bFinal = cat(1,bFinal,bAcquisitions);
%% Production DISCRETE z_00(j,k)+ Sum_i z_0(i,j,k)  <= 1 forall j,k
qProd0 = sparse(nt*nloc,dimTotal);
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

qFinal = cat(1,qFinal,qProd0);
bFinal = cat(1,bFinal,bProd0);
%% Production Discrete2 -z00(j,k) + Sum_i Sum_(1<=TT<=param.R_z(i,k)) z0(i,j-TT,k) == 0 forall j,k
% This constraint ensures that only one product can be "in production" at a
% time
qProd00 = sparse(nt*nloc,dimTotal);
bProd00 = sparse(nt*nloc,1); %Less than constraint
count = 1;
for j = 1:nt
    for k = 1:nloc
        q = zeros(1,dimTotal);
        q(J(7,1,j,k,1))=-1; %z00(j,k)
        for i = 1:nz
            for TT = 1:(param.R_z(i,k)-1) %%NOTE: sum runs to param.R_z(i,k)-1 instead of param.R_z(i,k) because
                %j-param.R_z(i,k) is the time-step of the production order, not an "in-production" time-step
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

qEqFinal = cat(1,qEqFinal,qProd00);
bEqFinal = cat(1,bEqFinal,bProd00);
%% Transportation DISCRETE t(i,j,k,k)=0 forall i,j,k && t_0(j,k,k)=0 forall j,k
%See upper bound
%% Hard-Constrained Linear Decisions
%Production of some FG happens throughout model duration

%J(5,i,j,k,1) = index of z0(i,j,k) = (1 or 0) whether i is produced at location k at time j
% 
% numConstraints = 1;
% qLinearDecision = sparse(numConstraints,dimTotal);
% bLinearDecision = sparse(numConstraints,1); %Less than constraint
% count = 1;
% q = zeros(1,dimTotal);
% for i = 1:nz%find(param.FG)
%     for j = 1:nt
%         for k = 1:nloc
%             q(J(2,i,j,k,1))=-1; %if i = 1:nz then SOMETHING is produced.
%             %if i = find(param.FB):
%             %Sum of all z variables corresponding to FG skus is >=1 (so sum of negatives is <= -1)
%             %This enforces that AT LEAST ONE burger or
%             %egg sandwich is produced throughout
%             %production.
%             %q(1,J(4,i,j,k,1)) = -1; %Sum of all w variables is >=1. SOMETHING must be sold.
%         end
%     end
% end
% bLinearDecision(count,1) = -1;
% qLinearDecision(count,:) = q;
% count = count+1;
% 
% 
% clear i; clear j; clear k; clear L; clear q; clear count; clear numConstraints;
% 
% qFinal = cat(1,qFinal,qLinearDecision);
% bFinal = cat(1,bFinal,bLinearDecision);

%% Bounds:
%Lower: >=0 for all (except <=0 for acquisitions - this requires NO constraint!)
%Upper Bounds: variable classification (and demand)


LB = zeros(dimTotal,1);
UB = Inf*ones(dimTotal,1);

%SALES: Can't sell beyond demand!
for i = 1:nz
    for k = 1:nloc
        if (param.wc(i,k)==1) %Can sell product i at location k
            %Sales default maximum is Inf given demand backlog constraints.
            for j = 1:nt
                UB(J(4,i,j,k,1)) = 0;
                for T = 1:j
                    UB(J(4,i,j,k,1)) = UB(J(4,i,j,k,1))+ param.mu_w(i,k)^(j-T)*param.dem(i,T,k); %Demand Backlog constraint in param.SALES equation.
                    %UB(J(4,i,j,k,1)) = dem(i,j,k); %w(i,j,k) <= dem(i,j,k)
                end
                UB(J(4,i,j,k,1)) = min(UB(J(4,i,j,k,1)),param.c_y(k)/param.s_y(i));
            end
        else %param.wc(i,k)==0 %Can't sell product i at location k
            for j = 1:nt
%                 LB(J(4,i,j,k,1))=0; %Already = 0
                UB(J(4,i,j,k,1))=0; %Can't sell at all if classified as non-retail
            end
        end
    end
end

%ACQUISITIONS:Can't purchase more than you can hold
%Binary tiered-purchasing decision variables are binary-constrained
for i = 1:nz
    for k = 1:nloc
        for j = 1:nt
            for m = 1:nTiers
                UB(J(8,i,j,k,m))=param.c_y(k)/param.s_y(i); %Can't purchase more than can be held
                UB(J(9,i,j,k,m))=1; %wm0(i,j,k)<=1
                if (param.f_w(i,j,k,m)<0)
                    UB(J(8,i,j,k,m))=0;
                    UB(J(9,i,j,k,m))=0;
                end
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
%z(i,j,k)<=z0(i,j,k)param.c_z(i,k),but may increase efficiency by reducing
%feasible region.
for i=1:nz
    for j = 1:nt
        for k = 1:nloc
            UB(J(2,i,j,k,1))=param.c_z(i,k); %z(i,j,k)<=param.c_z(i,k)
        end
    end
end

%Holding: Can't hold more than capacity
%NOTE: This is accounted for in the INVENTORY/HOLDING constraint
%Sum_i y(i,j,k) <= param.c_y(k) forall j,k
%since y(i,j,k)>=0 forall i,j,k
for i = 1:nz
    for j = 1:nt
        for k = 1:nloc
            UB(J(1,i,j,k,1))=param.c_y(k); %y(i,j,k)<=param.c_y(k)
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
%% Initialization: hard-constraining decisions, and vetoes
% Let init be a matrix the same size as J (i.e., accommodates all decision
% variables) consisting of all zeros EXCEPT for entries init(n,i,j,k,L)
% where j<=nt0, which means the time-step is included in initialization. In
% that case, init(n,i,j,k,L) contains the initialized value of the decision
% variable being referenced.
for j = 1:nt0
    for k = 1:nloc
        for n = [1,2,4,5] %Initialize y,z,t,w,z0,t0,z00 (i,j,k)
            for i = 1:nz
                UB(J(n,i,j,k,1))= param.init(n,i,j,k,1);
                LB(J(n,i,j,k,1))= param.init(n,i,j,k,1);
            end
        end
        for n = [8,9] %Initialize wm, wm0 (i,j,k)
            for i = 1:nz
                for m = 1:nTiers
                    UB(J(n,i,j,k,m))= param.init(n,i,j,k,m);
                    LB(J(n,i,j,k,m))= param.init(n,i,j,k,m);
                end
            end
        end
        for n = 3 %Initialize t(i,j,k,L)
            for i = 1:nz
                for L = 1:nloc
                    UB(J(n,i,j,k,L))= param.init(n,i,j,k,L);
                    LB(J(n,i,j,k,L))= param.init(n,i,j,k,L);
                end
            end
        end
        for n = 6 %Initialize t0(j,k,L)
            for L = 1:nloc
                UB(J(n,1,j,k,L))= param.init(n,1,j,k,L);
                LB(J(n,1,j,k,L))= param.init(n,1,j,k,L);
            end
        end
        for n = 7 %Initialize z00(j,k)
            UB(J(n,1,j,k,1))=param.init(n,1,j,k,1);
            LB(J(n,1,j,k,1))=param.init(n,1,j,k,1);
        end
    end
end
end