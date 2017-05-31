function pruner = Pruning(param,J)

%% Dimension determination
nz = param.nz;
nt = param.nt;
nloc = param.nloc;
nTiers = param.nTiers;
nt0 = param.nt0;


%% Pruning: Recording Indices
%TRANSPORT (cts & integer): %Note: these are just estimates, and their
%accuracy is irrelevant except to efficiency of this computationally
%cheap process.
nPrune = 2* nnz(param.zc<1)*nt + (nnz(param.tc<1)*nt + nnz(param.tc<1)*nt*nz) + nnz(param.wc==0)*nt +nnz((~any(param.c_z,1)))*nt; %number of y,z,t,w,z0,t0,z00 variables
intPrune =  nnz(param.zc<1)*nt + nnz(param.tc<1)*nt + nnz((~any(param.c_z,1)))*nt; %number of pruned z0, z00, and t0 variables
%{
If param.f_t(k,L) < 0 (set to -1 as convention), then remove columns
corresponding to each t(i,j,k,L) and t_0(j,k,L) forall i,j, which are
columns J(3,i,j,k,L) and J(6,1,j,k,L), respectively, forall i,j.
This corresponds to trucks NEVER travelling that route.
NOTE: Use tc for transportation classification
%}
pruner = zeros(nPrune,1);
count =1;
for k = 1:nloc
    for L = 1:nloc
        if (param.tc(k,L) ==0) %if transportation cannot happen between locations k and L
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
        if (param.wc(i,k) ==0) %if product i cannot be sold or bought at location k
            for j = 1:nt
                pruner(count) = J(4,i,j,k,1); %prune all unused w(i,j,k)
                count = count+1;
            end
        end
    end
end

%ACQUISITIONS: holding limit
for i = 1:nz
    for m = 1:nTiers
        for k = 1:nloc
            if (param.c_w(m,i)*param.s_y(i) >param.c_y(k)) %if minimum order size of tier m of product i is too large to store at location k
                for j = 1:nt
                    pruner(count) = J(8,i,j,k,m); %w_m(i,j,k)
                    count = count+1;
                    pruner(count) = J(9,i,j,k,m);%w_m0(i,j,k)
                    count = count +1;
                end
            end
        end
    end
end
%ACQUISITIONS: can only acquire raw materials
for i = 1:nz
    if (param.RM(i) ==0)
        for m = 1:nTiers
            for k = 1:nloc
                for j = 1:nt
                    pruner(count) = J(8,i,j,k,m); %w_m(i,j,k)
                    count = count+1;
                    pruner(count) = J(9,i,j,k,m);%w_m0(i,j,k)
                    count = count +1;
                end
            end
        end
    end
end



%PRODUCTION (cts & binary):

for k = 1:nloc
    if ((sum(param.c_z(:,k))==0))
        for j = 1:nt
            pruner(count) = J(7,1,j,k,1); %prune all unused z00(j,k) variables
            count = count+1;
        end
    end
    for i =1:nz
        if (param.zc(i,k) ==0)
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



end