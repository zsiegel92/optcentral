function [J,Jinv,intcon]= Indexer(param)
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

dimTotal = (4+nloc+2*nTiers)*nz*nloc*nt + (nloc^2)*nt + nt*nloc;
%dimReal = (nloc+3+2*nTiers)*nz*nloc*nt;
dimInt = nz*nloc*nt*(1+nTiers) + (nloc^2)*nt + nt*nloc;

%% Index Matrix J
J = zeros(9, nz, nt, nloc, max(nloc,param.nTiers)); %Index matrix
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
for nmat = 8
    for i = 1:nz
        for j = 1:nt
            for k = 1:nloc
                for m = 1:nTiers
                    J(nmat,i,j,k,m)=count; %Index all wm(i,j,k) [A8]
                    Jinv(count,:)=[nmat,i,j,k,m];
                    count=count+1;
                end
            end
        end
    end
end
for nmat = 9
    for i = 1:nz
        for j = 1:nt
            for k = 1:nloc
                for m = 1:nTiers
                    J(nmat,i,j,k,m)=count; %Index all wm0(i,j,k) [A9]
                    Jinv(count,:)=[nmat,i,j,k,m];
                    intcon(count2)=count;
                    count2=count2+1;
                    count=count+1;
                end
            end
        end
    end
end
clear count; clear nmat; clear i; clear j; clear k; clear L;
end