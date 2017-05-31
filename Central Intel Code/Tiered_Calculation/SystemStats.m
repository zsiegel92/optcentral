function [output_raw,output_tracked,stats] = SystemStats(x,J2,Jinv2,param)

%% System Variables
nz = param.nz;
nt = param.nt;
nloc = param.nloc;
nTiers = param.nTiers;

%% Storing Returned Values & System Statistics
Y = zeros(nz,nt,nloc); %(i,j,k)'th entry is quantity of i held over at loc. k at time-step j
Z = zeros(nz,nt,nloc);%(i,j,k)'th entry is quantity of i produced at loc. k at time-step j
T = zeros(nz,nt,nloc,nloc);%(i,j,k,L)'th entry is quantity of i transported from loc. k to L at time-step j
W = zeros(nz,nt,nloc);%(i,j,k)'th entry is quantity of i sold (bought if neg.) at loc. k at time-step j
Z0 = zeros(nz,nt,nloc);%(i,j,k)'th entry is WHETHER OR NOT product i went into production at location k at time-step j (0=no, 1=yes)
T0 = zeros(nt,nloc,nloc);%(i,k,L)'th entry is number of trucks was sent from loc. k to loc. L at time-step j
Z00 = zeros(nt,nloc);%(j,k)'th entry is WHETHER OR NOT location k is "in production" at time-step j
WM  = zeros(nz,nt,nloc,nTiers);
WM0 = zeros(nt,nt,nloc,nTiers);
track_Y = zeros(nnz(Y),4);%Row = [product, time-step,location,holdover]
track_Z = zeros(nnz(Z),4);%Row = [product,time-step,location,production order]
track_T = zeros(nnz(T),5);%Row = [product,time-step,location1,location2,shipment order]
track_W = zeros(nnz(W),4);%Row = [product,time-step,location,sales order]
track_Z0 = zeros(nnz(Z0),4);%Row = [product,time-step,location,"start production?" (1 or 0)]
track_T0 = zeros(nnz(T0),4);%Row = [time-step,location1,location2,number of trucks sent from 1 to 2]
track_Z00 = zeros(nnz(Z00),3);%Row = [time-step,location,"in production?" (1 or 0)]
track_WM = zeros(nnz(WM),5); %Row = [product,time-step,location,tier,quantity purchased]
track_WM0 = zeros(nnz(WM0),5);%Row = [product,time-step,location,tier,"Purchased?" (1 or 0)]

%Order of listing the components of a row of Jinv2 consisting of
%[n,i,j,k,L]
viewOrderY = [3,4,2]; %[j,k,i,y(i,j,k)]
viewOrderZ = [3,4,2]; %[j,k,i,z(i,j,k)]
viewOrderT = [3,4,5,2]; %[j,k,L,i,t(i,j,k,L)]
viewOrderW = [3,4,2]; %[j,k,i,w(i,j,k)]
viewOrderZ0 = [3,4,2]; %[j,k,i,z0(i,j,k)]
viewOrderT0 = [3,4,5]; %[j,k,L,t0(j,k,L)]
viewOrderZ00 = [3,4]; %[j,k,z00(j,k)]
viewOrderWM = [3,4,2,5];%[j,k,i,m,wm(i,j,k)]
viewOrderWM0 = [3,4,2,5];%[j,k,i,m,wm0(i,j,k)]

%Order of n: y,z,t,w,z0,t0,z00,wm,wm0
viewOrderGeneral = {viewOrderY,viewOrderZ,viewOrderT,viewOrderW,viewOrderZ0,viewOrderT0,viewOrderZ00,viewOrderWM,viewOrderWM0};

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
WMRow = (@(i) [arrayfun(@(n) Jinv2(i,n),viewOrderWM),x(i)]);
WM0Row = (@(i) [arrayfun(@(n) Jinv2(i,n),viewOrderWM0),x(i)]);

genRow = (@(i) [arrayfun(@(n) Jinv2(i,n),viewOrderGeneral(Jinv2(i,1))),x(i)]);

counts = ones(9,1); %Index of how many of each type of decision has been dealt with. [Y,Z,T,W,Z0,T0,Z00]
%NOTE: Jinv2(i,:) = [n,i,j,k,L], n maps [1,2,3,4,5]->[Y,Z,T,W,Z0,T0,Z00]
for i = 1:size(Jinv2,1)
    if (abs(x(i))>10^(-13))
        %track_Y(i,:)= [sku,time-step,location,holding]
        if (Jinv2(i,1)==1)
            Y(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4)) = x(i);
            %track_Y(counts(Jinv2(i,1)),:) = [Jinv2(i,2),Jinv2(i,3),Jinv2(i,4),x(i)];
            track_Y(counts(Jinv2(i,1)),:) = YRow(i);
            YCosts = YCosts + x(i)*param.s_y(Jinv2(i,2))*param.f_y(Jinv2(i,3),Jinv2(i,4)); %+= y(i,j,k)*param.s_y(i)*param.f_y(j,k)
            counts(Jinv2(i,1)) = counts(Jinv2(i,1))+1;
        end
        %track_Z(i,:)= [sku,time-step,location,production]
        if (Jinv2(i,1)==2)
            Z(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4)) = x(i);
            %track_Z(counts(Jinv2(i,1)),:) = [Jinv2(i,2),Jinv2(i,3),Jinv2(i,4),x(i)];
            track_Z(counts(Jinv2(i,1)),:) = ZRow(i);
            ZCosts = ZCosts + x(i)*param.f_z(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4)); %+= z(i,j,k)*param.f_z(i,j,k)
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
            salesOf(Jinv2(i,2)) = salesOf(Jinv2(i,2))+param.f_w(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4))*x(i); %+= f_w(i,j,k)*w(i,j,k)
            salesAtT(Jinv2(i,3)) = salesAtT(Jinv2(i,3)) + param.f_w(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4))*x(i);%+= f_w(i,j,k)*w(i,j,k)
            salesAtLoc(Jinv2(i,4)) = salesAtLoc(Jinv2(i,4))+ param.f_w(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4))*x(i);%+= f_w(i,j,k)*w(i,j,k)
            counts(Jinv2(i,1)) = counts(Jinv2(i,1))+1;
        end
        if (Jinv2(i,1)==5)
            Z0(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4)) = x(i);
            %track_Z0(counts(Jinv2(i,1)),:) = [Jinv2(i,2),Jinv2(i,3),Jinv2(i,4),x(i)];
            track_Z0(counts(Jinv2(i,1)),:) = Z0Row(i);
            Z0Costs = Z0Costs + param.f_z0(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4))*x(i); %+=param.f_z0(i,j,k)*z0(i,j,k)=param.f_z0(i,j,k)
            counts(Jinv2(i,1)) = counts(Jinv2(i,1))+1;
        end
        if (Jinv2(i,1)==6)
            T0(Jinv2(i,3),Jinv2(i,4),Jinv2(i,5)) = x(i);
            %track_T0(counts(Jinv2(i,1)),:) = [Jinv2(i,3),Jinv2(i,4),Jinv2(i,5),x(i)];
            track_T0(counts(Jinv2(i,1)),:) = T0Row(i);
            TCosts = TCosts + x(i)*param.f_t(Jinv2(i,4),Jinv2(i,5)); %+= T0(j,k,L)*param.f_t(k,L)
            counts(Jinv2(i,1)) = counts(Jinv2(i,1))+1;
        end
        if (Jinv2(i,1)==7)
            Z00(Jinv2(i,3),Jinv2(i,4)) = x(i);
            %track_Z00(counts(Jinv2(i,1)),:) = [Jinv2(i,3),Jinv2(i,4),x(i)];%[i,k,Z00(j,k)]
            track_Z00(counts(Jinv2(i,1)),:) = Z00Row(i);
            counts(Jinv2(i,1)) = counts(Jinv2(i,1))+1;
        end
        if (Jinv2(i,1)==8)
            WM(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4),Jinv2(i,5)) = x(i);
            %track_WM(counts(Jinv2(i,1)),:) = [Jinv2(i,2),Jinv2(i,3),Jinv2(i,4),Jinv2(i,5),x(i)];
            track_WM(counts(Jinv2(i,1)),:) = WMRow(i);
            counts(Jinv2(i,1)) = counts(Jinv2(i,1))+1;
        end
        if (Jinv2(i,1)==9)
            WM0(Jinv2(i,2),Jinv2(i,3),Jinv2(i,4),Jinv2(i,5)) = x(i);
            %track_T0(counts(Jinv2(i,1)),:) = [Jinv2(i,3),Jinv2(i,4),Jinv2(i,5),x(i)];
            track_WM0(counts(Jinv2(i,1)),:) = WM0Row(i);
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
track_WM = sortrows(track_WM);
track_WM0 = sortrows(track_WM0);

%% Writing Outputs

output_raw.Y = Y;
output_raw.Z = Z;
output_raw.T = T;
output_raw.W = W;
output_raw.Z0 = Z0;
output_raw.T0 = T0;
output_raw.Z00 = Z00;
output_raw.WM = WM;
output_raw.WM0=WM0;
output_tracked.track_Y=track_Y;
output_tracked.track_Z=track_Z;
output_tracked.track_T=track_T;
output_tracked.track_W=track_W;
output_tracked.track_Z0=track_Z0;
output_tracked.track_T0=track_T0;
output_tracked.track_Z00=track_Z00;
output_tracked.track_WM=track_WM;
output_tracked.track_WM0=track_WM0;

stats.salesOf = salesOf;
stats.unitsSold = unitsSold;
stats.salesAtLoc=salesAtLoc;
stats.salesAtT=salesAtT;
stats.TCosts=TCosts;
stats.ZCosts=ZCosts;
stats.Z0Costs=Z0Costs;
stats.YCosts=YCosts;


end