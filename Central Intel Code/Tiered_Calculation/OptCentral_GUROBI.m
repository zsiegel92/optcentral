close all
clear all


%% Generating Parameter Templates (call 'CSV_Parameterizer.m')
param = Manual_Parameterizer();

% load('prev_parameters.mat')
% param = params;
% clearvars -except param

%For frequent-use
nz = param.nz;
nt = param.nt;
nloc = param.nloc;
nTiers = param.nTiers;
%ENABLE THIS TO ENSURE THAT IN-PROCESS PROCESSES ARE RECORDED
%nt0 = max([max(max(param.R_t)),max(max(param.R_z))]); %Number of time-steps included as initial conditions.
nt0=param.nt0;

%% Load Demand Forecast
dem = demandGet(param);
param.dem = dem;

%% Dimension Specification
dimTotal = (4+nloc+2*nTiers)*nz*nloc*nt + (nloc^2)*nt + nt*nloc;
dimReal = (nloc+3+2*nTiers)*nz*nloc*nt;
dimInt = nz*nloc*nt*(1+nTiers) + (nloc^2)*nt + nt*nloc;
%% Generating J (Functionalized)
[J,Jinv,intcon]= Indexer(param);



%% Constraint Matrices
[qFinal,bFinal,qEqFinal,bEqFinal,qOb,LB,UB] = constraint_formalizer(param,J);

%% Pruning: calling "Pruning.m"
pruner = Pruning(param,J);

% Reparameterization

J2 = zeros(size(J)); %Index matrix
Jinv2 = Jinv;
Jinv2(pruner,:) = []; %Delete pruned indices
intcon2=zeros(dimInt,1);
count2=1; %Index of
for count = 1:size(Jinv2,1)
    J2(Jinv2(count,1),Jinv2(count,2),Jinv2(count,3),Jinv2(count,4),Jinv2(count,5)) = count;
    if (any(Jinv2(count,1)==[5,6,7,9]))
        intcon2(count2)=count;
        count2=count2+1;
    end
end

intcon2(logical(intcon2==0))=[];
clear J;
clear Jinv;
clear intcon;

% Pruning: Deleting Columns and rows
qFinal(:,pruner) = [];
qEqFinal(:,pruner) = [];
qOb(pruner) = [];
LB(pruner) = [];
UB(pruner) = [];


% bFinal(~any(qFinal,2))=[]; %Delete unused rows (unused constraints)
% qFinal(~any(qFinal,2),:) = []; %Delete unused rows (unused constraints)
% 

bFinal = bFinal(logical(sum(abs(qFinal),2)));
qFinal = qFinal(logical(sum(abs(qFinal),2)),:);
% disp('Number of zero constraint rows with nonzero constant')
% nnz(~any(qFinal,2).*any(bFinal,2)) %Should be zero

%% GUROBI!!

model = struct();
model.A = cat(1,qFinal,qEqFinal);
model.sense = [repmat('<',size(qFinal,1),1);repmat('=',size(qEqFinal,1),1)];
model.lb = LB;
model.ub = UB;
model.vtype = repmat('C',size(qEqFinal,2),1);%'B' or 'I' or 'S' or 'N' or 'C' for each variable
model.vtype(intcon2)='I';
model.vtype(intcon2(find(UB(intcon2)==1)))='B';
model.obj = full(qOb);
model.rhs = full(cat(1,bFinal,bEqFinal));
model.modelsense = 'min';

params = struct();
params.outputflag = 0;
params.resultfile = 'mip1.lp';
params.TimeLimit = 1000;


result = gurobi(model,params);
x = result.x;
% x = intlinprog(qOb,intcon2,qFinal,bFinal,qEqFinal,bEqFinal,LB,UB,options);

%% MILP

%The objective function (which is minimized) is total COST. Therefore, if
%the optimal value is negative, there is a PROFIT. If it is positive, there
%is a deficit.

%--'OutputFcn',@savemilpsolutions "puts the feasible points in a matrix named
%xIntSol in your base workspace, where each column is one integer feasible
%point. It saves the objective function values in a vector named fIntSol,"
%--'RootLPAlgorithm','primal-simplex' alternative is 'dual-simplex'
% 
% options = optimoptions('intlinprog','OutputFcn',@savemilpsolutions,'RootLPAlgorithm','primal-simplex');
% %options.MaxTime = 7200; %Default
% options.MaxTime = 15000; %4.16 hours
% %options.CutGeneration = 'advanced'; %Default is 'basic'
% options.CutGenMaxIter = 50; %Default is 10;
% x = intlinprog(qOb,intcon2,qFinal,bFinal,qEqFinal,bEqFinal,LB,UB,options);




%% Organizing output and System Statistics

[output_raw,output_tracked,stats] = SystemStats(x,J2,Jinv2,param);

%% Writing Returned Values to .csv
WriteResults(output_raw,output_tracked,stats);

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

oldSol.SOLUTION.W = output_raw.W;
oldSol.SOLUTION.Z = output_raw.Z;
oldSol.SOLUTION.T = output_raw.T;
oldSol.SOLUTION.T0 = output_raw.T0;
oldSol.SOLUTION.Z0 = output_raw.Z0;
oldSol.SOLUTION.Z00 = output_raw.Z00;
oldSol.SOLUTION.WM = output_raw.WM;
oldSol.SOLUTION.WM0 = output_raw.WM0;
oldSol.SOLUTION.Y = output_raw.Y;
oldSol.SOLUTION.track_T = output_tracked.track_T;
oldSol.SOLUTION.track_T0 = output_tracked.track_T0;
oldSol.SOLUTION.track_W = output_tracked.track_W;
oldSol.SOLUTION.track_Y = output_tracked.track_Y;
oldSol.SOLUTION.track_Z = output_tracked.track_Z;
oldSol.SOLUTION.track_Z0 = output_tracked.track_Z0;
oldSol.SOLUTION.track_Z00 = output_tracked.track_Z00;

% oldSol.oldParam.Jinv = Jinv;
% oldSol.oldParam.J = J;
% oldSol.oldParam.intcon = intcon;
% oldSol.oldParam.qEvF = qEvF;
% oldSol.oldParam.qInv = qInv;
% oldSol.oldParam.qProd = qProd;
% oldSol.oldParam.qTrans = qTrans;
% oldSol.oldParam.qSales = qSales;
% oldSol.oldParam.qAcquisitions = qAcquisitions;
% oldSol.oldParam.qProd0 = qProd0;
% oldSol.oldParam.bEvF = bEvF;
% oldSol.oldParam.bInv = bInv;
% oldSol.oldParam.bProd = bProd;
% oldSol.oldParam.bTrans = bTrans;
% oldSol.oldParam.bSales = bSales;
% oldSol.oldParam.bAcquisitions = bAcquisitions;
% oldSol.oldParam.bProd0 = bProd0;
% oldSol.oldParam.qEv = qEv;
% oldSol.oldParam.qProd00 = qProd00;
% oldSol.oldParam.bEv = bEv;
% oldSol.oldParam.bProd00 = bProd00;
% oldSol.oldParam.pruner = pruner;

% oldSol.oldTS.extraTimes = extraTimes;
% oldSol.oldTS.extraTS = extraTS;
% oldSol.oldTS.xlSkus = xlSkus;
% oldSol.oldTS.xlTimes= xlTimes;
% oldSol.oldTS.realTS = realTS;
% oldSol.oldTS.shiftFactor = shiftFactor;
% oldSol.oldTS.trainForesight = trainForesight;

%Save all parameters in a struct, for organization.
%oldSol.params = sysParams;
oldSol.params.param = param;
oldSol.param=param;
oldSol.params.Bom = param.Bom;
oldSol.params.Gbom = param.Gbom;
oldSol.params.c_truck = param.c_truck;
oldSol.params.c_y = param.c_y;
oldSol.params.c_z = param.c_z;
oldSol.params.d_t = param.d_t;
oldSol.params.dem = param.dem;
%oldSol.params.locations = locations;
oldSol.params.locNames = param.locNames;
oldSol.params.f_t = param.f_t;
oldSol.params.f_w = param.f_w;
oldSol.params.f_y = param.f_y;
oldSol.params.f_z = param.f_z;
oldSol.params.f_z0 = param.f_z0;
oldSol.params.fixedTruckCost = param.fixedTruckCost;
oldSol.params.hourlyTruckCost = param.hourlyTruckCost;
oldSol.params.mu_w = param.mu_w;
oldSol.params.perMileTruckCost = param.perMileTruckCost;
oldSol.params.R_t = param.R_t;
oldSol.params.R_z = param.R_z;
oldSol.params.s_y = param.s_y;


oldSol.params.nt = nt;
oldSol.params.nt0 = nt0;
oldSol.params.nz = nz;
oldSol.params.nloc = nloc;
oldSol.params.nTiers = nTiers;

filename = ['xold_',date, '.mat'];

save(filename);

