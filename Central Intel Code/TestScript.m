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

filename = '/Users/Zach/Desktop/Theo_K./Central Intel Code/Tiered_Calculation/xold_1_10.mat';

save(filename);
