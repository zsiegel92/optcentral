
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
oldSol.params.ntOld = ntOld;
oldSol.params.nz = nz;
oldSol.params.nloc = nloc;
oldSol.params.nTiers = nTiers;

filename = '~/Desktop/Theo_K./Central Intel Code/xold_1_9.mat';

save('prevSol_1_9.mat','oldSol');

