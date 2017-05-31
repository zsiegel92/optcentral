function dem = demandGet(param)
%[net,extraTS,extraTimes,xlTimes,realTS,xlSkus,trainForesight]=forecast();
load('forecast_workspace_12-Jan-2017.mat')
clearvars -except net extraTS extraTimes xlTimes realTS xlSkus trainForesight param


%% Demand Processing: Cutting Data per model size
extraTS = extraTS(:,1:(param.nt-param.nt0));
shiftFactor = abs(min(min(min(extraTS)),min(min(realTS))));
extraTS = extraTS + shiftFactor;%Imperfect solution to negative demand predictions
realTS = realTS + shiftFactor;%Imperfect solution to negative demand predictions
extraTS = extraTS(1:param.nz,:);%Optimize based on nz products only.
realTS = realTS(1:param.nz,:);
xlSkus = xlSkus(1:param.nz);

ntOld = size(realTS,2);


dem = zeros(param.nz,param.nt,param.nloc);
for i = 1:param.nz
    divFactor = 1/param.nloc;
    for j = 1:param.nt0
        for k = 1:param.nloc
            dem(i,j,k) = realTS(i,ntOld-param.nt0+j)*divFactor*(param.wc(i,k)>0);
        end
    end
    for j = (param.nt0+1):param.nt
        for k = 1:param.nloc
            dem(i,j,k) =  extraTS(i,j-param.nt0)*divFactor*(param.wc(i,k)>0);
        end
    end
end

end