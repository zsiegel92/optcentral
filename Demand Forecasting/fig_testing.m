function []=plotDemand(skus_to_plot)

symbolKey = trainForesight;%'o+*x^#w@$';
t0=datenum('5-Jan-2014');
xlTimes = xlTimes+t0;
newTotTime=[xlTimes,xlTimes(nt)+timeInc*(1:max(trainForesight))];

legendTerms = {'Historical Demand'};
for i = 1:dimOutput
    legendTerms{i+1}=[num2str(trainForesight(i)),'-Timestep Prediction'];
end

for m = 1:3:length(skus_to_plot)
    
plotters = [skus_to_plot(m)];%Can be any SKUs (indices of skus), but including more than 4 will make the figure crowded.
if (length(skus_to_plot)>m)
    plotters = [plotters,skus_to_plot(m+1)];
end
if (length(skus_to_plot)>m+1)
    plotters = [plotters,skus_to_plot(m+2)];
end
numPlots = length(plotters);
figure
for i = 1:numPlots
    skuInds = (1+nTest*(plotters(i)-1)):(nTest*plotters(i));
    subplot(numPlots,1,i)
    plot(xlTimes,realTS(plotters(i),:))
    xlabel('Date');
    ylabel(['Demand for Product ', num2str(realSkus(plotters(i)))]);
    datetick('x','keepticks','keeplimits')
    hold on
    for j = 1:dimOutput %(this displays all the predictions in trainForesight)
        scatter(newTotTime(predTimeSteps(skuInds)+trainForesight(j)), realPreds(j,skuInds),symbolKey(j))
    end
    hold off
    if (i==numPlots)
        legend(legendTerms);
        legend('Location','bestoutside');
    end
end

end