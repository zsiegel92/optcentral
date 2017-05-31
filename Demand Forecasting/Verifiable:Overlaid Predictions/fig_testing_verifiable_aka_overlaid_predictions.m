symbolKey = trainForesight;%'o+*x^#w@$';
t0=datenum('5-Jan-2014');
xlTimes = xlTimes+t0;
legendTerms = {'Historical Demand'};
for i = 1:dimOutput
    legendTerms{i+1}=[num2str(trainForesight(i)),'-Timestep Prediction'];
end

plotters = [4,5,6];%Can be any SKUs (indices of skus), but including more than 4 will make the figure crowded.
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
        scatter(xlTimes(predTimeSteps(skuInds)+trainForesight(j)), realPreds(j,skuInds),symbolKey(j))
    end
    hold off
    if (i==numPlots)
        legend(legendTerms);
        legend('Location','bestoutside');
    end
end

