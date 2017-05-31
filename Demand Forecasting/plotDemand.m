function []=plotDemand(skus_to_plot,realTS,realPreds,xlTimes,realSkus)
%PRE: this function assumes that trainForesight =1:trainAhead, and that 
%nTest=1 (in forecast.m)
%POST: Plot of the demand time-series and a scatterplot of predictions, for
%the skus indicated in skus_to_plot, given in 1:nz form (not real skus)
timeInc = xlTimes(2)-xlTimes(1);
nt = size(xlTimes,2);
dimOutput = size(realPreds,1);
trainForesight = 1:dimOutput;


symbolKey = trainForesight;%'o+*x^#w@$';
t0=datenum('5-Jan-2014');
xlTimes = xlTimes+t0;
newTotTime=[xlTimes,xlTimes(nt)+timeInc*(1:max(trainForesight))];

legendTerms = {'Historical Demand'};
for i = 1:dimOutput
    legendTerms{i+1}=[num2str(trainForesight(i)),'-Timestep Prediction'];
end

for m = 1:3:length(skus_to_plot)
    plotters = skus_to_plot(m:min(m+2,length(skus_to_plot)));
    
    
    numPlots = length(plotters);
    figure
    for i = 1:numPlots
        
        subplot(numPlots,1,i)
        plot(xlTimes,realTS(plotters(i),:))
        xlabel('Date');
        ylabel(['Demand for Product ', num2str(realSkus(plotters(i)))]);
        datetick('x','keepticks','keeplimits')
        hold on
        for j = 1:dimOutput %(this displays all the predictions in trainForesight)
            scatter(newTotTime(nt+trainForesight(j)), realPreds(j,plotters(i)),symbolKey(j))
        end
        hold off
        if (i==numPlots)
            legend(legendTerms);
            legend('Location','bestoutside');
        end
    end
end

end