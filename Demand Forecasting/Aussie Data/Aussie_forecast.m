clear all;
close all;
symbolKey = 'o+*x';

%% Importing Data
%[xlTS,xlSkus,xlTimes] = import_spreadsheet; %import_spreadsheet currently returns one product's time-series

[xlTS,xlSkus,xlTimes] = Aussie_data_test;
disp('Data imported');
%% Data Constants
nz = size(xlTS,1); %Number of products for which demand is recorded
nloc = 10; %Number of locations at which demand is recorded
%nt = 10000; %Number of time-steps at which demand is recorded
nt = size(xlTS,2);
%NOTE: at first, until complexity is increased, it will be assumed that
%each product is sold at each location.

%delay_line=[0,1,2,3,4,5,10,15,25];
delay_line = [0,1,2,10,25];
nDelay = length(delay_line);
delay_constant = max(delay_line);

nInfo = nDelay + 3; %Time-step, number of delay terms, first non-zero time-step, and (for now) sku

nFactors = 20;%Convolutional memory terms
dimInput = nFactors + nInfo; %Number of convolutional training terms, plus nDelay+1 for time-step and delayed demand, plus sku as inputs
trainForesight = [1,10,20,30,40]; %Training on demand(i+1), demand(i+10), demand(i+50), demand(i+100)
dimOutput = length(trainForesight);


maxTrain = nt - max(trainForesight); %The highest time-step at which training
%can happen. This is because for a training input, there needs to be an expected
%output consisting of future demand. If a late time-step is used, there IS
%no future demand. 

%%%%%NN Constants
nTest =20;
nTrain = floor((maxTrain-delay_constant-nTest)); %Number of training samples

%% Normalizing Input
productMeans = mean(xlTS,2);
xlTS = xlTS - repmat(productMeans,1,nt);
productStds = ((1/nt)*sum(xlTS.^2,2)).^((1/2));
xlTS = xlTS./repmat(productStds,1,nt);
%% Synthesizing Demand (OPTIONAL)
% ambientDim = 12; %Number of bases to choose from
% implicitDim = 6; %Each demand time-series consists a linear combination of implicitDim bases
% demand = zeros(nz,nloc,nt);
% bases = zeros(1,ambientDim,nt);
% for n = 1:5
%     %bases(1,n,1:nt) = reshape(sin(.0001*pi*n*(1:nt)),1,1,nt); %Frequencies of the form n/pi, n an integer
%     bases(1,n,1:nt) = reshape(sin((1/nt)*2*pi*(5/n)*(1:nt)),1,1,nt); %Frequencies of the form n/pi, n an integer
% end
% for n = 6:10
%     bases(1,n,1:nt) = reshape((1/30)*2*n.*log(1:nt),1,1,nt); %weird scaled polynomials
% end
% for n =11:ambientDim
%     bases(1,n,1:nt)=-reshape((1/30)*50*(1.02.^((-n/25).*(1:nt))),1,1,nt);
%     bases(1,n,1:nt)=bases(1,n,1:nt)+abs(min(bases(1,n,1:nt)));
% end
% %rando = randi([1,ambientDim],implicitDim,nz,nloc);
% rando1 = randi([1,5],3,nz,nloc);
% rando2 = randi([6,10],2,nz,nloc);
% rando3 = randi([11,ambientDim],1,nz,nloc);
% rando = [rando1;rando2;rando3]; %Each demand vector is a linear combination of 3 category-1 bases (sinusoids),
% %2 category-2 bases (logs), and 1 category-3
% %basis (exponential).
% for i = 1:nz
%     for k = 1:nloc
%         randSamp = rando(:,i,k);
%         %demand(i,k,:) = sum(reshape(reshape(bases(1,randSamp,:),nt,implicitDim)*rand([implicitDim,1]),1,1,nt),2);
%         demand(i,k,:) = sum(bases(1,randSamp,:),2);
%     end
% end
%
% timeSeries0 = squeeze(demand(1,1,:)); %For now, only train on one demand time-series, becomes nt x 1

%% Importing Demand
timeSeries0 = xlTS;

%% Optimizing Parameters
 nExp =5; %# Exponential Memory Terms
 nSin =3; %# Sinusoidal Memory Terms
 nWav = 7;%# Wavelet Memory Terms
 nGam = 5;%# Gamma Memory Terms

%Exponential Parameters
nnominees = floor(nExp*1.5);
takesome = @(x) x(1:nnominees);
candidates = .5:.5:100;
candidates = arrayfun(@(x) 1-(x/nt),candidates);

expNoms =  zeros(nz,nnominees);
for i=1:nz
    expNoms(i,:) = takesome(exp_param_output(timeSeries0(i,:),maxTrain,candidates));
end
[~,expParams] = kmeans(expNoms(:),nExp);
clear nnominees; clear takesome; clear candidates; clear expNoms;


% Sinusoidal Parameters
nnominees = nSin*2;
takesome = @(x) x(1:nnominees);
sinfreqs = zeros(nz,nnominees);
candidates = 1:.5:200; %Candidate period parameters
for i = 1:nz
    sinfreqs(i,:) = takesome(sine_freq_output(timeSeries0(i,:),maxTrain,candidates));
end
[~,sinParams] = kmeans(sinfreqs(:),nSin);
clear sinfreqs; clear takesome; clear candidates; clear nnominees;



% Wavelet Parameters
nnominees = nWav*2;
takesome = @(x) x(1:nnominees);
wavScales = zeros(nz,nnominees);
candidates = 1:.5:200; %Candidate period parameters
for i = 1:nz
    wavScales(i,:) = takesome(wavelet_param_output(timeSeries0(i,:),maxTrain,candidates));
end
[~,wavParams] = kmeans(wavScales(:),nWav);
clear wavScales; clear takesome; clear candidates; clear nnominees;

% Gamma Parameters
% nnominees = nGam*2;
% takesome = @(x) x(:,1:nnominees);
% gamParams = zeros(nz,nnominees);
% candidates1 = 1:.5:200; %Candidate period parameters
% candidates2 = [1,10];
% for i = 1:nz
%     gpo = gamma_param_output(timeSeries0(i,:),maxTrain,candidates1,candidates2);
%     gamParams(i,:) = takesome(gpo);
% end
% [~,gamParams] = kmeans(gamParams(:),nGam);
gamParams = [.05,.1,.5,.7,.9];
gamParams = [gamParams;ones(1,size(gamParams,2))];
nGam = size(gamParams,2);
clear takesome; clear candidates; clear nnominees;



%% Specifying Convolutional Bases

memTerms = zeros(nFactors,nt);

count = 1;

%Exponential Memory Terms
params = expParams;
j = 1;
for count = 1:nExp
    memTerms(count,:) = expTerm(nt,params(j));
    j = j +1;
end

%Sinusoidal Memory Terms
params = sinParams;%Period
j=1;
for count = count+1:count+nSin
    memTerms(count,:) = sin(2*pi*(1/params(j))*(1:nt));
    j=j+1;
end
%Wavelet Memory Terms
params = wavParams';
params = [params;zeros(1,size(params,2))];
j=1;
for count = count+1:count+nWav
    memTerms(count,:) = wavelet(nt,params(1,j),params(2,j));
    j=j+1;
end

%Gamma Memory Terms
params = gamParams;
j=1;
for count=count+1:count+nGam
    memTerms(count,:) = gammaTerm(nt,params(2,j),params(1,j));
    j=j+1;
end
clear count;


%% Generating Training data from Demand
%maxConvDepth = floor(maxTrain/3);%Only convolve over, at most, 1/3 of the
%historical time-steps
maxConvDepth=maxTrain; %Maximum
%memTerms = memTerms(:,1:maxTrain); %Only need memory terms on maxTrain time-steps
memTerms = memTerms(:,1:maxConvDepth); %Only need memory terms on maxConvDepth time-steps
trainSamples = ((delay_constant+1):(maxTrain-nTest))';%All possible trainable samples
fullTrainData = zeros(nTrain*nz,dimInput);
trainData = zeros(nTrain,dimInput); %Each sample input is a row

for i = 1:nz
    trainData(:,1) = trainSamples(:);%time-step
    timesWithDelays =repmat(trainSamples(:),1,nDelay)-repmat(delay_line,nTrain,1);
    trainData(:,2:(2+nDelay-1)) = arrayfun(@(n) timeSeries0(i,n),timesWithDelays);
    trainData(:,nInfo-1) = find(timeSeries0(i,:),1); %First time-step with non-zero demand
    trainData(:,nInfo) = xlSkus(i)*ones(nTrain,1);%sku
    
    %Center and normalize timeSeries
    timeSeries = timeSeries0(i,:) - mean(timeSeries0(i,:));
    timeSeries = timeSeries/std(timeSeries);
    %First two components of trainData samples are 1) time-step and 2) demand
    for j = 1:nFactors
        trainData(:,j+nInfo)=convolve_at_scaled(timeSeries,memTerms(j,:),trainSamples(:),maxConvDepth);
    end 
    fullTrainData(((i-1)*nTrain+1):(i*nTrain),:)=trainData;
end

%IF SETTING ASIDE TEST SAMPLES AT END


testSamples = ((maxTrain-nTest+1):maxTrain)'; %If setting aside test samples at end

testDataPrism = zeros(nz,dimInput,nTest);%If setting aside test samples


for i = 1:nz
    testDataPrism(i,1,:)= testSamples(:); %If setting aside test samples
    testTimesWithDelays = (repmat(testSamples,1,nDelay)-repmat(delay_line,nTest,1))';
    testDataPrism(i,2:(2+nDelay-1),:) = arrayfun(@(n) timeSeries0(i,n),testTimesWithDelays);
    testDataPrism(i,nInfo-1,:)=find(timeSeries0(i,:),1);%First time-step with non-zero demand
    testDataPrism(i,nInfo,:) = xlSkus(i)*ones(1,nTest);%sku
    
    for j = 1:nFactors
        testDataPrism(i,j+nInfo,:)= convolve_at_scaled(timeSeries,memTerms(j,:),testSamples,maxConvDepth);
    end
end

%% Specifying Desired Output for Training
%trainDesired = zeros(nTrain,dimOutput); %Each sample's desired output is a row

%trainDesired = repmat(timeSeries0(trainingSamples),1,nDelay)+repmat(delay_line,nTrain,1);
fullTrainDesired = zeros(nTrain*nz,dimOutput);
for i = 1:nz
    timesWithSkips = repmat(trainSamples(:),1,dimOutput)+repmat(trainForesight,nTrain,1);
    trainDesired = arrayfun(@(n) timeSeries0(i,n),timesWithSkips);
    %trainDesired = timeSeries0(i,repmat(trainSamples(:)',1,dimOutput)+repmat(trainForesight,nTrain,1));
    fullTrainDesired((i-1)*nTrain+1:i*nTrain,:)=trainDesired;
end


fullTestDesired = zeros(nz,dimOutput,nTest);
timesWithSkips = repmat(testSamples',dimOutput,1) +  repmat(trainForesight',1,nTest);
for i = 1:nz
    fullTestDesired(i,:,:) = arrayfun(@(n) timeSeries0(i,n),timesWithSkips);
end





%% Training the NN
trainData = fullTrainData';
trainDesired = fullTrainDesired';
%testDataPrism = reshape(testDataPrism,nz,dimInput, nTest); %If setting aside test data at end of time series
%fullTestDesired = reshape(fullTestDesired,nz,dimOutput,nTest); %If setting aside test data at end of time series

testDataStrip = zeros(dimInput,nz*nTest);
testDesiredStrip = zeros(dimOutput,nz*nTest);
for i = 1:nz
    testDataStrip(:,((i-1)*nTest+1):(i*nTest))=squeeze(testDataPrism(i,:,:));
    testDesiredStrip(:,((i-1)*nTest+1):(i*nTest)) = squeeze(fullTestDesired(i,:,:));
end

net = feedforwardnet([25,25]);%one hidden layer of 20 nodes
net.trainParam.max_fail = 6; %Default is 6!
net.divideFcn = 'divideind';
trainRatio = .7;
valRatio = .3;
%valRatio = .15;
%testRatio = .15;

net.divideParam.trainInd = sort(randperm(size(trainData,2),floor(size(trainData,2)*trainRatio)));
net.divideParam.valInd =setdiff(1:size(trainData,2),net.divideParam.trainInd);
net.divideParam.testInd = size(trainData,2)+1:size(trainData,2)+size(testDataStrip,2);


trainData = [trainData,testDataStrip];
trainDesired = [trainDesired,testDesiredStrip];
net = train(net,trainData,trainDesired);
netOutput = sim(net,trainData);

% Testing Net on Each Individual Product's Testing (latter) Time-Steps
% testOutput = zeros(nz,dimOutput,nTest);
% for i = 1:nz
%     testOutput(i,:,:) = reshape(sim(net,squeeze(testDataPrism(i,:,:))),1,dimOutput,nTest);
%     testPerformi = perform(net,squeeze(testOutput(i,:,:)),squeeze(fullTestDesired(i,:,:)));
%     disp(['Performance on test data, product ', int2str(i), ' is: ', int2str(testPerformi)]);
% end


%% Specifying Forecast

%Challenge: Include time-steps (ts-maxTrain-nTest):

predIndices = net.divideParam.testInd;
predTimeSteps = trainData(1,predIndices);
preds = netOutput(:,predIndices);
predSkus = trainData(nInfo,net.divideParam.testInd);


realTS = timeSeries0.*repmat(productStds,1,size(timeSeries0,2));
realTS = realTS + repmat(productMeans,1,size(realTS,2));

realPreds = preds.*repmat(productStds(predSkus)',dimOutput,1);
realPreds = realPreds + repmat(productMeans(predSkus)',dimOutput,1);

symbolKey = 'o+*x^#w@$';

t0=datenum('1-Jan-1956');
xlTimes = xlTimes+t0;
prodNames = {'Steel', 'Coal', 'Beer'};
legendTerms = {'Historical Demand'};
for i = 1:dimOutput
    legendTerms{i+1}=[num2str(trainForesight(i)),'-Timestep Prediction'];
end

figure
for i = 1:nz
    skuInds = (1+nTest*(i-1)):(nTest*i);
    subplot(3,1,i)
    plot(xlTimes,realTS(i,:))
    xlabel('Date');
    ylabel(['Demand for ', prodNames{i}]);
    datetick('x','keepticks','keeplimits')
    hold on
    for j = 1:dimOutput
        scatter(xlTimes(predTimeSteps(skuInds)+trainForesight(j)), realPreds(j,skuInds),symbolKey(j))
    end
    hold off
    if (i==nz)
        legend(legendTerms);
        legend('Location','bestoutside');
    end
end
