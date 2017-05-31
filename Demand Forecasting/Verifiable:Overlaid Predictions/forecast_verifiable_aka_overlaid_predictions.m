clear all;
close all;


%% Importing Data
%[xlTS,xlSkus,xlTimes] = import_spreadsheet; %import_spreadsheet currently returns one product's time-series

[xlTS,xlSkus,xlTimes] = import_spreadsheet_10_24_16;
timeInc = xlTimes(2)-xlTimes(1);
realSkus = xlSkus;
xlSkus = 1:length(xlSkus);
disp('Data imported');
%% Data Constants
%%Note that (delay_constant+1 <= maxTrain-nTest) must be TRUE (see trainSamples)

nz = size(xlTS,1); %Number of products for which demand is recorded
nloc = 10; %Number of locations at which demand is recorded
%nt = 10000; %Number of time-steps at which demand is recorded
nt = size(xlTS,2);


delay_line = [0,1,2,8];
nDelay = length(delay_line);
delay_constant = max(delay_line);

nInfo = nDelay + 3; %Time-step, number of delay terms, first non-zero time-step, and (for now) sku

nFactors = 20;%Convolutional memory terms
dimInput = nFactors + nInfo; %Number of convolutional training terms, plus nDelay+1 for time-step and delayed demand, plus sku as inputs
trainAhead = 55;
trainForesight = 1:trainAhead; %Training on demand(i+1), demand(i+10), demand(i+50), demand(i+100)
dimOutput = length(trainForesight);


maxTrain = nt - max(trainForesight); %The highest time-step at which training
%can happen. This is because for a training input, there needs to be an expected
%output consisting of future demand. If a late time-step is used, there IS
%no future demand. 

%%%%%NN Constants
nTest =3;
nTrain = floor((maxTrain-delay_constant-nTest)); %Number of training samples

%% Normalizing Input
productMeans = mean(xlTS,2);
xlTS = xlTS - repmat(productMeans,1,nt);
productStds = ((1/nt)*sum(xlTS.^2,2)).^((1/2));
xlTS = xlTS./repmat(productStds,1,nt);

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
    a=sine_freq_output(timeSeries0(i,:),maxTrain,candidates);
    if (length(a) > nnominees)
        sinfreqs(i,:) = takesome(a);
    else
        sinfreqs(i,:) = -1*ones(1,nnominees); %If no nominees are generated (peaks of convolution norm), placehold with negative frequencies    
    end
end
[~,sinParams] = kmeans(sinfreqs(find(sinfreqs>0)),nSin); %cluster non-zero frequencies, collect nSin of them.
clear sinfreqs; clear takesome; clear candidates; clear nnominees;



% Wavelet Parameters
nnominees = floor(nWav);
takesome = @(x) x(1:nnominees);
wavScales = zeros(nz,nnominees);
candidates = 1:.5:200; %Candidate period parameters
for i = 1:nz
    wavScales(i,:) = takesome(wavelet_param_output(timeSeries0(i,:),maxTrain,candidates));
end
wavScales = wavScales(:);
wavScales(wavScales==0)=[];%Padded with 10 zeros in 'wavelet_param_output.m'
[~,wavParams] = kmeans(wavScales,nWav);
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



%% Reshaping Data
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

trainRatio = .7;
valRatio = .3;


trainInd = sort(randperm(size(trainData,2),floor(size(trainData,2)*trainRatio)));
valInd =setdiff(1:size(trainData,2),trainInd);
testInd = size(trainData,2)+1:size(trainData,2)+size(testDataStrip,2);

trainData = [trainData,testDataStrip];
trainDesired = [trainDesired,testDesiredStrip];


%% Training the NN

net = feedforwardnet([15,35,15]);%one hidden layer of 20 nodes
net.trainParam.max_fail = 6; %Default is 6!
net.divideFcn = 'divideind';

net.divideParam.trainInd = trainInd;
net.divideParam.valInd =valInd;
net.divideParam.testInd = testInd;

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
symbolKey = 'o+*x^#w@$';
%Challenge: Include time-steps (ts-maxTrain-nTest):

predIndices = net.divideParam.testInd;
predTimeSteps = trainData(1,predIndices);
preds = netOutput(:,predIndices);
predSkus = trainData(nInfo,net.divideParam.testInd);


realTS = timeSeries0.*repmat(productStds,1,size(timeSeries0,2));
realTS = realTS + repmat(productMeans,1,size(realTS,2));

realPreds = preds.*repmat(productStds(predSkus)',dimOutput,1);
realPreds = realPreds + repmat(productMeans(predSkus)',dimOutput,1);

extraTimes = xlTimes(end) + ((trainForesight)*timeInc);
extraTS = realPreds(:,nTest*(1:nz)')'; %The forecast-series calculated at the final test sample of each product.

%%Saving Forecast Variables
%clearvars -except net extraTS extraTimes xlTimes realTS xlSkus nt nz nloc trainForesight
save('forecast_workspace_1_10.mat')