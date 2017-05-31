%% Training the NN

net = feedforwardnet([30,30,30,30]);%one hidden layer of 20 nodes
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
