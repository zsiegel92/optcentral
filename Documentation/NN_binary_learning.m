%All numbers that can be represented with 9 bits.
des_output = (0:2^(15)-1);
sum = 0;
trainData= arrayfun(@(x) str2num(x),dec2bin(des_output)');

net = feedforwardnet([10,10]);
net = train(net,des_output,trainData);
netOutput = sim(net,des_output);
%net = train(net,trainData,des_output);
%netOutput = sim(net,trainData);

%It appears to be easier for a neural network to learn to convert binary
%arrays to integers than to convert integers to binary arrays.