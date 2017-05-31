function [ausData, skus, times] = Aussie_data_test ()
ausData1 = xlsread('Beer_quarterly.xlsx');
ausData2_3 = xlsread('Steel_and_Gas.xlsx');
ausData1 = [ausData1,ausData1,ausData1];
ausData1 = ausData1';
ausData1= ausData1(:);
ausData = [ausData2_3,ausData1];
clear ausData1; clear ausData2_3;

times = ausData(:,1)';
ausData(:,1)=[];
ausData = ausData';
skus = (1:3)';
times = times - min(times);

end