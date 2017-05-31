function [timeSeries,skus,times] = import_spreadsheet_10_24_16( )
filename = 'Data_10-24-16.xlsx';
timeSeries = xlsread(filename);
skus = timeSeries(2:end,1); %Capture skus
timeSeries(:,1)=[]; %Delete product Numbers

skus(find(all(timeSeries==0,2)),:)=[];
timeSeries(find(all(timeSeries==0,2)),:)=[];

%array(1,1)=0; Delete Top-left corner nan
times = timeSeries(1,:);
times = times - min(times);
timeSeries(1,:) = [];

timeSeries = imputeData(timeSeries);%Fill in missing values

end

