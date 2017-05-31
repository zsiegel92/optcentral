function [ mu_over_d ] = gamma_param_output(ts,depth,param_domain1,param_domain2 ) %Domain 1 is scales, Domain 2 is centers
endTime = size(ts,2);
convolution_2Norms = zeros(length(param_domain1)*length(param_domain2),1);
indices = zeros(length(param_domain1)*length(param_domain2),2);
count = 1;
for i = 1:length(param_domain1)
    for j = 1:length(param_domain2)
        convo = convolve_at_scaled(ts,gammaTerm(endTime,param_domain2(j),param_domain1(i)),depth:endTime,depth);
        convolution_2Norms(count)=sum(convo.^2);
        indices(count,1) = i;
        indices(count,2)=j;
        count = count+1;
    end
end
[count_index,pks] = peak_finder(convolution_2Norms,1:(count-1));
mu_over_d = zeros(2,length(count_index));
disp(['size(indices): ', num2str(size(indices))]);
disp(['size(count_index): ', num2str(size(count_index))]);
for i = 1:size(mu_over_d,2)
    mu_over_d(:,i)=[indices(count_index(i),1);indices(count_index(i),2)];
end
disp(['Size of parameter suggestions: ', num2str(size(mu_over_d))]);
% 

% %FREQS AND PEAKS!!! :D
% 
% nplot = 9; %Choose a square integer
% figure
% subplot(2,1,1)
% plot(param_domain1,convolution_2Norms);
% hold on
% scatter(scales(1:nplot),pks(1:nplot),'o');
% hold off
% title(['2-Norm of (Scaled) Convolution. Peaks: ', mat2str(scales(1:nplot))])
% xlabel('Scale of Wavelet')
% ylabel('2-Norm of Convolution')
end

