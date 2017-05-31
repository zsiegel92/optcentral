function [ scales ] = wavelet_param_output(ts,depth,param_domain )
endTime = size(ts,2);
convolution_2Norms = zeros(length(param_domain),1);
for i = 1:length(param_domain)
    convo = convolve_at_scaled(ts,wavelet(endTime,param_domain(i),0),depth:endTime,depth);
    convolution_2Norms(i)=sum(convo.^2);
end
[scales,pks] = peak_finder(convolution_2Norms,param_domain);

scales = [scales, zeros(1,15)];

% %FREQS AND PEAKS!!! :D
% 
% nplot = 9; %Choose a square integer
% figure
% subplot(2,1,1)
% plot(param_domain,convolution_2Norms);
% hold on
% scatter(scales(1:nplot),pks(1:nplot),'o');
% hold off
% title(['2-Norm of (Scaled) Convolution. Peaks: ', mat2str(scales(1:nplot))])
% xlabel('Scale of Wavelet')
% ylabel('2-Norm of Convolution')
end

