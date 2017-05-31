%NOTE: This function outputs some awesome plots, just uncomment to see!

function [ freqs ] = sine_freq_output(ts, depth,param_domain)
endTime = size(ts,2);
convolution_2Norms = zeros(length(param_domain),1);
for i = 1:length(param_domain)
    convo = convolve_at_scaled(ts,sin(2*pi*(1/param_domain(i))*(1:endTime)),depth:endTime,depth);
    convolution_2Norms(i)=sum(convo.^2);
end
[freqs,pks] = peak_finder(convolution_2Norms,param_domain);

% 
% %FREQS AND PEAKS!!! :D
% 
% nplot = 9; %Choose a square integer
% figure
% subplot(2,1,1)
% plot(param_domain,convolution_2Norms);
% hold on
% scatter(freqs(1:nplot),pks(1:nplot),'o');
% hold off
% title(['2-Norm of (Scaled) Convolution. Peaks: ', mat2str(freqs(1:nplot))])
% xlabel('Period of Sine Wave')
% ylabel('2-Norm of Convolution')
end

