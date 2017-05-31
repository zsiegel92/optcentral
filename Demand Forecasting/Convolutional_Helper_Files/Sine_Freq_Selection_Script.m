close all


ts = timeSeries0(i,:);
endTime = size(ts,2);
depth = floor(nt/2); %Also minimum of time-series domain range
inc = 0.5;
param_domain = 5:inc:200;
param_domain = candidates;

convolution_maxima = zeros(length(param_domain),1);
convolution_1Norms = zeros(length(param_domain),1);
convolution_2Norms = zeros(length(param_domain),1);
convolution_2Norms_unscaled = zeros(length(param_domain),1);

%convolution_maxima(i) = max(abs(convo));
%convolution_1Norms(i) = sum(abs(convo));
    
for i = 1:length(param_domain)
    convo = convolve_at_scaled(ts,sin(2*pi*(1/param_domain(i))*(1:endTime)),depth:endTime,depth);
    convolution_2Norms(i)=sum(convo.^2);
    
    convo_unscaled = convolve_at(ts,sin(2*pi*(1/param_domain(i))*(1:endTime)),depth:endTime,depth);
    convolution_2Norms_unscaled(i)=sum(convo_unscaled.^2);
end



%OPTIMIZING PARAMETERS

%FREQS AND PEAKS!!! :D

nplot = 4; %Choose a square integer

%freqs = peak_finder(convolution_2Norms,param_domain);
[freqs,pks] = peak_finder(convolution_2Norms,param_domain);
figure
for i = 1:nplot
    subplot(floor(sqrt(nplot)),floor(sqrt(nplot)),i)
    %ind = ninth_step*i;
    plot(convolve_at(ts,sin(2*pi*(1/freqs(i))*(1:453)),depth:endTime,depth))
    title(['Convolution with sin of period: ', num2str(freqs(i))]);
end



[freqs_unscaled_conv,pks_unscaled] = peak_finder(convolution_2Norms_unscaled,param_domain);

figure
for i = 1:nplot
    subplot(floor(sqrt(nplot)),floor(sqrt(nplot)),i)
    %ind = ninth_step*i;
    plot(convolve_at(ts,sin(2*pi*(1/freqs_unscaled_conv(i))*(1:453)),depth:endTime,depth))
    title(['UN-Scaled conv. w/ sin of period: ', num2str(freqs_unscaled_conv(i))]);
end





figure
subplot(2,1,1)
plot(param_domain,convolution_2Norms);
hold on
scatter(freqs(1:nplot),pks(1:nplot),'o');
hold off
title(['2-Norm of Scaled Convolution. Peaks: ', mat2str(freqs(1:nplot))])
xlabel('Period of Sine Wave')
ylabel('2-Norm of Convolution')

subplot(2,1,2)
plot(param_domain,convolution_2Norms_unscaled);
hold on
scatter(freqs_unscaled_conv(1:nplot),pks_unscaled(1:nplot),'o');
hold off
title(['2-Norm of UN-Scaled Convolution. Peaks: ', mat2str(freqs_unscaled_conv(1:nplot))])
xlabel('Period of Sine Wave')
ylabel('2-Norm of Convolution')


