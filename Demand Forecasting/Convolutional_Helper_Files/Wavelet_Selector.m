close all


ts = timeSeries0(2,:);
endTime = size(timeSeries0,2);
depth = 200; %Also minimum of time-series domain range
inc = .5;
scales = 1:inc:200;
%ninth_step = floor(length(T)/9);
param1 = (1:6:54);
param2= zeros(1,9);


convolution_maxima = zeros(length(scales),1);
convolution_1Norms = zeros(length(scales),1);
convolution_2Norms = zeros(length(scales),1);
for i = 1:length(scales)
    convo = convolve_at(ts,wavelet(endTime-depth,scales(i),param2(1)),depth:endTime,depth);
    convolution_maxima(i) = max(abs(convo));
    convolution_1Norms(i) = sum(abs(convo));
    convolution_2Norms(i) = sum(convo.^2);
end
figure
plot(scales,convolution_maxima);
title('Maxima in Convolution')
xlabel('Period of Sine Wave')
ylabel('Absolute Maximum in Convolution')

figure
plot(scales,convolution_1Norms);
title('1-Norm of Convolution')
xlabel('Period of Sine Wave')
ylabel('1-Norm of Convolution')

figure
plot(scales,convolution_2Norms);
title('2-Norm of Convolution')
xlabel('Scale of Wavelet')
ylabel('2-Norm of Convolution')

[~,I] = sort(convolution_2Norms);
I = scales(I(end:-1:1));


figure
for i = 1:9
    subplot(3,3,i)
    %ind = ninth_step*i;
    plot(convolve_at(ts,wavelet(endTime-depth,I(i),param2(i)),depth:endTime,depth))
    title(['Convolution with wavelet of scale, offset: (', num2str(I(i)), ', ', num2str(param2(i)), ')']);
end
