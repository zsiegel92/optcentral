function [ mus ] = exp_param_output(ts, depth,param_domain)
endTime = size(ts,2);
convolution_2Norms = zeros(length(param_domain),1);
for i = 1:length(param_domain)
    convo = convolve_at_scaled(ts,expTerm(endTime,param_domain(i)),depth:endTime,depth);
    convolution_2Norms(i)=sum(convo.^2);
end
[mus,pks] = peak_finder(convolution_2Norms,param_domain);
if (length(mus) < 9)
    [pks,mus] = sort(convolution_2Norms,'descend');
    mus = param_domain(mus);
    [mus2,pks2] = peak_finder(convolution_2Norms,param_domain);
    mus = [mus2,mus(1:8)];
end
% plot(param_domain,convolution_2Norms)
% title('2-Norms of Convolution w/ Exp(mu) vs. mu')
end

