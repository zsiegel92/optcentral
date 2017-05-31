%Try: a = gammaTerm(1000,150,.5);
%A function that peaks at d, an
function [gammaFn,x0] = gammaTerm(T,d,mu)
x0=1:T;
% scaleIssues = [];
% scaleFacs = [];
% scalecap = 100000;
%k must be at least 1
gammaFn = zeros(T,1);
%gammaFn(1:d-1) = 0;
%gammaFn(d) = ((1-mu)^(d+1));
gammaFn(d)=1;
for i = (d+1):T
%     if (gammaFn(i-1)>scalecap)
%         gammaFn(i-1) = gammaFn(i-1)/scalecap;
%         scaleIssues = [scaleIssues, i-1];
%     end
    gammaFn(i) = gammaFn(i-1)*(i/(i-d))*mu;
end
%gammaFn = gammaFn*((1-mu)^(d+1));
gammaFn = normc(gammaFn);
end

