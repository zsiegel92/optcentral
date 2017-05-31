function [expTerm,x0] = expTerm(length,mu)
x0 = 0:(length-1);
expTerm = (1-mu)*mu.^x0;
end

