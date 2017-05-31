function [wave,x0] = wavelet(length,a,b)
%x0= (floor(-length/2):(floor(length/2)-1)); %Centering
x0 = 0:(length-1);
x = ((1/a)*(x0-b));
wave = (1/sqrt(a))*(2*sinc(2*pi*x)-sinc(pi*x));
%Sans x term
%wave = (1/sqrt(a))*(2*sinc(2*pi*((1/a)*(x0-b)))-sinc(pi*((1/a)*(x0-b))));
end

