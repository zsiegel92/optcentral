%dim(data)== [1,nt]
%dim(basis)== [nt,1]
%dim(x) == number of samples of convolution, should be sorted
%dim(maxDepth) == 1
%dim(convolution) == dim(x) == number of samples of convolution
function convolution = convolve_at(data, basis, x, maxDepth)
    x = sort(x); %input could simply be sorted
    data = reshape(data,1,length(data)); %Ensure input is row vector
    basis = reshape(basis,length(basis),1); %Ensure input is column vector
    convolution = zeros(length(x),1);
    depth = min(x(1)-1,maxDepth-1);
    convolution(1) = data(x(1):-1:(x(1)-depth))*basis(1:(depth+1));
    for i = 2:length(x) %separate from convolution(1) because x(i)>=2 forall i >=2!
        depth = min(x(i)-1,maxDepth-1);
        convolution(i) = data(x(i):-1:(x(i)-depth))*basis(1:(depth+1));
    end
end