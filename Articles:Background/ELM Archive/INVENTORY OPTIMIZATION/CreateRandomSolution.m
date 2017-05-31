%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPAP111
% Project Title: Inventory Control using PSO and Firefly Algorithm
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function xhat=CreateRandomSolution(model)

    K=model.K;
    H=model.H;
    
    xhat=rand(K,H);
    
end