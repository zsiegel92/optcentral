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

function [z, sol]=MyCost(xhat,model)

    sol=ParseSolution(xhat,model);
    
    SumAX=sol.SumAX;
    SumBI=sol.SumBI;
    VMIN=sol.VMIN;
    VMAX=sol.VMAX;
    
    alpha=100000;
    beta=10;
    z=((SumAX+SumBI)+alpha*VMIN)*(1+beta*VMAX);

end