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

function PlotSolution(sol,model)

    K=model.K;
    H=model.H;
    I0=model.I0;
    X0=zeros(K,1);
    u=model.u;
    Umax=model.Umax;
    UC0=sum(u.*I0);
    
    X=sol.X;
    I=sol.I;
    UC=sol.UC;

    subplot(3,1,1);
    stairs(0:H,[X X0]','LineWidth',2);
    xlabel('Time');
    ylabel('Order Amount');
    
    subplot(3,1,2);
    stairs(0:H,[I0 I]','LineWidth',2);
    xlabel('Time');
    ylabel('Inventory');

    subplot(3,1,3);
    stairs(0:H,[UC0 UC],'LineWidth',2);
    hold on;
    plot([0 H],[Umax Umax],'r:','LineWidth',2);
    hold off;
    xlabel('Time');
    ylabel('Used Capacity');
    
end