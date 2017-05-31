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

Choices = {'Particle Swarm Optimization', 'Firefly Algorithm'};

ANSWER = questdlg('Select the algorithm:', ...
                  'Inventory Control', ...
                  Choices{1}, Choices{2}, ...
                  Choices{1});

if strcmpi(ANSWER, Choices{1})
    pso;
    return;
end

if strcmpi(ANSWER, Choices{2})
    fa;
    return;
end

