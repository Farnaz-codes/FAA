function p = modelOrder(s)
dirCODE = 'C:\Users\oaa30\Documents\MATLAB\Omar\ARfit';
addpath(genpath(dirCODE));
%% model order of the MVAR model ------------------
pMin = 1;
pMax = 30;             
sNorm = s;
%% MVAR order selection ------------------
[~, ~, ~, sbc]= arfit(sNorm, pMin, pMax, 'zero');
[~,p] = min(sbc);