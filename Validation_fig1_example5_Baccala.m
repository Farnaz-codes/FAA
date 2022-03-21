%
% Example taken from:
% Baccala & Sameshima. Partial directed coherence: a new concept in neural 
% structure determination. 
% _Biol. Cybern._ *84*:463--474, 2001.
%
% <http://dx.doi.org/10.1007/PL00007990>
% 
% Example Five-dimensional VAR[2] with loop and feedback

%%
%% Data sample generation
% 
% nDiscard = 10000;    % number of points discarded at beginning of simulation
% nPoints  = 100;    % number of analyzed samples points


% u = fbaccala2001a_ex5( nPoints, nDiscard );
%-----------------------------------------------------------------%
Current_dir = pwd;
addpath(genpath([Current_dir,'\Fast Asymptotic Codes']));
addpath(genpath([Current_dir,'\Baccala Codes\routines']));
% The data are generated applying "fbaccala2001a_ex5.m" in asymp_package_v3 [24]
% The data are saved in 'u'
load u
disp('Running MVAR estimation routine...')
p=2;
% builtin least squares methods
[A,pf,~,~]=idMVAR(u1,p,0);
A2=reshape(A,5,5,p);
nFreqs = 32;
fs=64;
alpha=0.01;
metric='info';
cf=FastAsympAlg(u1,A2,pf,1:nFreqs,'ipdc',[],alpha);
df=FastAsympAlg(u1,A2,pf,1:nFreqs,'idtf',[],alpha);

cb = asymp_pdc(u1,A2,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics
db = asymp_dtf(u1,A2,pf,nFreqs,metric,alpha); % Estimate DTF and asymptotic statistics
validationplot(cb,cf,'pdc',fs)
validationplot(db,df,'dtf',fs)
%-------------------------------------------------------------------------%