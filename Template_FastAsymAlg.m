% Fast Asymptotic Algorithm 
% Template
%=============================================================================#
% Times series for analysis Generated with ESSimEEGSource 
% [37] E. Barzegaran et al., "EEGSourceSim: A framework for realistic simulation
% of EEG scalp data using MRI-based forward models and biologically plausible
% signals and noise," Journal of neuroscience method, vol.328, p. 108377, 2019.
%=============================================================================#
% Steps
% 1. Load signal S with dimension ns x K
% 2. Normalized signal S
% s=(S-mean(S))./std(S);
% 3. Find model order (using any knid of criterion) , e.g. P = 3
% 4. Define alpha, e.g. alpha = 0.01
% 5. Define Frequency band, e.g. Freq = 1:8;  %Freq: frequency band

%=============================================================================#
Current_dir = pwd;
addpath(genpath([Current_dir,'\Fast Asymptotic Codes']));
%=============================================================================#
% Estimation algorithm [22]:
% Modified Partial Correlation Estimation: Vieira-Morf [2,5] using unbiased covariance estimates.
[ARcoef,~,Ecov] = mvar(s2, P, 22); 
% Ecov: Covariance of noise
% ARCoef: MVAR Coeff
% alpha: significance level of statistical test
%Freq: frequency band
% FastAsympAlg: s2: time series, ,
% SF: sampling frequency
measure = 'dtf';
df = FastAsympAlg(s2,ARcoef,Ecov,Freq,measure,SF,alpha);
PHI = df.Phi;
CIupperbound = df.CIupperbound;
CIlowerbound = df.CIlowerbound;
Pval = df.Pval;
Threshold = df.Threshold;
Phi_th = df.Phi_th;
%=============================================================================#
