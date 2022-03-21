function c = FastAsympAlg(s,ARcoef,Ecov,Freq,measure,SR,alpha)
% inputs :
% s      : signal - [time x nchannel]
% ARcoef : Matrix of coefficients estimated by MVAR
% Ecov   : Covariance of noise
% SR	 : Sample Rate
% alpha  : significance level of statistical test
% Freq   :  frequency band
% measure:
%          'DTF, iDTF, gDTF ( with metric='gamma_group' )
%          'PDC, iPDC, gPDC (with metric 'pi_group')
%-------------------------------------------------------------------------%
% Outputs:
% c
%c.CIupperbound : Upper bound (confidence interval)
%c.CIlowerbound : Lower bound (confidence interval)
%c.Pval         : P value
% c.Threshold   : Threshold
%c.Phi          : Estimated connectivity measure
%c.Comment      : Comment
%c.Phi_th :     : (alpha-1)% statistically significant Connectivity measures
% c.SS          : Spectral density (SS)
%-------------------------------------------------------------------------%
% Other References :
%[1]Baccalá, L.A., De Brito, C.S., Takahashi, D.Y. and Sameshima, K., 2013.
%Unified asymptotic theory for all partial directed coherence forms. Phil. Trans. R. Soc. A, 371(1997), p.20120158.
%[2] Baccalá, L.A., Takahashi, D.Y. and Sameshima, K., 2016. Directed transfer function: Unified asymptotic theory and some of its implications. IEEE Transactions on Biomedical Engineering, 63(12), pp.2450-2460.
%-------------------------------------------------------------------------%
% Author : Farnaz Rezaei 10/4/2019
%          last version  7/3/2020              
%-------------------------------------------------------------------------%
if (nargin<7 || alpha==0); alpha=0;c.Comment=['Estimation ', upper(measure)];else
    c.Comment=['ASymptotic distribution of ', upper(measure), ' and its ', num2str(100*(1-alpha)),' % Confidence interval'];end

if isempty(SR); SR=2*Freq(end);end
disp(c.Comment)
measure=lower(measure);
c.measure=measure;
metric='pi_group';
if isempty(regexp(measure,'pdc', 'once'));metric='gamma_group';end
if (isempty(regexp(measure,'pdc', 'once')) && isempty(regexp(measure,'dtf', 'once')) || length(measure)>4)
    error('Unknown metric.')
end
if (length(measure)==4 && (strcmp(measure(1),'i')==0 && strcmp(measure(1),'g')==0))
    error('Unknown metric.')
end
[np, nChannels] = size(s);
if np<nChannels;s=s';[np, nChannels] = size(s);end
s = s-mean(s);
icdf_norm_alpha = norminv(1-alpha/2.0,0,1);
if ~issymmetric(Ecov)
   Ecov = triu(Ecov,1)+ triu(Ecov)';
end
%-------------------------------------------------------------------------%
p = size(ARcoef,3);                           % model order
AR=reshape(ARcoef,nChannels,nChannels*p);
Ecovdotproduct=Ecov.*Ecov';                  % Hadamard product of Covariance of noise
%-------------------------------------------------------------------------%
%% Sn and Sd table(4)
switch metric
    case 'gamma_group'
        evar_d = eye(nChannels).*Ecov;
        if strcmp(measure,'dtf');evar_d=eye(nChannels);end
        Sd=evar_d;
        if strcmp(measure,'idtf');Sd=Ecov;end
        
    case 'pi_group'
        evar_d = inv((eye(nChannels).*Ecov));
        if strcmp(measure,'pdc');evar_d=eye(nChannels);end
        Sd=evar_d;
        if strcmp(measure,'ipdc');Sd=pinv(Ecov);end
        ce = diag(evar_d).*diag(Ecov);
end
%-------------------------------------------------------------------------%
%% Initialization
Phi = zeros(nChannels, nChannels, length(Freq));
pvalues=zeros(nChannels, nChannels, length(Freq));
TH = zeros(nChannels, nChannels, length(Freq));
CIup=zeros(nChannels, nChannels, length(Freq));
CIlow=zeros(nChannels, nChannels, length(Freq));
SS=zeros(nChannels, nChannels, length(Freq));
gamma = zeros(nChannels*p, nChannels*p);
%%------------------------------------------------------------------------%
%% Estimation of gamma and its decomposition
for m = 1:p
    for n = m:p
        SlagOut = s(n-m+1:(end-m+1),:)'*s(1:(end-n+1),:)/np;
        gamma(((m-1)*nChannels+1):m*nChannels, ((n-1)*nChannels+1):n*nChannels)=SlagOut;
        gamma(((n-1)*nChannels+1):n*nChannels, ((m-1)*nChannels+1):m*nChannels)=SlagOut';
    end
end
invgamma=pinv(gamma);
%%------------------------------------------------------------------------%
%% Frequency loop
for ff = 1:length(Freq)
    f1 = Freq(ff);
    f = (f1-1)/(SR);  %f starting at 0
    C1 = cos(2*pi*f*(1:p));
    S1 = sin(2*pi*f*(1:p));
    C = [C1;-S1];
    Af = eye(nChannels)-AR*kron(C1-1i*S1,eye(nChannels)).';
    %%------------------------------------------------------------------------%
    Hf = pinv(Af);
    SS(:,:,ff)=Hf*Ecov*Hf';
    OmegaX=zeros(nChannels,nChannels);
    OmegaEcov=zeros(nChannels,nChannels);
    %% Estimation of measures
    % phi=x'*Tn*x/x'*Td*x;
    switch metric
        case 'gamma_group'
            X=Hf;
            num=abs(X).^2*evar_d;
            den=sum(num,2);
            if strcmp(measure,'idtf')
                den=diag(real(X*Ecov*X'));
            end
        case 'pi_group'
            X=Af;
            num=evar_d*abs(X).^2;
            den=sum(num);
            if strcmp(measure,'ipdc')
                den=diag(real(X'*pinv(Ecov)*X));
                den=den';
            end
    end
    Phi(:,:,ff) = num./den;
    %-------------------------------------------------------------------------%
    if alpha~=0
        % If the Null hypothesis is rejected, Threshold is: (I add this if later)
        switch metric
            case 'gamma_group'
                %------------------------------------------------------------------------%
                % Estimating df (degree of freedom and multiplier)
                % without eigenvalues estimation and cholesky factors
                % Defining gamma-group parameters 
                F1=0.5*[1,1j;-1j,1];
                HEcov =X*Ecov;
                miu1 = sum(HEcov.*X,2);
                miu2 = sum(HEcov.*conj(X),2);
                temp = kron(eye(p),X.')*invgamma;
                Gp_gamma_1= temp*kron(eye(p),X);
                Gp_gamma_2= temp*kron(eye(p),conj(X));
                %---------------------------------------------------------%
                T_th_1 = nan(1,nChannels);T_th_2=T_th_1;T_th_3=T_th_1;
                P2 = nan(2,2,nChannels);P1=P2;
                 for j=1:nChannels
                     tempj_1 = Gp_gamma_1(j:nChannels:nChannels*p,j:nChannels:nChannels*p);
                     tempj_2 = Gp_gamma_2(j:nChannels:nChannels*p,j:nChannels:nChannels*p);
                     P2(:,:,j) = F1*C*tempj_2*C'*F1;
                     P1(:,:,j) = F1*C*tempj_1*C'*F1.';
                     
                     T_th_1(j)=real(trace(P2(:,:,j))); %trace (P2)
                     T_th_2(j) = real(sum(P2(:,:,j).*P2(:,:,j).','all')); % trace (P2*P2)
                     T_th_3(j) = real(sum(P1(:,:,j).*P1(:,:,j)','all'));  % trace (P2*P2)
                 end
                 mult = abs(miu1).^2*T_th_3+abs(miu2).^2*T_th_2;
                 patdf = (2*(abs(miu2)*T_th_1).^2)./mult;
                 patmul=(mult*evar_d)./((den.*abs(miu2))*T_th_1);
         %-----------------------------------------------------------------%
            case 'pi_group'
                cgamma =nan(1,nChannels); nugamma = nan(1,nChannels);
                for j=1:nChannels
                    temp=C*invgamma(j:nChannels:nChannels*p,j:nChannels:nChannels*p)*C';
                    cgamma(j)=sum(temp(:).^2)/(den(j)*trace(temp));
                    nugamma(j) =trace(temp).^2/sum(temp(:).^2);
                end
                patmul = ce*cgamma;
                patdf = repmat(nugamma,nChannels,1);
        end    
        TH(:, :, ff) = real(patmul).*icdf('chi2',(1-alpha), real(patdf))./(np);
        pvalues(:,:,ff) = 1 - chi2cdf(Phi(:,:,ff)*np./real(patmul), real(patdf));
        %-------------------------------------------------------------------------%
        %% Confidence interval
        %% Covariance of measures with respect to the X :OmegaX 
        switch metric
            case 'gamma_group'
                % --------------------------------------------------------%
                 %Eidted -8/25/2020
                 Xbar = [real(X),imag(X)];
                 temp=Xbar*kron(F1*C,Sd);
                 gsigma1 = sum(temp*Gp_gamma_1.*temp,2);
                 gsigma2 = sum(temp*Gp_gamma_2.*conj(temp),2);
                 T2 = real(miu1.*gsigma1+miu2.*gsigma2);
                 for j=1:nChannels
                      temp1 = sum(Xbar(:,[j,j+nChannels])*P1(:,:,j).*Xbar(:,[j,j+nChannels]),2);
                      temp2 = sum(Xbar(:,[j,j+nChannels])*P2(:,:,j).*Xbar(:,[j,j+nChannels]),2);
                      T1 = real(miu1.*temp1+miu2.*temp2);
                      temp3 = Xbar(:,[j,j+nChannels])*F1*C;
                      temp3_1=sum((temp3*Gp_gamma_1(j:nChannels:nChannels*p,:)).*temp,2);
                      temp3_2=sum((temp3*Gp_gamma_2(j:nChannels:nChannels*p,:)).*conj(temp),2);
                      T3 = 2*real(miu1.*temp3_1+miu2.*temp3_2);                     
                      OmegaX(:,j) = 8*(evar_d(j,j)^2*T1 + Phi(:,j,ff).^2.*T2-evar_d(j,j)*Phi(:,j,ff).*T3)./(den.^2);
                 end
                 %--------------------------------------------------------%
                % updated 6/28/2020
            case 'pi_group'
                for j=1:nChannels
                    Gp =4*C*invgamma(j:nChannels:nChannels*p,j:nChannels:nChannels*p)*C'/(den(j)^2);
                    Xbar = [real(X(:,j)),imag(X(:,j))];
                    temp1 = Xbar*Gp*Xbar';
                    switch measure
                        case 'pdc'
                            temp2 = Ecov.*temp1;
                            T1 = diag(Ecov).*diag(temp1);
                            T2 = sum(temp2(:));
                            T3 = sum(temp2,2);
                            
                        case 'gpdc'
                            temp2 = Ecov.*temp1;
                            T1 = diag(evar_d).*diag(temp1);
                            T3 = diag(evar_d).*temp2*diag(evar_d);                       
                            T2 = sum(T3);
                            
                        case 'ipdc'
                            temp2 = Sd.*temp1;
                            T1 = diag(evar_d).*diag(temp1);
                            T2 = sum(temp2(:));
                            T3 = T1;
                    end
                    OmegaX(:,j)=T1+Phi(:,j,ff).^2.*T2-2*Phi(:,j,ff).*T3;    
                end
                %---------------------------------------------------------%
        end
        % end of Confidenc interval-part1
        %-------------------------------------------------------------------------%
        %% 2) Covariance of measures with respect to the noise :OmegaEcov  updated 6/19/2020
        if (~strcmp(measure,'pdc') || ~strcmp(measure,'dtf') )
            switch metric
                case 'pi_group'
                    T1=(evar_d^2)*(abs(X).^4);
                    switch measure
                        case 'gpdc'
                            AbsB2 = abs(X).^2;
                            temp2 = AbsB2'*(evar_d^2)*Ecovdotproduct*(evar_d^2)*AbsB2;
                            T2=repmat(diag(temp2)',nChannels,1).*Phi(:,:,ff).^2;
                            T3 = AbsB2.*Phi(:,:,ff).*((evar_d^2)*Ecovdotproduct*(evar_d^2)*AbsB2);
                            OmegaEcov=2*(T1+T2-2*T3)./den.^2;
                       %--------------------------------------------------%
                        case 'ipdc'
                            T3= T1.*Phi(:,:,ff);   
                            temp2=0.5*(abs(X'*Sd*X).^2+abs(X.'*Sd*X).^2);
                            T2= repmat(diag(temp2)',nChannels,1).*Phi(:,:,ff).^2;
                            OmegaEcov=2*(T1-2*T3+T2)./den.^2;
                    end
          %----------------------------------------------------------------%
                case 'gamma_group'
                    T1=abs(X).^4*(evar_d^2);
                    AbsH2 = abs(X).^2;
                    switch measure
                        case 'gdtf'
                            temp2 = diag(AbsH2*Ecovdotproduct*AbsH2');   
                            T2=repmat(temp2,1,nChannels).*Phi(:,:,ff).^2;
                            T3 = AbsH2.*Phi(:,:,ff).*(AbsH2*Ecovdotproduct);
                            OmegaEcov=2*(T1+T2-2*T3)./repmat(den,1,nChannels).^2;
                       %--------------------------------------------------%
                        case 'idtf'
                            T3=Phi(:,:,ff).*AbsH2.*(HEcov.*(conj(X)*Ecov)); 
                            temp2 = 0.5*(abs(miu1).^2+abs(miu2).^2);
                            T2=repmat(temp2,1,nChannels).*(Phi(:,:,ff).^2);
                            OmegaEcov=2*(T1+T2-2*T3)./repmat(den,1,nChannels).^2;
                    end
         %----------------------------------------------------------------%
            end
        end   % end of Confidenc interval-part2
        %-------------------------------------------------------------------------%
        OmegaT=(OmegaEcov+OmegaX)/np;
        CIup(:,:,ff)=Phi(:,:,ff)+sqrt(OmegaT)*icdf_norm_alpha;
        CIlow(:,:,ff)=Phi(:,:,ff)-sqrt(OmegaT)*icdf_norm_alpha;
    end
end
Phi_th=((abs(Phi)-abs(TH)) > 0).*Phi;
Phi_th(Phi_th==0)=nan;
%-------------------------------------------------------------------------%
if alpha~=0
    c.CIupperbound=CIup;
    c.CIlowerbound=CIlow;
    c.Pval=pvalues;
    c.Phi=Phi;
    c.Comment=['ASymptotic distribution of ', upper(measure), ' and its ', num2str(100*(1-alpha)),' % Confidence interval'];
    c.Threshold=TH;
    c.Phi_th=Phi_th;
    c.SS=SS;
else
    c.Phi=Phi;
    c.SS=SS;
end
%-------------------------------------------------------------------------%
%=========================================================================%
% The End
%=========================================================================%
