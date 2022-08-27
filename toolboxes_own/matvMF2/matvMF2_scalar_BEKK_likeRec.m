function [ nLogL, logLcontr, fit, fcst, varargout ] = ...
    matvMF2_scalar_BEKK_likeRec( param, R, dist )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.08.2020

m = 63;
L = 1;
t_ahead = 1;
p = size(R,1);

[k,~,T] = size(R);
k_ = k*(k+1)/2;
%% Parameters
intrcptTau = ivechchol(param(1:k_));
alphaTau = param(k_+1);
betaTau = param(k_+2);
alphaS = param(k_+3);
betaS = param(k_+4);
omega = param(k_+5);
theta = param(k_+6:end);

% Initialize recursions with backcast of the realized covariance matrices. 
% Inspired by Sheppard MFE Toolbox.
lag = ceil(sqrt(T));
w = .06 * .94.^(0:(lag-1));
w = reshape(w/sum(w),[1 1 lag]);
backCast = sum(bsxfun(@times,w, R(:,:,1:lag)),3);
iniTau = backCast;

%% Data Storage
SigmaS = NaN(k,k,T+t_ahead);
SigmaTau = NaN(k,k,T+t_ahead);
SigmaOverall = NaN(k,k,T+t_ahead);
Z = NaN(k,k,T+t_ahead);
I = eye(k);
logLcontr = NaN(T,1);
%% Recursion
SigmaTau(:,:,1) = iniTau;
Ctau = chol(SigmaTau(:,:,1),'lower');
SigmaOverall(:,:,1) = iniTau;
Z(:,:,1) = R(:,:,1);
SigmaS(:,:,1) = I;

[~, logLcontr(1)] = matvsLogLike( dist, SigmaOverall(:,:,1), theta, R(:,:,1) );
for tt=2:T
    
    SigmaS(:,:,tt) = (1-alphaS-betaS)*I ...
                        + alphaS*(Ctau\R(:,:,tt-1)/Ctau') ...
                        + betaS*SigmaS(:,:,tt-1);
    
    if tt-L*m < 1
% The decay function is coded below. It takes as arguments the parameters 
% omega, "numer of lags" L, "averaging length" m, "data average" X and in
% case L*m is larger than length of data you can put the optional argument
% how you want to initialize X.

%         rollMeanZ = decay(omega,L,m,R(:,:,1:tt-1),iniTau);
        rollMeanZ = decay(omega,L,m,Z(:,:,max(tt-L*m,1):tt-1),iniTau);
    else
%         rollMeanZ = decay(omega,L,m,R(:,:,tt-L*m:tt-1));
        rollMeanZ = decay(omega,L,m,Z(:,:,tt-L*m:tt-1));
    end
    
%     rollMeanZ = sum(Z(:,:,max(tt-m,1):tt-1),3)/min(m,tt-1) + I*length(tt-m:0)/m;
    SigmaTau(:,:,tt) = intrcptTau ...
                        + alphaTau*rollMeanZ ...
                        + betaTau*SigmaTau(:,:,tt-1);
    
    try
        Ctau = chol(SigmaTau(:,:,tt),'lower');
        Cs = chol(SigmaS(:,:,tt),'lower');
        C = Ctau*Cs;
        SigmaOverall(:,:,tt) = C*C';
        Z(:,:,tt) = Cs\R(:,:,tt)/Cs';
        % Likelihood Evaluation         
        [~, logLcontr(tt), score, ~, param_dist] = ...
            matvsLogLike( dist, SigmaOverall(:,:,tt), theta, R(:,:,tt) );   
    catch ME
%         tt
%         ME.message
%         [ME.stack.line]'
%         {ME.stack.name}'
        nLogL = inf;
        logLcontr = NaN;
        fit = NaN;
        fcst = NaN;
        varargout{1} = NaN;
        return
    end
    
end
%% Fcst
for tt=T+1:T+t_ahead
    
    if tt == T+1
        
        SigmaS(:,:,tt) = (1-alphaS-betaS)*I ...
                            + alphaS*(Ctau\R(:,:,tt-1)/Ctau') ...
                            + betaS*SigmaS(:,:,tt-1);
        
%         rollMeanZ = decay(omega,L,m,R(:,:,tt-L*m:tt-1));
        rollMeanZ = decay(omega,L,m,Z(:,:,tt-L*m:tt-1));
%         rollMeanZ = mean(Z(:,:,tt-m:tt-1),3);
        SigmaTau(:,:,tt) = intrcptTau ...
                            + alphaTau*rollMeanZ ...
                            + betaTau*SigmaTau(:,:,tt-1);

        Ctau = chol(SigmaTau(:,:,tt),'lower');
        Cs = chol(SigmaS(:,:,tt),'lower');
        C = Ctau*Cs;
        SigmaOverall(:,:,tt) = C*C';
        Z(:,:,tt) = I;
        
    else
        
        SigmaS(:,:,tt) = (1-alphaS-betaS)*I ...
                            + alphaS*SigmaS(:,:,tt-1) ...
                            + betaS*SigmaS(:,:,tt-1);
        
%         rollMeanZ = decay(omega,L,m,R(:,:,tt-L*m:tt-1));
        rollMeanZ = decay(omega,L,m,Z(:,:,tt-L*m:tt-1));
%         rollMeanZ = mean(Z(:,:,tt-m:tt-1),3);
        SigmaTau(:,:,tt) = intrcptTau ...
                            + alphaTau*rollMeanZ ...
                            + betaTau*SigmaTau(:,:,tt-1);

        Ctau = chol(SigmaTau(:,:,tt),'lower');
        Cs = chol(SigmaS(:,:,tt),'lower');
        C = Ctau*Cs;
        SigmaOverall(:,:,tt) = C*C';
        Z(:,:,tt) = I;
        
    end
end
%% Log-Likelihood
nLogL = -sum(logLcontr);
%% Output
fit.SigmaS = SigmaS(:,:,1:T);
fit.SigmaTau = SigmaTau(:,:,1:T);
fit.Sigma_ = SigmaOverall(:,:,1:T);
fit.Z = Z(:,:,1:T);

fcst.SigmaS = SigmaS(:,:,T+1:end);
fcst.SigmaTau = SigmaTau(:,:,T+1:end);
fcst.Sigma_ = SigmaOverall(:,:,T+1:end);
fcst.Z = Z(:,:,T+1:end);

if nargout >= 4
    param_out.intrcptTau = intrcptTau;
    param_out.alphaTau = alphaTau;
    param_out.betaTau = betaTau;
    param_out.alphaS = alphaS;
    param_out.betaS = betaS;
    param_out.omega = omega;
    
    if isfield(param_dist,'n')
        param_out.n = param_dist.n;
    end
    
    if isfield(param_dist,'nu')
        param_out.nu = param_dist.nu;
    end

    param_out.all = param;
    
    varargout{1} = param_out;
end

end

function out_ = decay(omega,L,m,X,varargin)
% The decay function takes as arguments the parameters omega, 
%"numer of lags" L, "averaging length" m, "data to average" X and in
% case L*m is larger than length of data you can put the optional argument
% how you want to initialize X.

length_ = size(X,3);
if nargin == 5
    X = cat(3,repmat(varargin{1},1,1,L*m-length_),X);
end

if size(X,3) ~= L*m
    error('Something went wrong.')
end
X = flip(X,3);

rho = NaN(L,1);
if L == 1
    rho = 1;
else
    for l = 1:L
        %weights
        rho(l) = (1-l/L)^(omega-1);
    end
    rho = rho./sum(rho);
end
%rvm
out_ = 0;
for l = 1:L
    out_ = out_ + rho(l)*mean(X(:,:,m*(l-1)+1:m*l),3);
end

end
