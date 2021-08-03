function [ nLogL, logLcontr, varargout ] = ...
    matvLaplaceWishlike( Sigma_, n, X, varargin)
%MATVLAPLACEWISHLIKE
%
% USAGE:
%   
%
% INPUTS:
%   
%
% OUTPUTS:
%   
%  See also 
%
% COMMENTS:
%   
% REFERENCES:
%      [1]                     
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 05.05.2021
%
% DEPENDENCIES:
%

narginchk(3,4);
nargoutchk(0,5);

[p,~,N] = size(X);
p_ = p*(p+1)/2;
%% Parameters
if nargin >= 4
    all_param = varargin{1};
    for ii = 1:size(all_param,2)
        try % Get rid of this try/catch as soon as you got used to T always begin last index regime!
        Sigma_(:,:,ii) = ivechchol(all_param(1:p_,ii));
        n(ii) = all_param(p_ + 1, ii);
        catch
        error('Have you input params in (#param by T) format? (T should always be last index.)')
        end
    end
end

if size(n,1) ~= 1 || size(n,1) ~= N
    error('Input parameters must be either single objects or of same length as data.')
end
param.Sigma_ = Sigma_;
param.n = n;
for ii = 1:size(n,1) % This serves as a p.d. check on the Sigma_'s. Optionally this param struct can be returned.
    param.all(:,ii) = [vechchol(Sigma_(:,:,ii)); n(ii)];
end
if nargout >= 5    
    varargout{3} = param;
end
%% Log-Likelihood
logLcontr = NaN(N,1);
 
for ii = 1:N
    A = X(:,:,ii);
    Sig = Sigma_(:,:,ii);
    n_ = n(ii);
    trQ = trace(Sig\A);
    
    term1 = (2-p*n_)/2*log(2);
    term2 = -n_/2*logdet(Sig);
    term3 = -mvgammaln(n_/2,p);
    log_norm_const = term1 + term2 + term3;    
    
    term4 = (n_-p-1)/2*logdet(A);
    term5 = (2-n_*p)/4*log(trQ/2);
    term6 = log(besselk((2-p*n_)/2,sqrt(2*trQ)));
    log_kernel = term4 + term5 + term6;
    
    logLcontr(ii) = log_norm_const + log_kernel;
end
nLogL = -sum(logLcontr);
%% Gradient
if nargout >= 3
    
    gradient.Sigma_ = NaN(N,p_);
    gradient.cholSigma = NaN;
    gradient.n = NaN;
    gradient.nu = NaN;
    gradient.all = NaN;
    
    for ii = 1:N        
        A = X(:,:,ii);
        Sig = Sigma_(:,:,ii);
        n_ = n(ii);
        trQ = trace(Sig\A);
        iSigAiSig = Sig\A/Sig;
        
        % General matrix derivative (ignoring symmetry of Sigma_):
        S = -n_*inv(Sig) ...
            -( (2-p*n_)/2/sqrt(2*trQ) - besselk( 2-p*n_/2, sqrt(2*trQ)) ...
                                      /exp(term6) ) ...
             /exp(term2)/sqrt(2*trQ)*iSigAiSig ...
            +(2-p*n_)/trQ*iSigAiSig;
        S = .5*S;
        
        % Accounting for symmetry of Sigma_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
        
        gradient.Sigma_(ii,:) = vech(S);

        % The gradient below are also easy to get quering wolframalpha.com 
        % with eg "d/da (log(Gamma(1/2 (a+ p n))))".
        % I am just too lazy to write them down right now.
        gradient.n(ii) = NaN;
        gradient.nu(ii) = NaN;
    
    end
    
    varargout{1} = gradient;
    
end
%% Hessian
if nargout >= 4
    
    hessian.Sigma_ = NaN;
    hessian.cholSigma = NaN;
    hessian.n = NaN;
    hessian.nu = NaN;
    hessian.all = NaN;    
    
    varargout{2} = hessian;
end
end
