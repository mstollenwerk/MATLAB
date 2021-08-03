function [ nLogL, logLcontr, varargout ] = matvWishlike( Sigma_, df, X, varargin )
%MATVWISHLIKE
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
% DEPENDENCIES:
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 14.08.2020

narginchk(3,4);
nargoutchk(0,5);

[p,~,N] = size(X);
p_ = p*(p+1)/2;
%% Parameters
if nargin == 4
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1:p_));
    df = all_param(p_ + 1);
end

if nargout == 5
    param.Sigma_ = Sigma_;
    param.df = df;
    param.all = [vechchol(Sigma_)' df];
    
    varargout{3} = param;
end
%% Input Checking
chol(Sigma_,'lower'); % Checking if Sigma_ is symmetric p.d.
%% Log-Likelihood
logLcontr = NaN(N,1);

log_norm_const = ...
  - df*p/2*log(2) ...
  - mvgammaln(df/2, p) ...
  - df/2*logdet(Sigma_);
	  
for ii = 1:N
    A = X(:,:,ii);
 
    log_kernel = ...
      -  trace(Sigma_\A) ...
      + (df-p-1)*logdet(A);
		
    logLcontr(ii) = log_norm_const + .5*log_kernel;
end
nLogL = -sum(logLcontr);
%% Scores
if nargout >= 3
    invSig = inv(Sigma_);    
    
    score.Sigma_ = NaN(N,p_);
    score.df = NaN(N,1);
    
    for ii = 1:N
        
        A = X(:,:,ii);
        
        % General matrix derivative (ignoring symmetry of Sigma_): [ enter to matrixcalculus.org: -df/2*log(det(Sigma))-tr(inv(Sigma)*A) ]
        S = - df/2*invSig ...
            + 1/2*(Sigma_\A/Sigma_);
        
        % Accounting for symmetry of Sigma_:
        S = 2*S - diag(diag(S));
        
        % Numerical inaccuracies make S not exactly symmetric
        S = (S+S')./2;
        
        score.Sigma_(ii,:) = vech(S);

        % The score below are also easy to get quering wolframalpha.com 
        % with eg "d/da (log(Gamma(1/2 (a+ p n))))".
        % I am just too lazy to write them down right now.
        score.df(ii) = NaN;
    
    end
    
    varargout{1} = score;
    
end
%% Fisher Info
if nargout >= 4
    
    G = Dmatrix(p);
    
    fisherinfo.Sigma_ = df/2*G'*kron(invSig,invSig)*G;
    
    varargout{2} = fisherinfo;
end
end
