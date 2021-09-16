function [ nLogL, logLcontr, varargout ] = matviWishlike( Sigma_, df, X, varargin )
%MATVIWISHLIKE
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
% Michael Stollenwerk
% michael.stollenwerk@live.com

narginchk(3,4);
nargoutchk(0,6);

[p,~,N] = size(X);
p_ = p*(p+1)/2;
%% Parameters
if nargin == 4
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1:p_));
    df = all_param(p_ + 1);
end
% Checking if Sigma_ is symmetric p.d.
param.Sigma_ = Sigma_;
param.df = df;
param.all = [vechchol(Sigma_); df];
%% Log-Likelihood
logLcontr = NaN(N,1);

log_norm_const = ...
  - df*p/2*log(2) ...
  - mvgammaln(df/2, p) ...
  + df/2*logdet(Sigma_);
	  
for ii = 1:N
    B = inv(X(:,:,ii));
 
    log_kernel = ...
      - trace(Sigma_*B) ...
      + (df+p+1)*logdet(B);
		
    logLcontr(ii) = log_norm_const + .5*log_kernel;
end
nLogL = -sum(logLcontr);
%% Scores
if nargout >= 3
    invSig = inv(Sigma_);    
    
    score.Sigma_ = NaN(N,p_);
    score.Sigma_WishFishScaling = NaN(p,p,N);
    score.df = NaN(N,1);
    for ii = 1:N
        
        invA = inv(X(:,:,ii));
       
        % General matrix derivative (ignoring symmetry of Sigma_):
        S = .5*( df*invSig - invA );
        
        score.Sigma_WishFishScaling(:,:,ii) = 2/df*Sigma_*S*Sigma_;

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
%% Hessian (Optional Output)

%% Optional Parameter Vector Output
if nargout >= 5
    varargout{3} = param;
end
%% Fisher Info (Optional Output)
if nargout >= 6
    
    G = Dmatrix(p);
    fisherinfo.Sigma_ = df/2*G'*kron(invSig,invSig)*G;
    
    fisherinfo.df = NaN;
    
    varargout{4} = fisherinfo;
    
end
end
