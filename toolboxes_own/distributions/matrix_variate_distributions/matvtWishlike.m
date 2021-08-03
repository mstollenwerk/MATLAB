function [ nLogL, logLcontr, varargout ] = matvtWishlike( Sigma_, n, nu, X, varargin )
%MATVTWISHLIKE
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
%   This is the distribution of the quadratic form XX', where the columns
%   of the p by df_n matrix X follow multivariate central t distributions
%   with the same degree of freedom df_t and the same dispersion matrix 
%   Sigma_.
%   
% REFERENCES:
%      [1]                    
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.08.2020
%
% DEPENDENCIES:

narginchk(4,6);
nargoutchk(0,6);

[p,~,N] = size(X);
p_ = p*(p+1)/2;
%% Parameters
if nargin >= 5
    all_param = varargin{1};
    for ii = 1:size(all_param,2)
        try % Get rid of this try/catch as soon as you got used to T always begin last index regime!
        Sigma_(:,:,ii) = ivechchol(all_param(1:p_,ii));
        n(ii) = all_param(p_ + 1, ii);
        nu(ii) = all_param(p_ + 2, ii);
        catch
        error('Have you input params in (#param by T) format? (T should always be last index.)')
        end 
    end
end
if size(n,1) ~= 1 || size(n,1) ~= N
    error('Input parameters must be either single objects or of same length as data.')
end

% Center distribution to have mean Sigma_?
if nargin >= 6
    centering = varargin{2};
end
if centering
    for ii = 1:size(n,1)
        c(ii) = n(ii)*nu(ii)/(nu(ii)-2);
        Sigma_(:,:,ii) = Sigma_(:,:,ii)./c(ii);
    end
end

param.Sigma_ = Sigma_;
param.n = n;
param.nu = nu;
for ii = 1:size(n,1) % This serves as a p.d. check on the Sigma_'s. Optionally this param struct can be returned.
    param.all(:,ii) = [vechchol(Sigma_(:,:,ii)); n(ii); nu(ii)];
end
%% Log-Likelihood
logLcontr = NaN(N,1);
trQ = NaN(N,1);

for ii = 1:N
    A = X(:,:,ii);
    Sig = Sigma_(:,:,ii);    
    n_ = n(ii);
    nu_ = nu(ii);
    trQ(ii) = trace(Sig\A);
    trQ_ = trQ(ii);

    term1 = gammaln( (nu_ + n_*p)/2 );
    term2 = - gammaln(nu_/2);
    term3 = - mvgammaln(n_/2, p);
    term4 = nu_/2*log(nu_);
    term5 = - n_/2*logdet(Sig);    
    log_norm_const = term1 + term2 + term3 + term4 + term5;
    
    term6 = (n_-p-1)/2*logdet(A);
    term7 = - (nu_ + n_*p)/2*log(nu_+trQ_);
    log_kernel = term6 + term7;    
    
    logLcontr(ii) = log_norm_const + log_kernel;
end
nLogL = -sum(logLcontr);
%% Gradient (Optional Output)
if nargout >= 3
    
    gradient.Sigma_ = NaN(N,p_);
    gradient.cholSigma = NaN;
    gradient.n = NaN;
    gradient.nu = NaN;
    gradient.all = NaN;
    gradient.Sigma_scaledByInvFisher = NaN(p,p,N);
    
    for ii = 1:N        
        A = X(:,:,ii);
        Sig = Sigma_(:,:,ii);
        invSig = inv(Sig);
        n_ = n(ii);
        nu_ = nu(ii); 
        trQ_ = trQ(ii);
       
        % General matrix derivative (ignoring symmetry of Sigma_): [ enter to matrixcalculus.org: -dfn/2*log(det(Sigma)) - (dft + dfn*p)/2*log(dft+tr(inv(Sigma)*A)) ]        
        score_Sig = - n_*invSig ...
            + (nu_ + p*n_) / (nu_ + trace(Sig\A)) * (Sig\A/Sig);
        score_Sig = .5*score_Sig;
        
        % Accounting for symmetry of Sigma_:
        score_vechSig = 2*score_Sig - diag(diag(score_Sig));
        
        % Numerical inaccuracies make S not exactly symmetric
        score_vechSig = vech((score_vechSig+score_vechSig')./2);
        
        
        if centering 
            gradient.Sigma_(ii,:) = score_vechSig./c(ii);
        else
            gradient.Sigma_(ii,:) = score_vechSig;
        end

        % Wolframalpha querie: d/da log(gamma((a + p*n)/2)) - p*n/2*log(a)
        % - log(gamma(a/2)) - (a + p*n)/2*log(1+q/a) 
        term1 = trQ_*(nu_ + n_*p)/nu_^2/(trQ_/nu_+1);
        term2 = -n_*p/nu_;
        term3 = psi(.5*(nu_+n_*p));
        term4 = log(trQ_/nu_+1);
        term5 = psi(nu_/2);
        gradient.nu(ii) = .5*(term1+term2+term3+term4+term5);
        
        gradient.n(ii) = NaN;
    
    end
    
    varargout{1} = gradient;
    
end
%% Hessian (Optional Output)
if nargout >= 4
    
    hessian.Sigma_ = NaN;
    hessian.cholSigma = NaN;
    hessian.n = NaN;
    hessian.nu = NaN;
    hessian.all = NaN;    
    
    varargout{2} = hessian;
end
%% Fisher Info (Optional Output)
if nargout >= 6
    
    G = Dmatrix(p);
    c1 = n_/2*(nu_+p*n_)/(nu_+p*n_+2);
    c2 = -n_^2/2/(nu_+p*n_+2);
    fisherInfo.Sigma_ = G'*(c1*kron2(invSig) + c2*vec2(invSig))*G;
    
    varargout{4} = fisherInfo;
end
%% Optional Parameter Vector Output
if nargout >= 5    
    varargout{3} = param;
end
end
