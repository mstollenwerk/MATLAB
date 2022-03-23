function [ nLogL, logLcontr, varargout ] = ...
    matvsFRieszlike( Sigma_, n, nu, X, varargin )
%MATVSFRIESZLIKE
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 11.02.2021

[p,~,N] = size(X);
I = eye(p);
p_ = p*(p+1)/2;
narginchk(4,5); %%%%%%%
nargoutchk(0,6);
%% Param
if nargin == 5 %%%%%%%
    if ~(isempty(Sigma_) && isempty(n) && isempty(nu))
        error('Cannot input all_param and any of the parameters individually!')
    end
    all_param = varargin{1};
    Sigma_ = ivechchol(all_param(1 : p_));
    n = all_param(p_ + 1 : p_ + p);
    nu = all_param(p_ + p + 1 : p_ + p + p);
end
% Checking if Sigma_ is symmetric p.d.
param.Sigma_ = Sigma_;
C = chol(Sigma_,'lower');
param.chol_Sigma_ = vech(C);
param.n = n;
param.nu = nu;
diaga = matvFRieszexpmat(n,nu);
param.all = [param.chol_Sigma_; n; nu];
%% Log-likelihood computation
logLcontr = NaN(N,1);

term1 = sum(n.*log(diag(diaga)))/2;
term2 = lgmvgammaln(flipud(n + nu)./2);%ugmvgammaln((n + nu)./2);
term3 = -lgmvgammaln(flipud(nu)./2);%-ugmvgammaln(nu./2);
term4 = -lgmvgammaln(n./2);
term5 = -loglpwdet([],nu./2,diag(C));

log_normalizing_constant = term1 + term2 + term3 + term4 + term5;

C_B = NaN(p,p,N);
for ii = 1:N
    
    R = X(:,:,ii);
    Cr = chol(R,'lower');
    B = I + C'\diaga/C*R;
    C_B(:,:,ii) = chol(B,'lower');
    
    term6 = loglpwdet([],(n-p-1)./2,diag(Cr));
    term7 = loglpwdet([], -(n+nu)./2, diag(C_B(:,:,ii)));
    
    log_kernel = term6 + term7;

    logLcontr(ii) = log_normalizing_constant + log_kernel;
end
nLogL = -sum(logLcontr);
%% Score computation
if nargout >= 3
    
    score.Sigma_ = NaN(p,p,N);    
    for ii = 1:N
        
        % General matrix derivative (ignoring symmetry of Sigma_):
        Nabla = - C'\trilHalfDiag(C'*tril( ...
            C'\diag(n) + C_B(:,:,ii)'\diag(nu+n)/C_B(:,:,ii)*C/diaga...
            ))/C; 

        % Accounting for symmetry of Sigma_:
        score.Sigma_(:,:,ii) = Nabla+Nabla' - diag(diag(Nabla));

    end
    
    varargout{1} = score;

end
%% Hessian (Optional Output)
%% Fisher Info (Optional Output)
if nargout >= 6
   
    %%% I don't have the exact closed form yet, there are some possible
    %%% approximation candidates.
    
    % 1. Approximate with matrix-F fisherinfo, replacing n and nu with the
    % mean(n) and mean(nu). Remember that the FRiesz becomes the matrix-F
    % if all n are equal and all nu are equal. This seems the fairest
    % approximation.
%     n_ = mean(n);
%     nu_ = mean(nu);    
%     c_1 = (n_^2*(nu_-p-2) + 2*n_)/((nu_-p)*(nu_-p-1)*(nu_-p-3));
%     c_2 = (n_*(nu_-p-2)+n_^2+n_)/((nu_-p)*(nu_-p-1)*(nu_-p-3));
%     c_4 = (n_-p-1)/((n_+nu_-1)*(n_+nu_+2))*((n_-p-2+1/(nu_+n_))*c_2-(1+(n_-p-1)/(n_+nu_))*c_1);
%     c_3 = (n_-p-1)/(n_+nu_)*((n_-p-2)*c_2 - c_1)-(n_+nu_+1)*c_4;
%     invSigF = inv(Sigma_);
%     
%     G = Dmatrix(p);
%     conskron2 = (nu_-(n_+nu_)*(c_3+c_4));
%     consvec2 = (n_+nu_)*c_4;
%     fishSig = 1/2*G'*(conskron2*kron2(invSigF) - consvec2*vec2(invSigF))*G;
    
     
%     invSig = inv(Sigma_);
%     [G, iG] = Dmatrix(p);
%     L = ELmatrix(p);
%     I = speye(p);
%     iC = inv(C);

      % 2. Approximate with part Riesz, part matrix-F fisherinfo. This is
      % motivated by the fact that there are two summands in the FRiesz
      % fisherinfo, the first one we have derived, it is equal to the Riesz
      % fisherinfo, the second one we have not derived its expectation yet.
      % The matrix-F has the same structure, but for the matrix-F we have
      % the expectation of the second summand as well. So we glue the first
      % summand of the FRiezs fisherinfo with the second summand of the
      % matrix-F fisherinfo. This seems a bit aribtrary and does not work
      % well in practice (measured by the fit).
%     n_ = mean(n);
%     nu_ = mean(nu); 
%     c_1 = (n_^2*(nu_-p-2) + 2*n_)/((nu_-p)*(nu_-p-1)*(nu_-p-3));
%     c_2 = (n_*(nu_-p-2)+n_^2+n_)/((nu_-p)*(nu_-p-1)*(nu_-p-3));
%     c_4 = (n_-p-1)/((n_+nu_-1)*(n_+nu_+2))*((n_-p-2+1/(nu_+n_))*c_2-(1+(n_-p-1)/(n_+nu_))*c_1);
%     c_3 = (n_-p-1)/(n_+nu_)*((n_-p-2)*c_2 - c_1)-(n_+nu_+1)*c_4;
% 
%     term_Riesz = 1/2*G'*kron(C'\diag(n),invSig)*L'/(iG*kron(C,I)*L');
%     term_Rest = -(n_ + nu_)/2*G'*((c_3+c_4)*kron2(invSig) + c_4*vec2(invSig))*G;
%     
%     fishSig = term_Riesz + term_Rest;
    
      % 3. Approximate by replacing the expectation of the second summand
      % of the hessian (See 2. Approximation for explanation) hessian,
      % where A is replaced by expectation of A ( so basically do f(E[A]) 
      % instead of E[f(A)] ). This provides good fit, but gives completely
      % out of the range paramter values for the dynamic parameters (~ -3
      % for the score parameter), which are not comparable to other GAS
      % models
%     fishSig = .5*G'*kron(iC',C'\diag(nu)/C)*L'/(G'*kron(C,I)*L')*(G'*G) ...
%               - G'*kron(iC',C'\diag(nu + n)/C)*L'/(G'*kron(C,I)*L')*(G'*G);
          
%     % 4. Approximate by Riesz fisherinfo since for nu->infinity, the
%     FRiesz goes to the Riesz.
%     fishSig = 1/2*G'*kron(C'\diag(n),invSig)*L'/(iG*kron(C,I)*L');
    
%     % 5. Approximate by Wishart fisherinfo since for nu->infinity, all n
%     equal and all nu equal the FRiezs converges to the Wishart. Replace n
%     by mean(n).
%     fishSig = mean(n)/2*G'*kron2(invSig)*G;
    
%     fisherinfo.Sigma_ = fishSig;
%     fisherinfo.df_1 = NaN;
%     fisherinfo.df_2 = NaN;
%     
%     varargout{4} = fisherinfo;
% end
%% Optional Parameter Vector Output
if nargout >= 5
    varargout{3} = param;
end
end
