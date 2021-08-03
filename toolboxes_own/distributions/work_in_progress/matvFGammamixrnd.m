function A = matvFGammamixrnd( Sigma_, df_1, df_2, lambda_, N, varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Parameters
if nargin == 7
    all_param = varargin{1};
    p = varargin{2};
    p_ = p*(p+1)/2;
    Sigma_ = ivechchol(all_param(1:p_));
    df_1 = all_param(p_ + 1);
    df_2 = all_param(p_ + 2);
    lambda_ = all_param(p_ + 3);
end

%% Random Draws
sqrtmSig = sqrtm(Sigma_);

A = NaN(p,p,N);
for ii=1:N
    
    V = wishrnd( eye(p), df_1 );
    sqrtmV = sqrtm( V );
    
    invT = iwishrnd( eye(p) , df_2);
       
    gamma_ = gamrnd( lambda_, 1/lambda_, 1 );
    
    A(:,:,ii) = gamma_.*sqrtmSig*sqrtmV*invT*sqrtmV*sqrtmSig;    
    
end

A = sym_covmat(A);
end

