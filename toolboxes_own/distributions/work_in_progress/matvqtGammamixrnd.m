function A = matvqtGammamixrnd( Sigma_, df_1, df_2, lambda_, N, varargin )
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
 
p = size(Sigma_,1);
%% Random Draws

A = NaN(p,p,N);
for ii=1:N
    
%     % This stochastic representation of the mvt is taken from wikipedia.
%     Y = mvnrnd( zeros(1,p), Sigma_, ceil(df_1) );
%     u = sqrt(chi2rnd( df_2, ceil(df_1), 1)./df_2);
%     
%     X = Y.*repmat(u,1,p);
%     
%     A_ = X(1:end-1,:)'*X(1:end-1,:) ...
%          + (df_1 - floor(df_1)).*X(end,:)'*X(end,:);
%     
%     A(:,:,ii) = gamrnd( lambda_, 1/lambda_, 1 ).*A_;

    W = wishrnd( Sigma_, df_1 );
    gam1 = gamrnd( lambda_, 1/lambda_ );
    gam2 = gamrnd( df_2/2, 2/df_2 );
    
    matvqTrnd = W./gam2;
    
    A(:,:,ii) = gam1.*matvqTrnd;
    
end

% for ii=1:N
%     
%     A(:,:,ii) = gamrnd( lambda_, 1/lambda_, 1 ).*matvtWishrnd( Sigma_, df_1, df_2, 1 ); 
%     
% end

A = sym_covmat(A);

end

