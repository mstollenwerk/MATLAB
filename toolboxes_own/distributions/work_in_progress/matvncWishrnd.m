function y = matvncWishrnd( df_, Sigma_, Theta_, varargin)
% as 4th element either enter N for the number of random matrices you want
% or enter Z a df_ceil*p by N array of random standard normal numbers.
df_ceil = ceil(df_);

p = size(Sigma_,1);
% Gupta, Nagar (1999) - Matrix Variate Distributions, p.114
MMt = Sigma_*Theta_;
if all(MMt(:)==0)
    M = zeros(p,df_ceil);
else
    M = [chol(MMt,'lower'), zeros(p,df_ceil-p)]; % M must be p by df_.
end
Mt = M';
vecMt = Mt(:);
% p.55 for matrix normal acc. to Gupta, Nagar.

if nargin == 3
    vecXt = mvnrnd( vecMt, kron( Sigma_, eye(df_ceil) ), 1 );
elseif nargin == 4
    if length(varargin{1})==1
        N = varargin{1};
        vecXt = mvnrnd( vecMt, kron( Sigma_, eye(df_ceil) ), N );
    elseif size(varargin{1},1) == p*df_ceil 
        Z = varargin{1};
        Sigma_large = kron( Sigma_, eye(df_ceil) );
        vecXt = vecMt + sqrtm(Sigma_large)*Z;
        vecXt = vecXt';
        N = size(Z,2);
    else
        error('Varagin must be N (number of random matrices) or Z (p*df_ceil by >1 array of random standard normal numbers.')
    end
else
    error('Nargin must be 4 or 5')
end   

y = NaN(p,p,N);
for ii = 1:N
    X = reshape(vecXt(ii,:),df_ceil,p)';
    % My heuristic way to deal with non-integer dfs.
    if df_ ~= df_ceil
        X(:,end) = X(:,end).*rem(df_,1);
    end
    y(:,:,ii) = X*X';
end

end