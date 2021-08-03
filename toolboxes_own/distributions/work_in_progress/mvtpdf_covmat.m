function [ y ] = mvtpdf_covmat( x, mu, cov_mat, df_t )
%MVTPDF_COVMAT give value of multivariate t-distribution pdf.
%
% Whats special is the input of the actual covariance matrix. Also since I
% ommited the pd checks this code is faster than the matlab provided
% mvtpdf.
%
% USAGE:
%   Y = mvtlike(x,mu,cov_mat,df_t)
%
% INPUTS:
%   X            - N by K array data points
%   MU           - N by K array of mean parameters.
%   COV_MAT      - K by K by N array of covariance matrices. 
%   DF_T         - Degree of freedom (dof) of the underlying multivariate 
%                  t-distribution
%
% OUTPUTS:
%   Y            - Value of the multivariate t-distribution pdf at X.
%
% COMMENTS:
%
%  See also 
%
% REFERENCES:
%      [1] 

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 26.06.2018

%warning('Not properly debugged yet')

%% Error checking
narginchk(4,4);

if df_t<=2
    warning('Degrees of freedom are smaller than 2. Covariance Matrix does not exist')
elseif df_t<=1
    warning('Degrees of freedom are smaller than 1. Mean and Covariance Matrix do not exist')
end
[k,k2,n] = size(cov_mat);
if k~=k2
    error('Covariance matrix parameters of multivariate t-distribution must be square.')
end

if size(x,1)~=n  || size(x,2)~=k
    error('Data and covariance parameter matrices must have appropriate dimensions')
end
%% Log-likelihood computation
x_c = x-mu; % center data
y=NaN(n,1);

if isnumeric(cov_mat) % Double input
    for i=1:n
        y(i)= gamma((df_t+k)/2)/gamma(df_t/2)/((df_t-2)*pi)^(k/2)/det(cov_mat(:,:,i))^.5*...
            (1+x_c(i,:)/cov_mat(:,:,i)*x_c(i,:)'/(df_t-2))^-((df_t+k)/2);
    end
elseif iscell(cov_mat) % Cell input
    for i=1:n
        y(i)= gamma((df_t+k)/2)/gamma(df_t/2)/((df_t-2)*pi)^(k/2)/det(cov_mat{i})^.5*...
            (1+x_c(i,:)/cov_mat{i}*x_c(i,:)'/(df_t-2))^-((df_t+k)/2);
    end
end

end
