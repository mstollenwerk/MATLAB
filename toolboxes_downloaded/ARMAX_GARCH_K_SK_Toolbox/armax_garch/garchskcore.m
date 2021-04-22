function [mu,h, sk, ku] = garchcore(parameters, data, p, q, m, T)
%-----------------------------------------------------------------------
% PURPOSE:
% Leon, A., Rubio, G., and Serna, G., (2004) "Autoregressive Conditional
% Volatility, Skewness and Kurtosis"
%-----------------------------------------------------------------------
% USAGE:
% [mu,h, sk, ku] = garchcore(parameters, data, p, q, m, T)
%
% INPUTS:
% parameters:	vector of parameters
% data:         (T x 1) vector of data
% p:            positive scalar integer representing the order of ARCH
% q:            positive scalar integer representing the order of GARCH
% m:            maximum of p and q
% T:            size of data
%
% OUTPUTS:
% mu:           conditional mean
% h:            conditional variance
% st:           vector of conditional skewness
% ku:           vector of conditional kurtosis
%-----------------------------------------------------------------------
% Author:
% Alexandros Gabrielsen, alexgabriel@live.com
% Date:     05/2016
%-----------------------------------------------------------------------

% Verifying that the vector of parameters is a column vector
[r,c] = size(parameters);
if c>r
    parameters = parameters';
    [r,c] = size(parameters);
end

% Initial parameters
mu = [];
h = [];
sk=[];
ku=[];
mu(1:m,1) = parameters(1);
h(1:m,1) = var(data); % another way is to set the initial value to the conditional
sk(1:m,1) = 0.5;
ku(1:m,1) = 10;
% Estimation of Conditional Mean and Variance
for t = (m+1):T; 
   mu(t,1) = parameters(1);
   % original conditional variance estimation for the Gaussian distribution
   % h(t,1) = parameters'*[zeros(size(parameters,1)-1,1) eye(size(parameters,1)-1)]'*[1;(data(t-(1:p)) - mu(t-(1:p))).^2; h(t-(1:q))]; 
   h(t,1) = parameters(2:2+p+q)'*[1;(data(t-(1:p)) - mu(t-1)).^2; h(t-(1:q))];
   sk(t,1) = parameters(3+p+q:5+p+q)'*[1;((data(t-(1:p)) - mu(t-1))./h(t-(1:p))).^3; sk(t-(1:q))];
   ku(t,1) = parameters(6+p+q:8+p+q)'*[1;((data(t-(1:p)) - mu(t-1))./h(t-(1:p))).^4; ku(t-(1:q))];
end


end
