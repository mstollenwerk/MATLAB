function [LLF,likelihoods,h, sk, ku, resid] = garchlik(parameters, data, garchtype, errortype, p, q, m, T)
%{
-----------------------------------------------------------------------
 PURPOSE:
 Leon, A., Rubio, G., and Serna, G., (2004) "Autoregressive Conditional
 Volatility, Skewness and Kurtosis"
-----------------------------------------------------------------------
 USAGE:
 [LLF,likelihoods,h sk,ku,resid] = garchlik(parameters, data, p, q, m, T)

 INPUTS:
 parameters:	vector of parameters
 data:          (T x 1) vector of data
 p:             positive scalar integer representing the order of ARCH
 q:             positive scalar integer representing the order of GARCH
 m:             maximum of p and q
 T:             size of data

 OUTPUT:
 LLF:          the value of the Log-likelihood Function
 likelihoods:  vector of log-likelihoods
 h:            vector of conditional variance
 st:           vector of conditional skewness
 ku:           vector of conditional kurtosis
 resid:        vector of residuals
-----------------------------------------------------------------------
 Author:
 Alexandros Gabrielsen, alexgabriel@live.com
 Date:     05/2016
 -----------------------------------------------------------------------
%}

% Verifying that the vector of parameters is a column vector
[r,c] = size(parameters);
if c>r
    parameters = parameters';
    [r,c] = size(parameters);
end

t = (m+1):T;
likelihoods = [];

% Estimation of conditional meanm, variance, skewness and kurtosis 
% given a set of parameters
[mu,h, sk, ku] = garchskcore(parameters,data, p, q, m, T);

% Gram-Charlier Expansion Series
% f = log((1 + (sk.*((data(t)-mu(t))./sqrt(h(t))).^3 - 3*(data(t)-mu(t))./sqrt(h(t)))./6 + (ku.*(((data(t)-mu(t)./sqrt(h(t))).^4-6*((data(t)-mu(t))./sqrt(h(t))).^2+3))./24)));
% likelihoods = -0.5*((log(h(t))) +(((data(t)-mu(t)).^2)./h(t)) + log(2*pi)) + f;  
% However for some skewness and kurtosis parameters the GC pdf becomes negative, 
% therefore another approach is to estimate the following specification proposed by:
% Leon, Rubio and Serna (2004) Autoregressive Conditional Variance, Skewness and Kurtosis and
% Brio, and Níguez, and Perote (2007) Multivariate Gram-Charlier Densities
f = log((1 + (sk(t).*((data(t)-mu(t))./sqrt(h(t))).^3 - 3*(data(t)-mu(t))./sqrt(h(t)))./6 + (ku(t).*(((data(t)-mu(t)./sqrt(h(t))).^4-6*((data(t)-mu(t))./sqrt(h(t))).^2+3))./24)).^2);
g = log(1+ (sk(t).^2)/6 + (ku(t).^2)/24);
likelihoods = -0.5*((log(h(t))) + (((data(t)-mu(t)).^2)./h(t)) + log(2*pi)) + f - g;      
likelihoods = - likelihoods;
LLF=sum(likelihoods);

if isnan(LLF)
    LLF=1e06;
end

resid=data(t)-mu(t);
h=h(t); 
sk=sk(t);
ku=ku(t);
end