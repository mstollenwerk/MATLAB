function [VF, SF, KF, VC] = garchskfor(data, p, q, max_forecast)
%{
-----------------------------------------------------------------------
 PURPOSE: 
 Mean, Volatility and Kurtosis Forecasting 
-----------------------------------------------------------------------
 USAGE:
 [MF, VF, SF, KF, MC, VC] = garchskfor(data, resids, ht, st, kt, parameters, model,ar, ma, p, q, max_forecast)
 
 INPUTS:
 data:          a vector of series
 p:             positive scalar integer representing the order of ARCH
 q:             positive scalar integer representing the order of GARCH
 max_forecasts: maximum number of forecasts (i.e. 1-trading months 22 days)
 
 OUTPUTS:
 VF:            a vector of volatility forecasts
 SF:            a vector of skewness forecasts
 KF:            a vector of kurtosis forecasts
 KF:            a vector of cumulative volatility forecasts
-----------------------------------------------------------------------
 Author:
 Alexandros Gabrielsen, alexgabriel@live.com
 Date:     05/2016
-----------------------------------------------------------------------
%}

if nargin == 0 
    error('Data, Residuals, Variance, Skewness, Kurtosis, GARCH Model, GARCH, Maximum Number of Forecasts') 
end

if size(data,2) > 1
   error('Data vector should be a column vector')
end

if (length(p) > 1) | (length(q) > 1) | p < 0 | q < 0
    error('P and Q should be positive scalars')
end

VF=[];
SF=[];
KF=[];

[parameters, stderrors, LLF, ht, sk, ku, resids, summary]  = garchsk(data, p, q);

 % Forecasting Volatility, Skewness and Kurtosis
 VF(1,1) = parameters(2:2+p+q)'*[1; resids(end-(0:(p-1))).^2;ht(end-(0:(q-1)))]; % 1-period ahead forecast
 SF(1,1) = parameters(3+p+q:5+p+q)'*[1;(resids(end-(0:(p-1)))./ht(end-(0:(p-1)))).^3; sk(end-(0:(q-1)))];
 KF(1,1) = parameters(6+p+q:8+p+q)'*[1;(resids(end-(0:(p-1)))./ht(end-(0:(q-1)))).^4; ku(end-(0:(q-1)))];
 
 for i = 2:max_forecast
     VF(i,1) = parameters(2) + ones(1,p+q)*parameters(3:2+p+q)*VF(i-1,1);
     SF(i,1) = parameters(3+p+q)+ ones(1,p+q)*parameters(4+p+q:5+p+q)*SF(i-1,1);
     KF(i,1) = parameters(6+p+q)+ ones(1,p+q)*parameters(7+p+q:8+p+q)*KF(i-1,1);
 end
clear i

% Estimating cumulative measures
VC = sqrt(cumsum(VF)); 

end
