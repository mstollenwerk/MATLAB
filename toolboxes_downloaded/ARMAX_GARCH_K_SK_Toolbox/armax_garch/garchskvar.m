function VaR = garchskvar(data, p, q, max_forecast, alpha)
%{
-----------------------------------------------------------------------
 PURPOSE: 
 Value-at-Risk Estimation for both long and short positions (Model Estimation)
-----------------------------------------------------------------------
 USAGE:
 VaR = garchvar(data, model, ar, ma, p, q, max_forecast, alpha)
 
 INPUTS:
 data:          (T x 1) vector of data
 p:             positive scalar integer representing the order of ARCH
 q:             positive scalar integer representing the order of GARCH
 max_forecasts: maximum number of forecasts (i.e. 1-trading months 22 days)
 alpha:         confidence level
 
 OUTPUTS:
 VaR:           a vector of VaR forecasts
-----------------------------------------------------------------------
 Author:
 Alexandros Gabrielsen, alexgabriel@live.com
 Date:     05/2016
-----------------------------------------------------------------------
%}

if nargin == 0 
    error('Data,GARCH, Maximum Number of Forecasts, Confidence Level') 
end

if size(data,2) > 1
   error('Data vector should be a column vector')
end

if (length(p) > 1) | (length(q) > 1) | p < 0 | q < 0 | alpha < 0 | alpha > 1
    error('P,Q and alpha should be positive scalars')
end

% Estimate the model
[parameters, stderrors, LLF, ht, sk, ku, resids, summary] = garchsk(data, p, q);

% Calling garchkfor2 function to pass the parameters
[VF, SF, KF, VC] = garchskfor2(data, resids, ht, sk, ku, parameters,p, q, max_forecast);

% Estimating Inverce-CDF as a function of the forecasted conditional skewness, kurtosis
cdfqtile1 = norminv(alpha,zeros(size(VF,1),1),ones(size(VF,1),1)).*(1+(SF/6).*(norminv(alpha,zeros(size(VF,1),1),ones(size(VF,1),1)).^2-1)+((KF-3)/24).*(norminv(alpha,zeros(size(VF,1),1),ones(size(VF,1),1)).^3 - 3*norminv(alpha,zeros(size(VF,1),1),ones(size(VF,1),1))));
cdfqtile2 = norminv(1-alpha,zeros(size(VF,1),1),ones(size(VF,1),1)).*(1+(SF/6).*(norminv(alpha,zeros(size(VF,1),1),ones(size(VF,1),1)).^2-1)+((KF-3)/24).*(norminv(alpha,zeros(size(VF,1),1),ones(size(VF,1),1)).^3 - 3*norminv(alpha,zeros(size(VF,1),1),ones(size(VF,1),1))));

% Estimating VaR
VaR = [ cdfqtile1.*VC, cdfqtile2.*VC];

end
