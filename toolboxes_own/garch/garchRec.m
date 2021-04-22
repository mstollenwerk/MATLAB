function [ vars, fcst ] = garchRec( params, vardata, P, Q )
%GARCHREC recursion of variances given innovations acc. to garch(p,q)
%   VARS = GARCHREC(PARAMS,DATA,P,Q)  returns a time series of
%   variances according to a GARCH(p,q) specification, given the DATA of
%   innovations and given the parameters PARAMS. The number of lagged
%   innovations is P, lagged variances is Q. The recursion is initialized
%   with the mean of the data.
%
%   [VARS, FCST] = also returns the T+1 forecast from garch recursion.

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 28.08.2016
%% Parameter Transformations
intrcpt = params(1);
a = params(2:1+P);  
b = params(2+P:end);
%% Initialization
initial=mean(vardata);
vars = NaN(1,max(size(vardata))+1);
%% Recursion
for t=1:max(size(vardata))+1
    v = intrcpt;
    for p=1:P
        if t-p<1
            v = v + a(p)*initial;
        else
            v = v + a(p)*vardata(t-p);
        end
    end
    for q=1:Q
        if t-q<1
            v = v + b(q)*initial;
        else
            v = v + b(q)*vars(t-q);
        end
    end
    vars(t) = v;
end
fcst = vars(end);
vars = vars(1:end-1);
end