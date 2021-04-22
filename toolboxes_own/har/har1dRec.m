function [ vars, fcst ] = har1dRec( param, vardata )
%HARREC1D recursion of variances acc. to HAR process
%   VARS = HARREC1D(PARAMS,DATA)  returns a time series of
%   variances according to a HAR specification, given the DATA of
%   innovations and given the parameters PARAMS. The recursion is 
%   initialized with the mean of the vardata.
%
%   [ VARS, FCST ] = HARREC1D(PARAMS,DATA)  also returns the one period 
%   ahead forcast.

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 29.08.2016  
%% Initialization and Storage
vars = repmat(mean(vardata),1,max(size(vardata)));
%% Recursion
for t=21:max(size(vardata))+1
    vars(t) = param(1) + param(2)*vardata(t-1) + param(3)*mean(vardata(t-5:t-1)) + param(4)*mean(vardata(t-10:t-1)) + param(5)*mean(vardata(t-20:t-1));
end
fcst = vars(end);
vars = vars(1:end-1);
end