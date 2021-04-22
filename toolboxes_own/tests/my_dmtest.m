function [DM,pval] = my_dmtest(loss1, loss2)
% Diebold-Mariano Test.
%
% Usage: [DM,pval] = my_dmtest(loss1, loss2)
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
%% Input Checking
if nargin < 2
    error('dmtest:TooFewInputs','At least two arguments are required');
end
if size(loss1) ~= size(loss2)
    error('dmtest:InvalidInput','Vectors should be of equal length');
end
if size(loss1,2) > 1 || size(loss2,2) > 1
    error('dmtest:InvalidInput','Input should have T rows and 1 column');
end
%%
d = loss1-loss2;
d = d(~isnan(d));
T = size(d,1);
dMean = mean(d);
var_d = covnw(d);
% Retrieve the diebold mariano statistic DM ~N(0,1)
    DM = (dMean / sqrt ( (var_d/T) ));    
%P_VALUE is calculated    
    pval = 2*(1-tcdf(abs(DM),T-1));
end
    
