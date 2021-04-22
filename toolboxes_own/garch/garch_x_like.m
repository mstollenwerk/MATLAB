function [nLogL, lls, innov, condVar, fcst_condVar] = garch_x_like(param, returndata, X, X_g, P, Q)

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 03.05.2017

[T,~]=size(returndata);
% Parameters
paramMean = param(1:size(X,2));
intrcpt = param(size(X,2)+1);
a = param(size(X,2)+2:size(X,2)+P+1);
b = param(size(X,2)+P+2:size(X,2)+P+Q+1);
c = param(size(X,2)+P+Q+2:size(X,2)+P+Q+size(X_g,2)+1);
% Mean Recursion
innov = NaN(size(returndata));
if ~isempty(X)
    for t=1:T
        innov(t) = returndata(t) - paramMean*X(t,:)';
    end
else
    innov = returndata;
end
% Variance Recursion
% m_innov_sq = mean(innov.^2); % Initialization GARCH condVars.
m_innov_sq = var(innov); % Initialization GARCH condVars.
lls = NaN(T,1);
condVar = NaN(T+1,1);
for t=1:T+1
    v = intrcpt;
    for p=1:P
        if t-p<1
            v = v + a(p)*m_innov_sq; 
        else
            v = v + a(p)*innov(t-p).^2;
        end
    end
    for q=1:Q 
        if t-q<1
            v = v + b(q)*m_innov_sq; 
        else
            v = v + b(q)*condVar(t-q);
        end
    end
    for i=1:size(X_g,2)
            v = v + c(i)*X_g(t,i);
    end         
    condVar(t) = v;
end
fcst_condVar = condVar(T+1);
condVar = condVar(1:end-1);
%
% if any(condVar<0)
%     nLogL = 1e10;
%     return
% end    
%Likelihood Calculation
for t=1:T
    lls(t) = -0.5*(log(2*pi) + log(condVar(t)) + innov(t).^2./condVar(t));
end
nLogL = -sum(lls);

end