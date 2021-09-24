function [nLogL, logLcontr, expData, S, Q, L, d, diagLRL, fcstData ] =...
    cs_dpcLikeRec( param, M, data, Pl, Ql, Pc, Qc)
%CS_DPCLIKEREC calculates the likelihood of the consisten DPC model
%
% USAGE:
%       [ eM, eparam, logL, eS, eQ, e_d, fcstS, enu ] = cs_dpc( data, Pl, Ql, Pc, Qc, x0 )
%
% INPUTS:
%       param      - ARCH and GARCH parameters of loadings and components
%       M          - Centering Matrix (parameter)
%       data       - Data Input, i.e. realized cov-matrices
%       Pl         - Number of lagged covMats in Loadings Recursion
%       Ql         - Number of lagged Q Matrices in Loadings Recursion
%       Pc         - Number of lagged diagLRL in Component Recursion
%       Qc         - Number of lagged d_t in Component Recursion
%
% OUTPUTS:
%       nLogL      - value of negative log-likelihood 
%
% COMMENTS:
%
%  See also 
%
% REFERENCES:
%   [1] Gribisch and Stollenwerk

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 08.08.2017

[n,~,T] = size(data);
% Centering data
C = chol(M,'lower');
data_cen = NaN(n,n,T);
for t=1:T
    data_cen(:,:,t) = C\data(:,:,t)/C';
end
% Parameter Transformations
paramLoad = param(1:Pl+Ql);
paramComp = param(Pl+Ql+1:end);
paramComp = reshape(paramComp,n,Pc+Qc);
a = paramLoad(1:Pl);
b = paramLoad(Pl+1:Pl+Ql);
intrcpt = (1-sum(a)-sum(b))*eye(n);
% Storage
Q = cell(1,T+1);
L = NaN(n,n,T);
diagLRL = NaN(n,T);
d = ones(n,T);
fcstd = NaN(n,1);
S = NaN(n,n,T);
expData = NaN(n,n,T); 
% Recursion Q
initial = eye(n);
for t=1:T+1
    Q_ = intrcpt;
    for p=1:Pl
        if t-p<1
            Q_ = Q_ + a(p)*initial;
        else
            Q_ = Q_ + a(p)*data_cen(:,:,t-p);
        end
    end
    for q=1:Ql
        if t-q<1
            Q_ = Q_ + b(q)*initial;
        else
            Q_ = Q_ + b(q)*Q{t-q};
        end
    end
    Q{t} = Q_;
end
fcstQ = Q{end};
Q = Q(1:end-1);
% Recursion L/diagLRL
for t=1:T
    [L_,~]=dpceig(Q{t});
    L(:,:,t)=L_;
    diagLRL(:,t)=diag(L_'*data_cen(:,:,t)*L_);
end
fcstL = dpceig(fcstQ);
% Recursion d
for i=1:n
    [d(i,:), fcstd(i)] = garchRec([(1-sum(paramComp(i,:))) paramComp(i,:)],diagLRL(i,:),Pc,Qc);
end
% Recursion S/expData
for t=1:T
    S(:,:,t) = L(:,:,t)*diag(d(:,t))*L(:,:,t)';
    expData(:,:,t) = C*S(:,:,t)*C';
end
fcstS = fcstL*diag(fcstd)*fcstL';
fcstData = C*fcstS*C';
% Likelihood
likecons = log(det(M));    
logLcontr = -.5*( likecons + sum(log(d) + diagLRL./d) );
nLogL = -sum(logLcontr);
end
