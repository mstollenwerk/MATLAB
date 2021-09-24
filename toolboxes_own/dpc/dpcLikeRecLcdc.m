function [ nLogL, nLogLcontr, S, L, d ] = dpcLikeRecLcdc( params, data, p, q )
%DPCCAWREC calcualtes recursions of DPCCAW model, plugs these into the
%dpccaw-log-likelihood function and returns the negative log-like value
%as well as the recursions.

[n,~,T] = size(data);
dataCell = mat2cell(data(:,:),n,ones(1,T)*n);
%% Parameter Transformations
paramLoad = params(1:n^2+3*n);
paramComp = params(n^2+3*n+1:end);
Lc=reshape(params(1:n^2),n,n);
dc=params(n^2+1:n^2+n)';
%% Loadings Recursion
Q = rbekkRecLcdc(paramLoad, dataCell);
%% Create Components Recursion Input and store L_t's
diagLRL = NaN(n,T);
L = NaN(n,n,T);
for t=1:T
    [L_,~]=dpceig(Q{t},Lc);
    L(:,:,t)=L_;
    diagLRL(:,t)=diag(L_'*dataCell{t}*L_); % diagLRL is not in cell format because so far garchRec requires double input format. In future this might be changed.
end
%% Components Recursion
paramComp = reshape(paramComp,n,q+p); % Each paramComp row now corresponds to parameters for one component.
d = repmat(dc,1,T); % Storage initialize d to unconditional mean s.th. all components that are not included in the Factor-DPCCAW are set to the constant proper value (dc).
for i=1:n
    d(i,:) = garchRec([(1-sum(paramComp(i,:)))*dc(i) paramComp(i,:)],diagLRL(i,:),p,q);
end
%% Get S 
S = NaN(n,n,T);
for t=1:T
    S(:,:,t) = L(:,:,t)*diag(d(:,t))*L(:,:,t)';
end
%% Get Likelihood
[ nLogL, nLogLcontr ] = dpcqlike( d, diagLRL );
end
