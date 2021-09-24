function [ nLogL, nLogLcontr, S, L, d ] = sdpcLikeRecLcdc( params, data, lP, lQ, cP, cQ )
%DPCCAWREC calcualtes recursions of DPCCAW model, plugs these into the
%dpccaw-log-likelihood function and returns the negative log-like value
%as well as the recursions.

%% Input Checking
if isnumeric(data)
    [n,~,T] = size(data);
    data = mat2cell(data(:,:),n,ones(1,T)*n);
elseif iscell(data)
    n = size(data{1},1);
    T = max(size(data));
else
    error('sdcpLikeRecLcdc input data must have cell or 3-dim array format')
end
%% Parameter Transformations
Lc=reshape(params(1:n^2),n,n);
dc=sort(params(n^2+1:n^2+n),'descend');
paramLoad = [params(1:n^2) dc params(n^2+n+1:n^2+n+lP+lQ)];
paramComp = params(n^2+n+lP+lQ+1:end);
%% Loadings Recursion
lQ = sbekkRecLcdc(paramLoad, lP, lQ, data);
%% Create Components Recursion Input and store L_t's
diagLRL = NaN(n,T);
L = NaN(n,n,T);
for t=1:T
    try % Very dirty fix for when Q recursion has nan or inf values due to non-stationary garch parameters.
        [L_,~]=dpceig(lQ{t},Lc);
    catch 
        nLogL = 1e7;
        return
    end
    L(:,:,t)=L_;
    diagLRL(:,t)=diag(L_'*data{t}*L_); % diagLRL is not in cell format because so far garchRec requires double input format. In future this might be changed.
end
%% Components Recursion
paramComp = reshape(paramComp,n,cQ+cP); % Each paramComp row now corresponds to parameters for one component.
d = repmat(dc',1,T); % Storage initialize d to unconditional mean s.th. all components that are not included in the Factor-DPCCAW are set to the constant proper value (dc).
for i=1:n
    d(i,:) = garchRec([(1-sum(paramComp(i,:)))*dc(i) paramComp(i,:)],diagLRL(i,:),cP,cQ);
end
%% Get S 
S = NaN(n,n,T);
for t=1:T
    S(:,:,t) = L(:,:,t)*diag(d(:,t))*L(:,:,t)';
end
%% Get Likelihood
[ nLogL, nLogLcontr ] = dpcqlike( d, diagLRL );
end
