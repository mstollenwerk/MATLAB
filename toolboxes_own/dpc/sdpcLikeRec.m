function [ nLogL, nLogLcontr, S, L, d, Q ] = sdpcLikeRec( params, data, p, q )
%
[n,~,T] = size(data);
%% Parameter Transformations
paramLoad = params(1:n*(n+1)/2+2);
paramComp = params(n*(n+1)/2+3:end);
Sc = ivech(params(1:n*(n+1)/2),'lower');
Sc = Sc*Sc';
[Lc,dc]=dpceig(Sc);
%% Loadings Recursion
Q = sbekkRec_s(paramLoad, data, 1, 1);
%% Create Components Recursion Input and store L_t's
diagLRL = NaN(n,T);
L = NaN(n,n,T);
for t=1:T
    try % Very dirty fix for when Q recursion has nan or inf values due to non-stationary garch parameters.
        [L_,~]=dpceig(Q{t},Lc);
    catch 
        nLogL = 1e7;
    end
    L(:,:,t)=L_;
    diagLRL(:,t)=diag(L_'*data(:,:,t)*L_); % diagLRL is not in cell format because so far garchRec requires double input format. In future this might be changed.
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
