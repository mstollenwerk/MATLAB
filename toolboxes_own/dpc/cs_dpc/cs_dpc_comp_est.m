function [ eparamComp, e_d, fcstd ] = cs_dpc_comp_est( diagLRL, P, Q, x0Comp )
%DPCCAWESTIMATORCOMP estimates components recursion in DPCCAW
%   DPCCAWESTIMATORCOMP 

[N,T] = size(diagLRL);
%% Storage and Reshape x0 for individual deployment
e_d = NaN(N,T);
fcstd = NaN(N,1);
x0Comp = reshape(x0Comp,N,[]); % Each x0Comp row now corresponds to starting values for one component.
%% Minimizations
% options = optimoptions('fminunc','Display','off','Algorithm','quasi-newton');
options = optimoptions('fmincon','Display','off');    
%Restrictions
if P>1 || Q>1
    error('Only P=Q=1 supported in the components recursion so far, because I dunno the restrictions for bigger P/Q.')
end
AAgarch = ones(1,2); % Stationarity 
bbgarch = 1;
lb = zeros(2,1);
%Minimization
eparamComp = NaN(size(x0Comp));
for i=1:N
    eparamComp(i,:) = fmincon(@(params)cawRVqlike(garchRec([(1-sum(params)) params],diagLRL(i,:),P,Q),...
        diagLRL(i,:)), x0Comp(i,:), AAgarch, bbgarch,[],[],lb,[],[],options);
    [e_d(i,:), fcstd(i)] = garchRec([(1-sum(eparamComp(i,:))) eparamComp(i,:)],diagLRL(i,:),P,Q);
end
eparamComp = eparamComp(:)';
end
