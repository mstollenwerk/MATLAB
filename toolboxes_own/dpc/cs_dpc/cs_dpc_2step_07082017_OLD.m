function [ eparam, eM, nLogL, logLcontr, estExpData, eS, eQ, eL, e_d, e_diagLRL ] = cs_dpccaw_est_07082017( data, x0 )
%CAW_CHOL_EST estimates a CAW model, where the matrix recursion is a scalar
% BEKK with identity matrix intercept, scaled up by pre- and 
% postmultiplying the recursion matrix with the lower- and upper cholesky
% factor of a parameter matrix M.
%
% USAGE:
%  [eparam,eM,nLogL,logLcontr,expData,S] = caw_chol_est(data)
%
% INPUTS:
%   DATA         - K by K by T array of sym pd matrices (e.g. realized covariance)
%
% OUTPUTS:´
%   
% COMMENTS:
%
%  See also 
%
% REFERENCES:
%      [1] Golosnoy, Gribisch and Liesenfeld (2012) - The conditional 
%      autoregressive Wishart model for multivariate stock volatility

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 02.08.2017
%% Step 1
eM = mean(data,3);
%% Optimization
if isempty(x0)
    x0 = [.05 .92 .05 .05 .05 .92 .92 .92];
end
options = optimoptions('patternsearch','Display','iter'); %,'StepTolerance',1e-5
eparam = patternsearch(@(param)likerec( param, eM, data ), x0,...
    [1 1 0 0 0 0 0 0; zeros(3,2) eye(3) eye(3)], ones(4,1),[],[],...
    zeros(8,1),ones(8,1),[],options );

[ nLogL, logLcontr, estExpData, eS, eQ, eL, e_d, e_diagLRL ] = likerec( eparam, eM, data);
%% Likelihood
function [ nLogL, logLcontr, expData, S, Q, L, d, diagLRL ] = likerec( param, M, data )
    [n,~,T] = size(data);
    C = chol(M,'lower');
    data_cen = NaN(n,n,T);
    for t=1:T
        data_cen(:,:,t) = C\data(:,:,t)/C';
    end
    a = param(1);
    b = param(2);
    alphas = param(3 : 2+n);
    betas = param(2+n+1 : 2+n+n);
    % Initialization
    Q = NaN(n,n,T);    
    L = NaN(n,n,T);
    diagLRL = NaN(n,T);
    d = NaN(n,T);
    S = NaN(n,n,T);
    expData = NaN(n,n,T);    
    eye_n = eye(n);    
    % Initialization Q
    Q(:,:,1) = eye_n;
    L(:,:,1) = eye_n;
    diagLRL(:,1) = diag(data_cen(:,:,1));
    % Recursion Q
    for t=2:T
        Q(:,:,t) = (1-a-b)*eye_n + a*data_cen(:,:,t-1) + b*Q(:,:,t-1);
        [L(:,:,t),~] = dpceig(Q(:,:,t));
        diagLRL(:,t) = diag(L(:,:,t)'*data_cen(:,:,t)*L(:,:,t));        
    end
    % Initialization diagLRL
    for i=1:n
        d(i,1) = (1-alphas(i)-betas(i)) + (alphas(i) + betas(i))*mean(diagLRL(i,:),2);
    end
    % Recursion d 
    for t=2:T
        for i=1:n
            d(i,t) = (1-alphas(i)-betas(i)) + alphas(i)*diagLRL(i,t-1) + betas(i)*d(i,t-1);
        end
    end
    % Recursion S/expData
    for t=1:T % Starting from 1!
        S(:,:,t) = L(:,:,t)*diag(d(:,t))*L(:,:,t)';
        expData(:,:,t) = C*S(:,:,t)*C';
    end
    % Likelihood
    likecons = log(det(M));    
	logLcontr = -.5*( likecons + sum(log(d) + diagLRL./d) );
    nLogL = -sum(logLcontr);
end

end
