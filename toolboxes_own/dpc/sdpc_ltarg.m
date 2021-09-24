function [ eparam, logL, eS, eL, e_d, enu ] = sdpc_ltarg( data, p, q, x0 )
%SDPC estimates the SDPC(-CAW) model
%
% Reference: Golosnoy et al. (2010); Aielli, Gian Piero and Caporin, 
%            Massimiliano, Dynamic Principal Components: A New Class of 
%            Multivariate GARCH Models (February 3, 2015).
%
% INPUTS:
%
% data       - Data Input, i.e. realized cov-matrices
% x0         - Starting point for the minimization
%
% OUTPUTS:
%
% eparam     - estimated parameters
% logL       - value of log-like @ eparas
% eS         - estimated expectations for Cov-Matrices "S"
% eQ         - estimated loadings recursion
% e_d        - estimated components recursion
[n,~,~]=size(data);
%% x0
[~,dc]=dpceig(mean(data,3));
if isempty(x0)   
    x0 = [dc' ones(1,2+(p+q)*n)*.1];
end
%% Restrictions !!!Notice that the restriction of the eigenvalues is not really necessary and actually appears to be harmful.
% eig_desc_restr= [-eye(n-1) zeros(n-1,1)]+[zeros(n-1,1) eye(n-1)];
AA = [-eye(2+(1+p+q)*n); % all paramters >= 0 except Sc; actually at least one a_i and b_i is required to be strictly >0, but ok... see HVMA p.51
      zeros(1,n) 1 1 zeros(1,(p+q)*n); % stationarity (a_i+b_i <=0) Assumption 2.6.1 DPC paper Aielli and Caporin (2016)
      zeros(n,n+2) repmat(eye(n),1,p+q)]; % stationarity (sum(alpha_i's+beta_i's) <=0) see HVMA p.51
%       eig_desc_restr zeros(n-1,2+(p+q)*n)]; % restriction of eigenvalues eigenvalue i>i-1
bb = [zeros(1,2+(1+p+q)*n) ones(1,1+n)*(1-1e-8)];% zeros(1,n-1)];
%% Loading Targeting
[Lc,~]=dpceig(mean(data,3));
%% Minimization
% options = optimset('Display','iter','MaxFunEvals',1e10,'MaxIter',1e10); %,'UseParallel','always'
options = optimoptions('patternsearch','Display','iter','UseParallel','always','MaxTime',30000); %,'InitialMeshSize',.1,'UseCompletePoll',true,'UseCompleteSearch',true,'MaxMeshSize',.5,'AccelerateMesh',true);
% eparam = fmincon(@(param)sdpcLikeRecLcdc( [Lc(:)' param], data, 1, 1, p, q ), x0, AA, bb,[],[],[],[],[], options);
eparam = patternsearch(@(param)sdpcLikeRecLcdc( [Lc(:)' param], data, 1, 1, p, q ), x0, AA, bb,[],[],[],[],[], options);

options = optimoptions('fmincon','Display','iter');
[~,nLogLcontr_2ndpart,eS,eL,e_d] = sdpcLikeRecLcdc( [Lc(:)' eparam], data, 1, 1, p, q );
[enu, nLogL] = fmincon(@(nu)cawlike_givn2ndpart( nu, data, nLogLcontr_2ndpart ), 2*n, -1, 0,[],[],[],[],[], options);
logL=-nLogL;
end