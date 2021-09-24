function [ eparam, logL, eS, eL, e_d, enu ] = sdpc_ftarg( data, p, q, Lc, dc, x0 )
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
if isempty(x0)   
    x0 = ones(1,2+(p+q)*n)*.1;
end
%% Restrictions 
AA = [-eye(2+(p+q)*n); % all paramters >= 0 except Sc; actually at least one a_i and b_i is required to be strictly >0, but ok... see HVMA p.51
      1 1 zeros(1,(p+q)*n); % stationarity (a_i+b_i <=0) Assumption 2.6.1 DPC paper Aielli and Caporin (2016)
      zeros(n,2) repmat(eye(n),1,p+q)]; % stationarity (sum(alpha_i's+beta_i's) <=0) see HVMA p.51
bb = [zeros(1,2+(p+q)*n) ones(1,1+n)*(1-1e-8)];
%% Minimization
options = optimset('Display','off','MaxFunEvals',1e10,'MaxIter',1e10); %,'UseParallel','always'
% options = optimoptions('patternsearch','Display','iter','UseParallel','always','MaxTime',300); %,'InitialMeshSize',.1,'UseCompletePoll',true,'UseCompleteSearch',true,'MaxMeshSize',.5,'AccelerateMesh',true);
eparam = fmincon(@(param)sdpcLikeRecLcdc( [Lc(:)' dc' param], data, p, q ), x0, AA, bb,[],[],[],[],[], options);
% eparam = patternsearch(@(param)sdpcLikeRecLcdc( [Lc(:)' param], data, p, q ), x0, AA, bb,[],[],[],[],[], options);

options = optimoptions('fmincon','Display','off');
[~,LogLcontr_2ndpart,eS,eL,e_d] = sdpcLikeRecLcdc( [Lc(:)' dc' eparam], data, p, q );
[enu, nLogL] = fmincon(@(nu)cawlike_givn2ndpart( nu, data, -LogLcontr_2ndpart ), 2*n, -1, 0,[],[],[],[],[], options);
logL=-nLogL;
end