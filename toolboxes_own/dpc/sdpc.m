function [ eparam, logL, eS, eL, e_d, enu ] = sdpc( data, p, q, x0 )
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
    x0 = [vech(chol(mean(data,3),'lower'))' ones(1,2+(p+q)*n)*.1];
end
%% Restrictions
AA = [zeros(2+(p+q)*n,n*(n+1)/2) -eye(2+(p+q)*n); % all paramters >= 0 except Sc; actually at least one a_i and b_i is required to be strictly >0, but ok... see HVMA p.51
      zeros(1,n*(n+1)/2) 1 1 zeros(1,(p+q)*n); % stationarity (a_i+b_i <=0) Assumption 2.6.1 DPC paper Aielli and Caporin (2016)
      zeros(n,n*(n+1)/2+2) repmat(eye(n),1,p+q)]; % stationarity (sum(alpha_i's+beta_i's) <=0) see HVMA p.51
bb = [zeros(1,2+(p+q)*n) ones(1,1+n)-1.1e-3];
%% Minimization
% patternsearch for 5 min prior to local minimum search with fmincon -
% update x0 so to say
% % options = optimoptions('patternsearch','Display','iter','UseParallel','always','MaxTime',60); %,'MaxTime',300,'InitialMeshSize',.1,'UseCompletePoll',true,'UseCompleteSearch',true,'MaxMeshSize',.5,'AccelerateMesh',true);
% % eparam = patternsearch(@(param)sdpcLikeRec( param, data, p, q ), x0, AA, bb,[],[],[],[],[], options);
% find local minimum w fmincon
options = optimset('Display','iter','MaxFunEvals',1e10,'MaxIter',1e10); %,'UseParallel','always'
eparam = fmincon(@(param)sdpcLikeRec( param, data, p, q ), x0, AA, bb,[],[],[],[],[], options);
% if sdpcLikeRec(x0,data,p,q)<sdpcLikeRec(eparam,data,p,q)
%     warning('fmincon did not minimize function further. output is patternsearch result')
%     eparam=x0;
% end

options = optimoptions('fmincon','Display','iter');
[~,LogLcontr_2ndpart,eS,eL,e_d] = sdpcLikeRec( eparam, data, p, q );
[enu, nLogL] = fmincon(@(nu)cawlike_givn2ndpart( nu, data, -LogLcontr_2ndpart ), 2*n, -1, 0,[],[],[],[],[], options);
logL=-nLogL;
end