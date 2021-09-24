function [ eparam, logL, eS, eL, e_d, enu ] = sdpc_procedure( data, p, q, x0 )
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
%% Loading Targeting
[Lc,dc]=dpceig(mean(data,3));
%% x0
if isempty(x0)   
    x0 = [Lc(:)' dc' ones(1,2+(p+q)*n)*.1];
end

%% Iterative Optimization Algorithm - It is not completed, it keeps looping through the while loop b/c optimizer tolerances are too low
% Options
options_fmin = optimset('Display','iter','MaxFunEvals',1e10,'MaxIter',1e10); %,'UseParallel','always'
options_pat = optimoptions('patternsearch','Display','iter','UseParallel','always','MaxTime',300); %,'InitialMeshSize',.1,'UseCompletePoll',true,'UseCompleteSearch',true,'MaxMeshSize',.5,'AccelerateMesh',true);

% Minimization - Lc/GARCH parameters
eparamLcGARCH = fminunc(@(eparamLcGARCH)sdpcLikeRecLcdc( [eparamLcGARCH(1:n^2) dc' eparamLcGARCH(n^2+1:end)], data, p, q ), [x0(1:n^2) x0(n^2+n+1:end)], options_fmin);

% Minimization - Eigenvalues
edc = patternsearch(@(edc)sdpcLikeRecLcdc( [eparamLcGARCH(1:n^2) edc eparamLcGARCH(n^2+1:end)], data, p, q ), x0(n^2+1:n^2+n), [], [], [], [], [], [], [], options_pat);

%% Minimization - Degrees of freedom "nu"
options = optimoptions('fmincon','Display','iter');
[~,LogLcontr_2ndpart,eS,eL,e_d] = sdpcLikeRecLcdc( [eparamLcGARCH(1:n^2) edc eparamLcGARCH(n^2+1:end)], data, p, q );
[enu, nLogL] = fmincon(@(nu)cawlike_givn2ndpart( nu, data, -LogLcontr_2ndpart ), 2*n, -1, 0,[],[],[],[],[], options);
logL=-nLogL;
eparam = [eparamLcGARCH(1:n^2) edc eparamLcGARCH(n^2+1:end)];
end