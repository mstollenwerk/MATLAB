function [ eparam, logL, eS, eL, e_d, enu ] = sdpc_alter_withoutrestrictions( data, p, q, x0 )
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
    x0 = ones(1,n*(n+1)/2+2+(p+q)*n)*.1;
end
%% Minimization
options_fmin = optimset('Display','iter','MaxFunEvals',1e10,'MaxIter',1e10); %,'UseParallel','always'
options_pat = optimoptions('patternsearch','Display','iter','UseParallel','always','MaxTime',300); %,'InitialMeshSize',.1,'UseCompletePoll',true,'UseCompleteSearch',true,'MaxMeshSize',.5,'AccelerateMesh',true);

v_x0 = -sdpcLikeRec( x0, data, p, q )-1;
x1 = x0;
v_x1 = -sdpcLikeRec( x1, data, p, q );
while (v_x1-v_x0) > 1e-4
    x0 = x1;
    % Minimization - All Parameters w non-gradient based patternsearch
    x1 = patternsearch(@(x1)sdpcLikeRec( x1, data, p, q ),x0, [], [], [], [], [], [], [], options_pat);
    if -sdpcLikeRec( x1, data, p, q ) <= -sdpcLikeRec( x0, data, p, q )
        warning('Full parameter pattersearch optimization did not improve likeval')
        x1 = x0;
    end
    % Minimization - All Parameters w gradient based fminunc
    x1_ = fminunc(@(x1_)sdpcLikeRec( x1_, data, p, q ), x1, options_fmin);
    if -sdpcLikeRec( x1_, data, p, q ) > -sdpcLikeRec( x1, data, p, q )
        x1 = x1_;
    else
        warning('Full parameter pattersearch optimization did not improve likeval')
    end
    v_x0 = -sdpcLikeRec( x0, data, p, q );
    v_x1 = -sdpcLikeRec( x1, data, p, q );
    disp(v_x1 - v_x0)
end
eparam=x1;
options = optimoptions('fmincon','Display','iter');
[~,LogLcontr_2ndpart,eS,eL,e_d] = sdpcLikeRec( eparam, data, p, q );
[enu, nLogL] = fmincon(@(nu)cawlike_givn2ndpart( nu, data, -LogLcontr_2ndpart ), 2*n, -1, 0,[],[],[],[],[], options);
logL=-nLogL;
end