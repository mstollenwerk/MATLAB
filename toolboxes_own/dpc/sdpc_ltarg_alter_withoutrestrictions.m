function [ eparam, logL, eS, eL, e_d, enu ] = sdpc_ltarg_alter_withoutrestrictions( data, p, q, x0 )
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
    x0 = [dc' ones(1,2+(p+q)*n)*.1];
end

%% Iterative Optimization Algorithm - It is not completed, it keeps looping through the while loop b/c optimizer tolerances are too low
% Options
options_fmin = optimset('Display','iter','MaxFunEvals',1e10,'MaxIter',1e10); %,'UseParallel','always'
options_pat = optimoptions('patternsearch','Display','iter','UseParallel','always','MaxTime',300); %,'InitialMeshSize',.1,'UseCompletePoll',true,'UseCompleteSearch',true,'MaxMeshSize',.5,'AccelerateMesh',true);

v_x0 = -sdpcLikeRecLcdc( [Lc(:)' x0], data, p, q )-1;
x1 = x0;
v_x1 = -sdpcLikeRecLcdc( [Lc(:)' x1], data, p, q );
% while (v_x1-v_x0) > 1e-4
    x0 = x1;
    % Minimization - GARCH parameters
    x1_garch = fminunc(@(x1_garch)sdpcLikeRecLcdc( [Lc(:)' x0(1:n) x1_garch], data, p, q ), x0(n+1:end), options_fmin);
    if -sdpcLikeRecLcdc( [Lc(:)' x0(1:n) x1_garch], data, p, q ) <= -sdpcLikeRecLcdc( [Lc(:)' x0], data, p, q )
        warning('GARCH parameter optimization did not improve likeval')
        x1_garch = x0(n+1:end);
    end
    % Minimization - Eigenvalues
    x1_eig = patternsearch(@(x1_eig)sdpcLikeRecLcdc( [Lc(:)' x1_eig x1_garch], data, p, q ), x0(1:n), [], [], [], [], [], [], [], options_pat);
    if -sdpcLikeRecLcdc( [Lc(:)' x1_eig x1_garch], data, p, q ) <= -sdpcLikeRecLcdc( [Lc(:)' x0(1:n) x1_garch], data, p, q )
        warning('Eigenvalue parameter optimization did not improve likeval')
        x1_eig = x0(1:n);
    end
%     % Minimization - All Parameters w gradient based fminunc
%     x1_fmin = fminunc(@(x1_fmin)sdpcLikeRecLcdc( [Lc(:)' x1_fmin], data, p, q ), [x1_eig x1_garch], options_fmin);
%     if -sdpcLikeRecLcdc( [Lc(:)' x1_fmin], data, p, q ) <= -sdpcLikeRecLcdc( [Lc(:)' [x1_eig x1_garch]], data, p, q )
%         warning('Full parameter fminunc optimization did not improve likeval')
        x1_fmin = [x1_eig x1_garch];
%     end
%     % Minimization - All Parameters w non-gradient based patternsearch
%     x1 = patternsearch(@(x1)sdpcLikeRecLcdc( [Lc(:)' x1], data, p, q ),x1_fmin, [], [], [], [], [], [], [], options_pat);
%     if -sdpcLikeRecLcdc( [Lc(:)' x1], data, p, q ) <= -sdpcLikeRecLcdc( [Lc(:)' x1_fmin], data, p, q )
%         warning('Full parameter pattersearch optimization did not improve likeval')
        x1 = x1_fmin;
%     end
%     v_x0 = -sdpcLikeRecLcdc( [Lc(:)' x0], data, p, q );
%     v_x1 = -sdpcLikeRecLcdc( [Lc(:)' x1], data, p, q );
%     disp(v_x1 - v_x0)
% end

%% Minimization - Degrees of freedom "nu"
options = optimoptions('fmincon','Display','iter');
[~,LogLcontr_2ndpart,eS,eL,e_d] = sdpcLikeRecLcdc( [Lc(:)' x1], data, p, q );
[enu, nLogL] = fmincon(@(nu)cawlike_givn2ndpart( nu, data, -LogLcontr_2ndpart ), 2*n, -1, 0,[],[],[],[],[], options);
logL=-nLogL;
eparam = x1;
end