function [ eM, eparam, nLogL, logLcontr, e_expData, eS, eQ, eL, e_d, e_diagLRL, fcstData ] =...
    cs_dpc_3step( data, Pl, Ql, Pc, Qc, x0 )
%DPCCAW estimates the consistent DPC-CAW model.
%
% USAGE:
%       [ eM, eparam, logL, eS, eQ, e_d, fcstS, enu ] = cs_dpc( data, Pl, Ql, Pc, Qc, x0 )
%
% INPUTS:
%       data       - Data Input, i.e. realized cov-matrices
%       Pl         - Number of lagged covMats in Loadings Recursion
%       Ql         - Number of lagged Q Matrices in Loadings Recursion
%       Pc         - Number of lagged diagLRL in Component Recursion
%       Qc         - Number of lagged d_t in Component Recursion
%       x0         - [Optional] Starting point for the minimization
%
% OUTPUTS:
%       eM         - estimated M matrix of consistent DPC-CAW (see paper)
%       eparam     - estimated parameters of centered process
%       logL       - value of log-like of centered process @ eparam 
%       eS         - estimated expectations for Cov-Matrices "S"
%       eQ         - estimated loadings recursion matrices
%       e_d        - estimated components recursion
%       fcstS      - Forecast of S in T+1
%       enu        - estimated dof parameter of centered process
%
% COMMENTS:
%
%  See also 
%
% REFERENCES:
%   [1] Gribisch and Stollenwerk

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 12.06.2017
% 08.08.2017: removed helper functions
%             included eM into likelihood

[n,~,~] = size(data);
%% Step 1 - center data
eM = mean(data,3);
%% x0
if isempty(x0)
    x0 = [ones(1,Pl)*.05/Pl ones(1,Ql)*.92/Ql ones(1,Pc*n)*.05/Pc ones(1,Qc*n)*.92/Qc];
end
%% Minimization 
% Restrictions for Loadings Parameters
if Pl>1 || Ql>1
    error('Only P=Q=1 supported in the loadings recursion so far, because I dunno the restrictions for bigger P/Q.')
end
AA = [1 1]; % Stationarity of Loadings Recursion (a+b <=1) 
bb = 1; %*(1-1e-8);
lb = [0 .7];
ub = [.25 1];
% options = optimoptions('fmincon','Display','iter'); %,'UseParallel','always','MaxTime',3600,'MeshContractionFactor',.77,'InitialMeshSize',.1,'UseCompletePoll',true,'UseCompleteSearch',true,'MaxMeshSize',.5,'AccelerateMesh',true);
options = optimoptions('patternsearch','MaxTime',3600);
like_diff = 1;
step_ = 1;
eparam = NaN(size(x0));
while like_diff > 1e-4
% MULTISTART
%     problem = createOptimProblem('fmincon','objective',...
%      @(param)cs_dpcLikeRec([param x0(sum(Pl+Ql)+1:end)],eM,data,Pl,Ql,Pc,Qc),...
%      'x0',x0(1:sum(Pl+Ql)),'Aineq',AA,'bineq',bb,'lb',lb,'options',options);
%     ms = MultiStart;
%     eparam(1:sum(Pl+Ql)) = run(ms,problem,10);

% GRIDSEARCH
%     eparam(1:sum(Pl+Ql)) = ab_gridsearch( .01, @(param)cs_dpcLikeRec([param x0(sum(Pl+Ql)+1:end)],eM,data,Pl,Ql,Pc,Qc) );

% PATTERNSEARCH
%     eparam(1:sum(Pl+Ql)) = patternsearch(@(param)cs_dpcLikeRec([param x0(sum(Pl+Ql)+1:end)],eM,data,Pl,Ql,Pc,Qc),...
%         x0(1:sum(Pl+Ql)), AA, bb,[],[],lb,ub,[], options);

% GRIDSEARCH -> PATTERNSEARCH
    grid_x0 = ab_gridsearch( .01, @(param)cs_dpcLikeRec([param x0(sum(Pl+Ql)+1:end)],eM,data,Pl,Ql,Pc,Qc) );
    eparam(1:sum(Pl+Ql)) = patternsearch(@(param)cs_dpcLikeRec([param x0(sum(Pl+Ql)+1:end)],eM,data,Pl,Ql,Pc,Qc),...
        grid_x0, AA, bb,[],[],lb,ub,[], options);
    
    [~,~,~,~,~,~,~,diagLRL,~] = cs_dpcLikeRec( eparam, eM, data, Pl, Ql, Pc, Qc);
    
    eparam(sum(Pl+Ql)+1:end) = cs_dpc_comp_est( diagLRL, Pc, Qc, x0(sum(Pl+Ql)+1:end) );

    like_diff = cs_dpcLikeRec(x0, eM, data, Pl, Ql, Pc, Qc) - cs_dpcLikeRec(eparam, eM, data, Pl, Ql, Pc, Qc);
    disp('===================================================')
    disp(['Quasi-Likelihood Difference in Step ',num2str(step_),': ',num2str(like_diff)])
    disp('===================================================')
    step_ = step_ + 1;
    if like_diff < 0
        eparam = x0; % Reverse last optimization if it worsened likeval
    else % update x0
        x0 = eparam;
    end
end

% options = optimoptions('fmincon','Display','iter');
[nLogL, logLcontr, e_expData, eS, eQ, eL, e_d, e_diagLRL, fcstData ] =...
    cs_dpcLikeRec( eparam, eM, data, Pl, Ql, Pc, Qc);
% [enu, nLogL] = fmincon(@(nu)cawlike_givn2ndpart_( nu, data, LogLcontr_2ndpart ), 2*n,[],[],[],[],n-1,[],[], options);
% logL = -nLogL;
end
