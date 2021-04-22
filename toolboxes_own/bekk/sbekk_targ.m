function [ eparam, logL, eS, fcstS ] = sbekk_targ( covMdata, interceptmatrix, P, Q, x0 )
%SBEKK_TARG estimates a restricted scalar BEKK model.
%   SBEKK_TARG  restricts recursion to be covariance 
%   stationary and uses the covariance matrix targeting estimation
%   procedure.
%
%   [eparas]=SBEKK_TARG(covMdata,interceptmatrix,P,Q,x0)  returns only the 
%   coefficient parameters. Inputs are data of covariance matrix
%   estimators, COVMDATA, the (targeted) interceptmatrix, lagged covMdata P,
%   lagged conditional covariance matrices Q, and the starting point
%   for maximum likelihood optimization, X0. X0 can only contain 
%   coefficient starting points.
%
%   [eparas,logL,eS]=SBEKK_TARG  also returns value of the
%   log-likelihood, LOGL, and the recursion of estimated conditional 
%   expected (realized) covariance matrices.

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 29.12.2016

[N,~,T] = size(covMdata);
%% x0 - Starting Values if not provided
if isempty(x0)
    if P==1 && Q==1
        x0=[.03 .96];
    else
        x0=.1*ones(1,P+Q);
    end
end
%% Step 1 - Targeting
if isempty(interceptmatrix)
    interceptmatrix = mean(covMdata,3);
end    
cellcovMdata = mat2cell(covMdata(:,:),N,ones(1,T)*N); % Transforming R to cellformat for recursions (cawqlike and sbekkRec require cell input for quicker evaluation)
%% Step 2 - Maximum Likelihood Estimation
AA = [-eye(2); % Constraints: a,b >= 0 
      1 1];    %              a+b <= 1
bb = [0 0 1];
options = optimoptions('fminunc','Display','off','Algorithm','quasi-newton'); % ,'UseParallel','always'
warning('off','MATLAB:nearlySingularMatrix') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off','MATLAB:illConditionedMatrix') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [eparas,logL] = fmincon(@(params)cawqlike(sbekkRec([vech(chol(interceptmatrix.*(1-sum(params)),'lower'))', params],covMdata,1,1),cellcovMdata ), x0, AA, bb,[],[],[],[],[], options);
if Q==0
    [eparam,logL] = fminunc(@(param_P)cawqlike(sbekkRec_s([vech(chol(interceptmatrix,'lower'))', param_P, 0],covMdata,P,0),cellcovMdata ), x0, options); % Without parameter restrictions, interceptmatrix.*(1-sum(params)) could become neg def., This is why I created sbekkRec_s, to use fminunc.
else
    [eparam,logL] = fminunc(@(param)cawqlike(sbekkRec_s([vech(chol(interceptmatrix,'lower'))', param],covMdata,P,Q),cellcovMdata ), x0, options); % Without parameter restrictions, interceptmatrix.*(1-sum(params)) could become neg def., This is why I created sbekkRec_s, to use fminunc.
end
warning('on','MATLAB:illConditionedMatrix') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('on','MATLAB:nearlySingularMatrix') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logL = -logL;
%% Output
if Q==0
    [eS, fcstS] = sbekkRec_s([vech(chol(interceptmatrix,'lower'))', eparam, 0],covMdata,P,0);
else
    [eS, fcstS] = sbekkRec_s([vech(chol(interceptmatrix,'lower'))', eparam],covMdata,P,Q);
end
eS = reshape(cell2mat(eS),N,N,T);
end