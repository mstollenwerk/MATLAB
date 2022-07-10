function [ eparam, tstats, logL, optimoutput ] = matvRCest(R,dist,x0,varargin)
%MATVRCEST Static estimation of all distributions in Stollenwerk (2022).
%   Detailed explanation goes here
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 08.06.2022

narginchk(2,inf);
[p,~,T] = size(R);
p_ = p*(p+1)/2;

%% Optimization Starting Points and Lower Bounds for Degrees of Freedom
% The lower bounds below are for existence of the distributions. For
% existence of the expected value matrix Sigma_, restrictions on the dfs
% are stricter. The inverse Wishart, for example, exists for nu > p-1, but
% its mean only exists for nu > p+1.
lb_df_barlett = p-1;
lb_df_barlettL = 0:p-1;
lb_df_barlettU = flip(0:p-1);
lb_df_chi = 0;
x0_df_barlett = 2*p;
x0_df_barlettL = 2*p*ones(1,p);
x0_df_barlettU = 2*p*ones(1,p);
x0_df_chi = 5;

if strcmp( dist, 'Wish' )
    x0_df = x0_df_barlett;
    lb_df = lb_df_barlett;
elseif strcmp( dist, 'iWish' )
    x0_df = x0_df_barlett;  
    lb_df = lb_df_barlett;
elseif strcmp( dist, 'tWish' )
    x0_df = [ x0_df_barlett, x0_df_chi ];
    lb_df = [ lb_df_barlett, lb_df_chi ];
elseif strcmp( dist, 'itWish' )
    x0_df = [ x0_df_chi, x0_df_barlett ];
    lb_df = [ lb_df_chi, lb_df_barlett ];
elseif strcmp( dist, 'F' )
    x0_df = [ x0_df_barlett, x0_df_barlett ];   
    lb_df = [ lb_df_barlett, lb_df_barlett ]; 
elseif strcmp( dist, 'Riesz' )
    x0_df = x0_df_barlettL;
    lb_df = lb_df_barlettL;
elseif strcmp( dist, 'iRiesz2' )
    x0_df = x0_df_barlettU;
    lb_df = lb_df_barlettU;
elseif strcmp( dist, 'tRiesz' )
    x0_df = [ x0_df_barlettL, x0_df_chi ]; 
    lb_df = [ lb_df_barlettL, lb_df_chi ];
elseif strcmp( dist, 'itRiesz2' )
    x0_df = [ x0_df_chi, x0_df_barlettU ];
    lb_df = [ lb_df_chi, lb_df_barlettU];
elseif strcmp( dist, 'FRiesz' )
    x0_df = [ x0_df_barlettL, x0_df_barlettU ];
    lb_df = [ lb_df_barlettL, lb_df_barlettU];
elseif strcmp( dist, 'iFRiesz2' )
    x0_df = [ x0_df_barlettL, x0_df_barlettU ];
    lb_df = [ lb_df_barlettL, lb_df_barlettU ];
end
if ~isempty(x0)
    x0_df = x0;
end

if ~isempty(varargin) && strcmp(varargin{1},'EstimationStrategy:MethodOfMoments')
    %% Optimization

    meanR = mean(R,3);
    
    if contains(dist,'Riesz')
        
        obj_fun = @(dfs,perm_) matvsLogLike_wrapper( ...
            dist, meanR(perm_,perm_,:), dfs, R(perm_,perm_,:) ...
         );        
        [eparam,optimoutput] = ...
            fmincon_Rieszperm(...
                p, ...
                obj_fun, ...
                x0_df', ...
                [],[],[],[],lb_df',[],[],varargin{2:end} ...
            );             
        obj_fun = @(dfs) obj_fun(dfs,optimoutput.perm_); 
        
    else
        
        obj_fun = @(dfs) matvsLogLike_wrapper( dist, meanR, dfs, R);        
        [eparam,optimoutput] = ...
            my_fmincon(...
                obj_fun, ...
                x0_df', ...
                [],[],[],[],lb_df',[],[],varargin{2:end} ...
            );
        
    end
    %% tstats
    %[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
    [VCV,scores,gross_scores] = vcv(obj_fun, eparam);
    tstats = eparam./sqrt(diag(VCV));
    tstats = struct(...
        'Sigma_', NaN(p_,1), ...
        'dfs', tstats, ...               
        'all', [NaN(p_,1); tstats] ...
    );
    %% nLogL, logLcontr and eparam
    [nLogL, logLcontr, ~, ~, eparam ] = obj_fun( eparam );

    aic = 2*nLogL + 2*(numel(x0) + p_);
    bic = 2*nLogL + log(T)*(numel(x0) + p_);
    logL = struct(...
        'logL', -nLogL, ...
        'logLcontr', logLcontr, ...
        'bic', bic, ...
        'aic', aic ...
    );
else
    error('Only EstimationStrategy:MethodOfMoments supported so far.')
end

end

function [ nLogL, logLcontr, score, hessian, param ] = ...
    matvsLogLike_wrapper(dist, Sigma_, dfs, R)

    try
        [ nLogL, logLcontr, score, hessian, param ] = ...
            matvsLogLike(dist, Sigma_, dfs, R);
    catch
        nLogL = NaN;
        logLcontr = NaN;
        score = NaN;
        hessian = NaN;
        param = NaN;
        return
    end

end

