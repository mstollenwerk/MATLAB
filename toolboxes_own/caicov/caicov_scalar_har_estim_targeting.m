function [eparam, tstats, logL, fit, fcst, optimoutput] = ...
	caicov_scalar_har_estim_targeting( X, dist, x0, varargin )
%CAICOV_LIKEREC estimates conditional autoregressive integrated covariance 
%   model using intercept targeting.
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 16.02.2021

%% Input Checking
% Will be added later
[k,~,T] = size(X);
k_ = k*(k+1)/2;
if ~isempty(varargin) && strcmp(varargin{1},'OptimizeOrdering')
    optimize_ordering = 1;
    if numel(varargin)~=1
        varargin = varargin(2:end);
    end
else
    optimize_ordering = 0;    
end
%% Optimization
% Objective Function-------------------------------------------------------
% obj_fun = @(param) caicov_scalar_har_likeRec( ...
%     [vechchol( mean(X,3)*(1-sum(param(1:3))) ); param], ...
%     X, ...
%     dist ...
% );
obj_fun = @(param) obj_fun_wrapper(param, X, dist);
% x0-----------------------------------------------------------------------
if isempty(x0)
    x0 = [.1, .1, .1];
    if strcmp( dist, 'Wish' )
		x0 = [ x0, 2*k ]';
    elseif strcmp( dist, 'iWish' )
		x0 = [ x0, 2*k ]';  
    elseif strcmp( dist, 'tWish' )
		x0 = [ x0, 2*k, 2*k ]';
    elseif strcmp( dist, 'itWish' )
		x0 = [ x0, 2*k, 2*k ]';        
    elseif strcmp( dist, 'F' )
		x0 = [ x0, 2*k, 2*k ]';        
    elseif strcmp( dist, 'Riesz' )
		x0 = [ x0, ones(1,k).*2*k ]';
    elseif strcmp( dist, 'iRiesz' )
		x0 = [ x0, ones(1,k).*2*k ]';         
    elseif strcmp( dist, 'tRiesz' )
		x0 = [ x0, ones(1,k).*2*k, 2*k ]'; 
    elseif strcmp( dist, 'itRiesz' )
		x0 = [ x0, ones(1,k).*2*k, 2*k ]';    
    elseif strcmp( dist, 'FRiesz' )
		x0 = [ x0, ones(1,k).*2*k, ones(1,k).*2*k ]';        
    end
end
% Restrictions-------------------------------------------------------------
if strcmp( dist, 'Wish' )
    lb = [-inf(3,1); k-1];
elseif strcmp( dist, 'iWish' )
    lb = [-inf(3,1); k+1];
elseif strcmp( dist, 'F' )
    lb = [-inf(3,1); k+1; k-1];
elseif strcmp( dist, 'tWish' )
    lb = [-inf(3,1); k-1; 2];
elseif strcmp( dist, 'itWish' )
    lb = [-inf(3,1); k+1; 0];    
elseif strcmp( dist, 'Riesz' )
    lb = [-inf(3,1)', 0:k-1]';  
elseif strcmp( dist, 'iRiesz' )
    lb = [-inf(3,1)', 2:k+1]';  
elseif strcmp( dist, 'tRiesz' )
    lb = [-inf(3,1)', 0:k-1, 2]'; 
elseif strcmp( dist, 'itRiesz' )
    lb = [-inf(3,1)', 0:k-1, 2]';
elseif strcmp( dist, 'FRiesz' )
    lb = [-inf(3,1)', (0:k-1), (k-(1:k)+2)]';    
end
A = [ ones(1,3) zeros(1, length(x0)-3) ];   % Stationarity
b = 1;                                      % Stationarity
% Optimization-------------------------------------------------------------
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp(strcat('Starting caicov scalar har ',dist,' estimation, targeting.'))
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
[eparam,optimoutput] = ...
	my_fmincon_robust(...
		obj_fun, ...
		x0, ...
        lb,[], ...
        [3, numel(x0)]', ...
        1e-2,...
		varargin{:} ...
 	);
% Implementing an algo for permuting order of assets akin to Blasques et
% al.
if optimize_ordering
    
    perm_improvement = inf;
    while perm_improvement > 0
        % Try Loglike values for 100 random asset permutations, find Min.  
        rng(1)
        for ii = 1:100
            perm_(ii,:) = randperm(k);
        end
        perm_ = unique(perm_, 'rows', 'stable');
        for ii = 1:size(perm_,1)
            nLogL_permuted_assets(ii) = ...
                obj_fun_wrapper( eparam, X(perm_(ii,:),perm_(ii,:),:), dist);
        end
        [perm_min_nLogL,min_ii] = min(nLogL_permuted_assets);
        % If min is smaller than current nLogL, permute assets and let
        % solver run again, update current nLogL, save permutation.
        perm_improvement = -perm_min_nLogL + optimoutput{end}.history.fval(end);
        if perm_improvement > 0
            perm_improvement        
            obj_fun = @(param) obj_fun_wrapper( param, X(perm_(min_ii,:),perm_(min_ii,:),:), dist);
            [eparam,optimoutput] = ...
                my_fmincon_robust(...
                    obj_fun, ...
                    eparam, ...
                    lb,[], ...
                    [3, numel(x0)]', ...
                    1e-2,...
                    varargin{:} ...
                );
            optimoutput{end}.perm_ = perm_(min_ii,:);            
        end
        clear perm_ nLogL_permuted_assets
    end
    
end
if ~isfield(optimoutput{end},'perm_')
    optimoutput{end}.perm_ = 1:k;
end
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp(strcat('Finished caicov scalar har ',dist,' estimation, targeting.'))
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
% [eparam,optimoutput] = ...
% 	my_fmincon(...
% 		obj_fun, ...
% 		x0, ...
%         A,b,[],[],lb,[],[], ...
% 		varargin{:} ...
% 	);
% while strcmp(optimoutput.message(50:111),'fmincon stopped because the gradient calculation is undefined.')
%     x0_new = eparam.all;
%     x0_new(1:3) = x0_new(1:3)*.9;
%     [eparam,optimoutput] = ...
%         my_fmincon(...
%             obj_fun, ...
%             x0_new, ...
%             A,b,[],[],lb,[],[], ...
%             varargin{:} ...
%         );
% end
% Tstats-------------------------------------------------------------------
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv( obj_fun, eparam );
tstats = eparam./sqrt(diag(VCV));

% Output Creation----------------------------------------------------------
[ nLogL, logLcontr, Sigma_, eparam ] = obj_fun( eparam );

aic = 2*nLogL + 2*(numel(x0)+k_);
bic = 2*nLogL + log(T)*(numel(x0)+k_); % see Yu, Li and Ng (2017) 
logL = struct(...
    'logL', -nLogL,...
    'aic', aic,...
    'bic', bic,...
    'logLcontr', logLcontr...
);

iperm_ = iperm(optimoutput{end}.perm_);
Sigma_ = Sigma_(iperm_,iperm_,:);
fit = struct('Sigma_', Sigma_(:,:,1:T));
fcst = struct('Sigma_', Sigma_(:,:,T+1:end));

end

function [ nLogL, logLcontr, Sigma_, param ] = obj_fun_wrapper(param, X, dist) 

    if sum(param(1:3)) >= 1
        nLogL = inf;
        return
    end
    
    vechcholSig = vechchol( mean(X,3)*(1-sum(param(1:3))) );
    
    [ nLogL, logLcontr, Sigma_, param ] = caicov_scalar_har_likeRec( ...
        [vechcholSig; param], ...
        X, ...
        dist ...
    );

end