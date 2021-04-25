function [eparam, tstats, logL, fit, fcst, optimoutput] = ...
	gas_scalar_BEKK_estim( X, p, q, dist, x0, varargin )
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 25.03.2021

%% Input Checking
% Will be added later
[k,~,T] = size(X);
k_ = k*(k+1)/2;
meanX = mean(X,3);
if ~isempty(varargin) && strcmp(varargin{1},'OptimizeOrdering')
    optimize_ordering = 1;
    if numel(varargin)~=1
        varargin = varargin(2:end);
    end
else
    optimize_ordering = 0;    
end
%% Optimization
obj_fun = @(param) gas_scalar_BEKK_likeRec(param, p, q, X, dist);
% x0-----------------------------------------------------------------------
if isempty(x0)
    x0 = [vechchol(meanX.*.9)' ones(1,p)*.01/p ones(1,q)*.75/q];
    if strcmp( dist, 'Wish' )
		x0 = [ x0, 2*k ]';
    elseif strcmp( dist, 'iWish' )
		x0 = [ x0, 2*k ]';  
    elseif strcmp( dist, 'tWish' )
		x0 = [ x0, 2*k, 5 ]';
    elseif strcmp( dist, 'itWish' )
		x0 = [ x0, 2*k, 5 ]';
    elseif strcmp( dist, 'F' )
		x0 = [ x0, 2*k, 2*k ]';        
    elseif strcmp( dist, 'Riesz' )
		x0 = [ x0, ones(1,k).*(2*k) ]';
    elseif strcmp( dist, 'iRiesz' )
		x0 = [ x0, ones(1,k).*(2*k) ]';         
    elseif strcmp( dist, 'tRiesz' )
		x0 = [ x0, ones(1,k).*(2*k), 5 ]'; 
    elseif strcmp( dist, 'itRiesz' )
		x0 = [ x0, ones(1,k).*(2*k), 5 ]';    
    elseif strcmp( dist, 'FRiesz' )
		x0 = [ x0, ones(1,k).*(2*k), ones(1,k).*(2*k) ]';        
    end
end
% Restrictions-------------------------------------------------------------
if strcmp( dist, 'Wish' )
    lb = [-inf(k_+p+q,1); k-1];
elseif strcmp( dist, 'iWish' )
    lb = [-inf(k_+p+q,1); k+1];
elseif strcmp( dist, 'F' )
    lb = [-inf(k_+p+q,1); k+1; k-1];
elseif strcmp( dist, 'tWish' )
    lb = [-inf(k_+p+q,1); k-1; 2];
elseif strcmp( dist, 'itWish' )
    lb = [-inf(k_+p+q,1); k+1; 0];    
elseif strcmp( dist, 'Riesz' )
    lb = [-inf(k_+p+q,1)', 0:k-1]';  
elseif strcmp( dist, 'iRiesz' )
    lb = [-inf(k_+p+q,1)', 2:k+1]';  
elseif strcmp( dist, 'tRiesz' )
    lb = [-inf(k_+p+q,1)', 0:k-1, 2]'; 
elseif strcmp( dist, 'itRiesz' )
    lb = [-inf(k_+p+q,1)', 0:k-1, 2]';
elseif strcmp( dist, 'FRiesz' )
    lb = [-inf(k_+p+q,1)', (0:k-1), (k-(1:k)+2)]';    
end
A = [ zeros(1,k_+p) ones(1,q) zeros(1, length(x0)-k_-p-q) ];   % Stationarity
b = 1;                                                         % Stationarity
% Optimization-------------------------------------------------------------
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp(strcat("Starting gas scalar BEKK(",num2str(p),",",num2str(q),")-",dist," estimation."))
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
[eparam,optimoutput] = ...
	my_fmincon(...
		obj_fun, ...
		x0, ...
        A,b,[],[],lb,[],[], ...
		varargin{:} ...
	);
% Implementing an algo for permuting order of assets akin to Blasques et
% al.
if optimize_ordering
    
    perm_improvement = inf;
    while perm_improvement > 0
        rng(1)        
        for ii = 1:100
            perm_(ii,:) = randperm(k);
        end
        perm_ = unique(perm_, 'rows', 'stable');
        for ii = 1:size(perm_,1)
            nLogL_permuted_assets(ii) = ...
                gas_scalar_BEKK_likeRec(eparam, p, q, X(perm_(ii,:),perm_(ii,:),:), dist);
        end
        [perm_min_nLogL,min_ii] = min(nLogL_permuted_assets);
        % If min is smaller than current nLogL, permute assets and let
        % solver run again, update current nLogL, save permutation.
        perm_improvement = -perm_min_nLogL + optimoutput{end}.history.fval(end);
        if perm_improvement > 0
            perm_improvement
            obj_fun = @(param) gas_scalar_BEKK_likeRec(param, p, q, X(perm_(min_ii,:),perm_(min_ii,:),:), dist);
            [eparam,optimoutput] = ...
                my_fmincon(...
                    obj_fun, ...
                    eparam, ...
                    A,b,[],[],lb,[],[], ...
                    varargin{:} ...
                );
            optimoutput.perm_ = perm_(min_ii,:);
        end
        clear perm_ nLogL_permuted_assets
    end
    
end
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp(strcat("Finished gas scalar BEKK(",num2str(p),",",num2str(q),")-",dist," estimation."))
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
% Tstats-------------------------------------------------------------------
%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv( obj_fun, eparam );
tstats = eparam./sqrt(diag(VCV));

% Output Creation----------------------------------------------------------
[ nLogL, logLcontr, Sigma_, ScaledScore, eparam ] = obj_fun( eparam );

aic = 2*nLogL + 2*(numel(x0));
bic = 2*nLogL + log(T)*(numel(x0)); % see Yu, Li and Ng (2017) 
logL = struct(...
    'logL', -nLogL,...
    'aic', aic,...
    'bic', bic,...
    'logLcontr', logLcontr...
);

fit = struct( ...
    'Sigma_', Sigma_(:,:,1:T), ...
    'ScaledScore', ScaledScore ...
);
fcst = struct('Sigma_', Sigma_(:,:,T+1:end));

end