function [ eparam, tstats, logL, optimoutput ] = matvest(R,dist,x0_df,varargin)
%MATVEST Simple wrapper around my matrix-variate distribution estimations

[p,~,T] = size(R);
p_ = p*(p+1)/2;
meanR = mean(R,3);

obj_fun = @(dfs, perm) matvsLogLike_errorINF(...
    dist, chol(meanR(perm,perm,:),'lower'), dfs, chol3d(R(perm,perm,:))...
);

if isempty(x0_df)
    x0_df = matvX0df(dist,p)';
end
lb_df = matvLBdf(dist,p);

if contains(dist,'Riesz')
    perm_optim = 1;
else
    perm_optim = 0;
end
if ~isempty(varargin) && strcmp(varargin{1},'noperm')
    perm_optim = 0;
    varargin = varargin(2:end);
end
if perm_optim
    [eparam,optimoutput] = ...
        fmincon_Rieszperm(...
            p, ...
            obj_fun, ...
            x0_df, ...
            [],[],[],[], ...
            lb_df,[],[], ...
            varargin{:} ...
        );
else
    [eparam,optimoutput] = ...
        my_fmincon(...
            @(param) obj_fun(param,1:p), ...
            x0_df, ...
            [],[],[],[], ...
            lb_df,[],[], ...
            varargin{:} ...
        );
    optimoutput.perm_ = 1:p;
end

% Output Creation----------------------------------------------------------
[ nLogL, logLcontr, ~, ~, eparam ] = obj_fun( eparam, optimoutput.perm_ );
eparam.perm_ = optimoutput.perm_;

%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv( obj_fun, eparam.all(p_+1:end), optimoutput.perm_ );
tstats = eparam.all(p_+1:end)./sqrt(diag(VCV));


aic = 2*nLogL + 2*(numel(x0_df)+p_);
bic = 2*nLogL + log(T)*(numel(x0_df)+p_); % see Yu, Li and Ng (2017) 
logL_detR_part = -(p+1)/2*sum(logdet3d(R));
logL = struct(...
    'logL', -nLogL,...
    'logL_detR_part', logL_detR_part, ...
    'aic', aic,...
    'bic', bic,...
    'logLcontr', logLcontr...
);

end

function [ nLogL, logLcontr, varargout ] = ...
    matvsLogLike_errorINF(dist, C, dfs, Cr) 

    [p,~,N] = size(Cr);
    Cz = NaN(p,p,N);
    I = eye(p);
    
    M = matvExpMat(dist,dfs,p);
    sqrtM = sqrt(M);
    COm = C/sqrtM;
    
    try
        for ii = 1:N
            Cz(:,:,ii) = COm\Cr(:,:,ii);
            if any(strcmp( dist, {'iWish','iRiesz2','itWish','itRiesz2'} ))
                X(ii,:) = sum(sum(inv(Cz(:,:,ii)).^2));
            elseif strcmp( dist, 'F' )
                X(ii,:) = logdet(I + Cz(:,:,ii)*Cz(:,:,ii)');
            elseif strcmp( dist, 'FRiesz' )
                X(ii,:) = diag(chol(I + Cz(:,:,ii)*Cz(:,:,ii)','lower'));
            elseif strcmp( dist, 'iFRiesz2' )
                C_Z_twisted_t = Cr(:,:,ii)'/COm';
                B = I + C_Z_twisted_t*C_Z_twisted_t';
                X(ii,:) = diag(cholU(B));
%                 X(ii,:) = diag(chol(inv(),'lower'));
%                 C_iZ = COm'/Cr(:,:,ii)';
%                 X(ii,:) = diag(chol(inv(I + C_iZ*C_iZ'),'lower'));
            else
                X = [];
            end        
        end
    catch
        nLogL = inf;
        logLcontr = NaN;
        return
    end
    
    if nargout == 5
        R = NaN(p,p,N);
        for ii = 1:N
            R(:,:,ii) = Cr(:,:,ii)*Cr(:,:,ii)';
        end
        [ nLogL, logLcontr, score, hessian, param ] = ...
            matvsLogLike(dist, C*C', dfs, R);
        varargout{1} = score;
        varargout{2} = hessian;
        varargout{3} = param;
    elseif nargout <= 2
        [nLogL, logLcontr] = matvLogLikeCz( dist, Cz, dfs', X );
    end


end

