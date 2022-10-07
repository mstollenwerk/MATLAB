function [ eparam, tstats, logL, optimoutput ] = matvDFest(C,Cr,x0,varargin)
%MATVEST Simple wrapper around my matrix-variate distribution estimations

[p,~,T] = size(Cr);
p_ = p*(p+1)/2;

% obj_fun = @(dfs, perm) matvsLogLike_errorINF(...
%     dist, C, dfs, C_R(perm,perm,:)...
% );
obj_fun = @(dfs) matvsLogLike_errorINF(...
    C, dfs, Cr...
);

% x0_df_barlett = 2*p;
x0_df_barlettL = 2*p*ones(1,p);
% x0_df_barlettU = 2*p*ones(1,p);
% x0_df_chi = 5;

% if strcmp( dist, 'Wish' )
%     x0_df = x0_df_barlett;
% elseif strcmp( dist, 'iWish' )
%     x0_df = x0_df_barlett;
% elseif strcmp( dist, 'tWish' )
%     x0_df = [ x0_df_barlett, x0_df_chi ];
% elseif strcmp( dist, 'itWish' )
%     x0_df = [ x0_df_chi, x0_df_barlett ];
% elseif strcmp( dist, 'F' )
%     x0_df = [ x0_df_barlett, x0_df_barlett ];
% elseif strcmp( dist, 'Riesz' )
    x0_df = x0_df_barlettL;
% elseif strcmp( dist, 'iRiesz2' )
%     x0_df = x0_df_barlettU;
% elseif strcmp( dist, 'tRiesz' )
%     x0_df = [ x0_df_barlettL, x0_df_chi ];
% elseif strcmp( dist, 'itRiesz2' )
%     x0_df = [ x0_df_chi, x0_df_barlettU ];
% elseif strcmp( dist, 'FRiesz' )
%     x0_df = [ x0_df_barlettL, x0_df_barlettU ];
% elseif strcmp( dist, 'iFRiesz2' )
%     x0_df = [ x0_df_barlettL, x0_df_barlettU ];
% end

% if contains(dist,'Riesz')
%     perm_optim = 1;
% else
%     perm_optim = 0;
% end
% if ~isempty(varargin) && strcmp(varargin{1},'noperm')
%     perm_optim = 0;
%     varargin = varargin(2:end);
% end
% if perm_optim
%     [eparam,optimoutput] = ...
%         fmincon_Rieszperm(...
%             p, ...
%             obj_fun, ...
%             x0_df', ...
%             [],[],[],[], ...
%             [],[],[], ...
%             varargin{:} ...
%         );
% else
    [eparam,optimoutput] = ...
        my_fmincon(...
            @(param) obj_fun(param), ...
            x0_df', ...
            [],[],[],[], ...
            [],[],[], ...
            varargin{:} ...
        );
%     optimoutput.perm_ = 1:p;
% end

% Output Creation----------------------------------------------------------
[ nLogL, logLcontr, ~, ~, eparam ] = obj_fun( eparam );
% eparam.perm_ = optimoutput.perm_;

%[VCV,A,B,scores,hess,gross_scores] = robustvcv(fun, eparam, 3);
[VCV,scores,gross_scores] = vcv( obj_fun, eparam.all(p_+1:end) );
tstats = eparam.all(p_+1:end)./sqrt(diag(VCV));


aic = 2*nLogL + 2*(numel(x0)+p_);
bic = 2*nLogL + log(T)*(numel(x0)+p_); % see Yu, Li and Ng (2017) 
% logL_detR_part = -(p+1)/2*sum(logdet3d(C_R));
logL = struct(...
    'logL', -nLogL,...%'logL_detR_part', logL_detR_part, ...
    'aic', aic,...
    'bic', bic,...
    'logLcontr', logLcontr...
);

end

function [ nLogL, logLcontr, varargout ] = ...
    matvsLogLike_errorINF(C, dfs, Cr) 
    
    if nargout == 5
        try
            [ nLogL, logLcontr, score, hessian, param ] = ...
                matvsRieszlikeC(C, dfs, Cr);
            varargout{1} = score;
            varargout{2} = hessian;
            varargout{3} = param;
        catch
            nLogL = inf;
            logLcontr = NaN;
            varargout{1} = NaN;
            varargout{2} = NaN;
            varargout{3} = NaN;
        end
    elseif nargout <= 2
        try
            [ nLogL, logLcontr ] = ...
                matvsRieszlikeC(C, dfs, Cr);
        catch
            nLogL = inf;
            logLcontr = NaN;
        end
    end


end

