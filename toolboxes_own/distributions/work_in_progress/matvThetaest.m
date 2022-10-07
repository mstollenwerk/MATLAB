function [ eparam, logL, optimoutput ] = matvThetaest(Cr,C,dist,x0,varargin)

[p,~,T] = size(Cr);
p_ = p*(p+1)/2;
Cz = NaN(p,p,T);
for ii = 1:T
    Cz(:,:,ii) = C(:,:,ii)\Cr(:,:,ii);
end

if isempty(x0)
    x0 = matvX0df(dist,p);
end

varargin = [varargin, {'Algorithm','trust-region', ...
                       'StepTolerance', 1e-10, ...
                       'FunctionTolerance', 1e-10, ...
                       'SpecifyObjectiveGradient', true}               ];
[eparam,~] = ...
    my_fminunc(...
        @(dfs) obj_fun_wrapper_ini(...
            dist, Cz, dfs ...
        ), ...
        x0, ...
        varargin{:} ...
    );

nLogL = obj_fun_wrapper_ini(dist,Cz,eparam);

tol = 1;
while tol > 1e-10
    
    E = matvExpMat(dist, eparam, p);   
    Cz_Om = NaN(p,p,T);
    for ii = 1:T
        Cz_Om(:,:,ii) = sqrt(E)*Cz(:,:,ii);
    end
    [eparam_,~] = ...
        my_fminunc(...
            @(dfs) obj_fun_wrapper(...
                dist, Cz_Om, eparam ...
            ), ...
            x0, ...
            varargin{:} ...
        );
    nLogL_ = obj_fun_wrapper(dist,Cz_Om,eparam_);
    tol = nLogL - nLogL_;
    if tol > 0 
        eparam = eparam_;
        nLogL = nLogL_;
    end
        
end
E = matvExpMat(dist, eparam, p);   
Cz_Om = NaN(p,p,T);
for ii = 1:T
    Cz_Om(:,:,ii) = sqrt(E)*Cz(:,:,ii);
end

[ ~, ~, ~, logLcontr, g, H ] = obj_fun_wrapper(dist,Cz_Om,eparam_);
logL_detR_part = -(p+1)*sum(log(diag3d(Cr)),2);
logLcontr = logLcontr + logL_detR_part;
nLogL = -sum(logLcontr);

aic = 2*nLogL + 2*(numel(x0)+p_);
bic = 2*nLogL + log(T)*(numel(x0)+p_); % see Yu, Li and Ng (2017) 
logL = struct(...
    'logL', -nLogL,...
    'logLcontr', logLcontr,...
    'logL_detR_part', sum(logL_detR_part), ...
    'aic', aic,...
    'bic', bic...
);

optimoutput.gradient = g;
optimoutput.Hessian = H;

end

function [ nLogL, sumNg, sumN_H, logLcontr, g, H ] = ...
    obj_fun_wrapper(dist, Cz_Om, dfs) 
    
    [p,~,T] = size(Cz_Om);

%     E = matvExpMat(dist, dfs, p);   
%     Cz_Om = NaN(p,p,T);
%     for ii = 1:T
%         Cz_Om(:,:,ii) = sqrt(E)*Cz(:,:,ii);
%     end
    
    X = [];
    I = eye(p);
    if strcmp( dist, 'iWish' ) || strcmp( dist, 'itWish' ) || strcmp( dist, 'iRiesz2' ) || strcmp( dist, 'itRiesz2' )
        tr_iZ_Om = reshape(sum(sum(inv3d(Cz_Om).^2)),[],1);
%         tr_iZ = reshape(sum(sum(inv3d(Cz).^2)),[],1);
        X = tr_iZ_Om;
    elseif strcmp( dist, 'F' )
        logdet_IpZ_Om = NaN(T,1);
        for ii = 1:T
            logdet_IpZ_Om(ii) = logdet(I + Cz_Om(:,:,ii)*Cz_Om(:,:,ii)');
        end
        X = logdet_IpZ_Om;
    elseif strcmp( dist, 'FRiesz' )
        diag3d_chol_IpZ_Om = NaN(T,p);
        for ii = 1:T
            diag3d_chol_IpZ_Om(ii,:) = diag(chol(I + Cz_Om(:,:,ii)*Cz_Om(:,:,ii)','lower'));
        end       
        X = diag3d_chol_IpZ_Om;
    elseif strcmp( dist, 'iFRiesz2' )
        diag3d_chol_invIpInvZ_Om = NaN(T,p);
        iCz = inv3d(Cz_Om);
        for ii = 1:T
            diag3d_chol_invIpInvZ_Om(ii,:) = ...
                diag(chol(inv(I + iCz(:,:,ii)'*iCz(:,:,ii)),'lower'));
        end         
        X = diag3d_chol_invIpInvZ_Om;
    end
    
    try
        [nLogL, logLcontr] = matvLogLikeCz(dist, Cz_Om, dfs, X);
    catch
        nLogL = NaN;
        logLcontr = NaN;
        g = NaN;
        H = NaN;
        sumNg = NaN;
        sumN_H = NaN;
        return
    end
    g = matvLogLikeGradientThetaCz(dist, Cz_Om, dfs, X);
%     g = g + dfs*p/2/(dfs-p-1) - 1/2*tr_iZ;
    H = matvLogLikeHessianThetaCz(dist, Cz_Om, dfs, X);
%     H = H + p/2/(dfs-p-1) + p/2/(dfs-p-1) - dfs*p/2/(dfs-p-1)^2;
    sumNg = -sum(g,1);
    sumN_H = -sum(H,3);   

end

function [ nLogL, sumNg, sumN_H, logLcontr, g, H ] = ...
    obj_fun_wrapper_ini(dist, Cz, dfs) 
    
    [p,~,T] = size(Cz);

    E = matvExpMat(dist, dfs, p);   
    Cz_Om = NaN(p,p,T);
    for ii = 1:T
        Cz_Om(:,:,ii) = sqrt(E)*Cz(:,:,ii);
    end
    
    X = [];
    I = eye(p);
    if strcmp( dist, 'iWish' ) || strcmp( dist, 'itWish' ) || strcmp( dist, 'iRiesz2' ) || strcmp( dist, 'itRiesz2' )
        tr_iZ_Om = reshape(sum(sum(inv3d(Cz_Om).^2)),[],1);
%         tr_iZ = reshape(sum(sum(inv3d(Cz).^2)),[],1);
        X = tr_iZ_Om;
    elseif strcmp( dist, 'F' )
        logdet_IpZ_Om = NaN(T,1);
        for ii = 1:T
            logdet_IpZ_Om(ii) = logdet(I + Cz_Om(:,:,ii)*Cz_Om(:,:,ii)');
        end
        X = logdet_IpZ_Om;
    elseif strcmp( dist, 'FRiesz' )
        diag3d_chol_IpZ_Om = NaN(T,p);
        for ii = 1:T
            diag3d_chol_IpZ_Om(ii,:) = diag(chol(I + Cz_Om(:,:,ii)*Cz_Om(:,:,ii)','lower'));
        end       
        X = diag3d_chol_IpZ_Om;
    elseif strcmp( dist, 'iFRiesz2' )
        diag3d_chol_invIpInvZ_Om = NaN(T,p);
        iCz = inv3d(Cz_Om);
        for ii = 1:T
            diag3d_chol_invIpInvZ_Om(ii,:) = ...
                diag(chol(inv(I + iCz(:,:,ii)'*iCz(:,:,ii)),'lower'));
        end         
        X = diag3d_chol_invIpInvZ_Om;
    end
    
    try
        [nLogL, logLcontr] = matvLogLikeCz(dist, Cz_Om, dfs, X);
    catch
        nLogL = NaN;
        logLcontr = NaN;
        g = NaN;
        H = NaN;
        sumNg = NaN;
        sumN_H = NaN;
        return
    end
    g = matvLogLikeGradientThetaCz(dist, Cz_Om, dfs, X);
%     g = g + dfs*p/2/(dfs-p-1) - 1/2*tr_iZ;
    H = matvLogLikeHessianThetaCz(dist, Cz_Om, dfs, X);
%     H = H + p/2/(dfs-p-1) + p/2/(dfs-p-1) - dfs*p/2/(dfs-p-1)^2;
    sumNg = -sum(g,1);
    sumN_H = -sum(H,3);   

end

