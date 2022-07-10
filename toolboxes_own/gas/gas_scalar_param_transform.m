function [param_out] = gas_scalar_param_transform(dist,param,p)
%GAS_PARAM_TRANSFORM Transforms param vector into individual parameters.
%   Detailed explanation goes here

p_ = p*(p+1)/2;

param_out.Sig.intrcpt = ivechchol(param(1:p_));
param_out.Sig.score1 = param(p_+1);
param_out.Sig.score2 = param(p_+2);
param_out.Sig.garchD = param(p_+3);
param_out.Sig.garchW = param(p_+4);
param_out.Sig.garchM = param(p_+5);
theta = param(p_+6:end);

if strcmp(dist,'Wish')
    param_out.n.intrcpt = theta(1);
    param_out.n.score = theta(2);
    param_out.n.garch = theta(3);
elseif strcmp(dist,'iWish')
    param_out.nu.intrcpt = theta(1);
    param_out.nu.score = theta(2);
    param_out.nu.garch = theta(3);
elseif strcmp(dist,'tWish') || strcmp(dist,'itWish') || strcmp(dist,'F')
    param_out.n.intrcpt = theta(1);
    param_out.n.score = theta(2);
    param_out.n.garch = theta(3);
    param_out.nu.intrcpt = theta(4);
    param_out.nu.score = theta(5);
    param_out.nu.garch = theta(6);
elseif strcmp(dist,'Riesz')
    param_out.n.intrcpt = theta(1:p);
    param_out.n.score = theta(p+1);
    param_out.n.garch = theta(p+2);
elseif strcmp(dist,'iRiesz2')
    param_out.nu.intrcpt = theta(1:p);
    param_out.nu.score = theta(p+1);
    param_out.nu.garch = theta(p+2);  
elseif strcmp(dist,'tRiesz')
    param_out.n.intrcpt = theta(1:p);
    param_out.n.score = theta(p+1);
    param_out.n.garch = theta(p+2); 
    param_out.nu.intrcpt = theta(p+3);
    param_out.nu.score = theta(p+4);
    param_out.nu.garch = theta(p+5);
elseif strcmp(dist,'itRiesz2')
    param_out.n.intrcpt = theta(1);
    param_out.n.score = theta(2);
    param_out.n.garch = theta(3); 
    param_out.nu.intrcpt = theta(4:p+3);
    param_out.nu.score = theta(p+4);
    param_out.nu.garch = theta(p+5);    
elseif strcmp(dist,'FRiesz') || strcmp(dist,'iFRiesz2')
    param_out.n.intrcpt = theta(1:p);
    param_out.n.score = theta(p+1);
    param_out.n.garch = theta(p+2);
    param_out.nu.intrcpt = theta(p+3:p+2+p);
    param_out.nu.score = theta(p+2+p+1);
    param_out.nu.garch = theta(p+2+p+2);
end

end

