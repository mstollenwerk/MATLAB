function [param_out] = gas_scalar_har_mean_param_transform_opschoor_scaling(dist,param,p)
%GAS_PARAM_TRANSFORM Transforms param vector into individual parameters.
%   Detailed explanation goes here

p_ = p*(p+1)/2;

param_out.Sig.intrcpt = ivechchol(param(1:p_));
param_out.Sig.score1 = param(p_+1);
param_out.Sig.garchD = param(p_+2);
param_out.Sig.garchW = param(p_+3);
param_out.Sig.garchM = param(p_+4);
theta = param(p_+5:end);

if strcmp(dist,'Wish')
    param_out.n = theta(1);
elseif strcmp(dist,'iWish')
    param_out.nu = theta(1);
elseif strcmp(dist,'tWish') || strcmp(dist,'itWish') || strcmp(dist,'F')
    param_out.n = theta(1);
    param_out.nu = theta(2);
elseif strcmp(dist,'Riesz')
    param_out.n = theta(1:p);
elseif strcmp(dist,'iRiesz2')
    param_out.nu = theta(1:p);
elseif strcmp(dist,'tRiesz')
    param_out.n = theta(1:p);
    param_out.nu = theta(p+1);
elseif strcmp(dist,'itRiesz2')
    param_out.n = theta(1);
    param_out.nu = theta(2:p+1);  
elseif strcmp(dist,'FRiesz') || strcmp(dist,'iFRiesz2')
    param_out.n = theta(1:p);
    param_out.nu = theta(p+1:p+p);
end

end

