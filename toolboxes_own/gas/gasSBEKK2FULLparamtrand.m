function [param_out] = gasSBEKK2FULLparamtrand(param,dist,p)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
p_ = p*(p+1)/2;

param_out = [param(1:p_); 
             ones(p,1)*sqrt(param(p_+1)); 
             ones(p,1)*sqrt(param(p_+2)); 
             ones(p,1)*sqrt(param(p_+3));
             zeros(2*p,1)];

zeros_ = zeros(4,1);
if strcmp(dist,'Wish') || strcmp(dist,'iWish')
    df = [param(p_+4); zeros_];
elseif strcmp(dist,'tWish') || strcmp(dist,'itWish') || strcmp(dist,'F')
    df = [param(p_+4); zeros_; param(p_+5); zeros_];
elseif strcmp(dist,'Riesz') || strcmp(dist,'iRiesz')
    df = [param(p_+4:p_+3+p); zeros_];
elseif strcmp(dist,'tRiesz')
    df = [param(p_+4:p_+3+p); zeros_; param(p_+4+p); zeros_];
elseif strcmp(dist,'itRiesz')
    df = [param(p_+4); zeros_; param(p_+5:p+4+p); zeros_];
elseif strcmp(dist,'FRiesz') || strcmp(dist,'iFRiesz')
    df = [param(p_+4:p_+3+p); zeros_; param(p_+4+p:p_+3+p+p); zeros_];
end
param_out = [param_out;df];

end

