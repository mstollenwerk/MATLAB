function [Omega_, Sigma_, n, nu] = matvParamTrans( dist, param, p )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 24.08.2022
p_ = p*(p+1)/2;

Omega_ = [];
Sigma_ = [];
n = [];
nu = [];

if strcmp( dist, 'Wish' )
elseif strcmp( dist, 'iWish' ) 
elseif strcmp( dist, 'tWish' ) 
elseif strcmp( dist, 'itWish' )
elseif strcmp( dist, 'F' ) 
elseif strcmp( dist, 'Riesz' )
elseif strcmp( dist, 'iRiesz2' )
elseif strcmp( dist, 'tRiesz' )
elseif strcmp( dist, 'itRiesz2' )
elseif strcmp( dist, 'FRiesz' )
    for ii = 1:size(param,2)
        Omega_(:,:,ii) = ivechchol(param(1:p_, ii));
        n(:,ii) = param(p_ + 1 : p_ + p, ii);
        nu(:,ii) = param(p_ + p + 1 : p_ + p + p, ii);  
    end
elseif strcmp( dist, 'iFRiesz2' )
    for ii = 1:size(param,2)
        Omega_(:,:,ii) = ivechchol(param(1:p_, ii));
        n(:,ii) = param(p_ + 1 : p_ + p, ii);
        nu(:,ii) = param(p_ + p + 1 : p_ + p + p, ii);  
    end    
elseif strcmp( dist, 'sWish' ) 
elseif strcmp( dist, 'siWish' ) 
elseif strcmp( dist, 'stWish' ) 
elseif strcmp( dist, 'sitWish' )
elseif strcmp( dist, 'sF' ) 
elseif strcmp( dist, 'sRiesz' )
elseif strcmp( dist, 'siRiesz2' )
elseif strcmp( dist, 'stRiesz' )
elseif strcmp( dist, 'sitRiesz2' )
elseif strcmp( dist, 'sFRiesz' )
    for ii = 1:size(param,2)
        Sigma_(:,:,ii) = ivechchol(param(1:p_, ii));
        n(:,ii) = param(p_ + 1 : p_ + p, ii);
        nu(:,ii) = param(p_ + p + 1 : p_ + p + p, ii);  
    end    
elseif strcmp( dist, 'siFRiesz2' )
    for ii = 1:size(param,2)
        Sigma_(:,:,ii) = ivechchol(param(1:p_, ii));
        n(:,ii) = param(p_ + 1 : p_ + p, ii);
        nu(:,ii) = param(p_ + p + 1 : p_ + p + p, ii);  
    end        
end

end
