function [corrmat] = cov2corr(covmat)
%COV2COR returns correlation matrices given the input covariance matrices.
%
% USAGE:
%  [corrmat] = cov2corr(covmat)
%
% INPUTS:
%   COVMAT 
%
% OUTPUTS:
%   CORRMAT 
%
% COMMENTS:
%
% REFERENCES:
%
% DEPENDENCIES:
%
%  See also 
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 06.05.2019

[k,~,T] = size(covmat);

corrmat = NaN(k,k,T);

for ii = 1:T
    covmat_isnan = isnan(covmat(:,:,ii));
    if ~any(covmat_isnan(:))        
        sdevs = sqrt(diag(covmat(:,:,ii)));
        corrmat(:,:,ii) = covmat(:,:,ii)./(sdevs*sdevs');
    elseif all(covmat_isnan(:))
        corrmat(:,:,ii) = NaN(k);
    else
        error(['covmat ', num2str(ii), ' contains NaN values and is not all NaN'])
    end
end

end

