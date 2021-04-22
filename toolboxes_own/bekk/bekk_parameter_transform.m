function [C,A,B] = bekk_parameter_transform(param,p,q,type,k)
%BEKK_PARAMETER_TRANSFORM Parameter transformation for BEKK(p,q) multivariate volatility model.
%
% USAGE:
%  [C,A,B] = bekk_parameter_transform(PARAM,P,Q,K,TYPE)
%
% INPUTS:
%   PARAMETERS - Vector of parameters governing the dynamics.  See BEKKREC
%   P          - Positive, scalar integer representing the number of symmetric innovations
%   Q          - Non-negative, scalar integer representing the number of conditional covariance lags
%   K          - Number of assets
%   TYPE       - String indicating type, 
%                  Scalar (Default, if input is [])
%                  Diagonal
%                  Full
%
% OUTPUTS:
%   C - K by K covariance model intercept
%   A - K by K by P matrix of symmetric innovation parameters
%   B - K by K by Q matrix of smoothing parameters
%
% COMMENTS:
%
% See also BEKKREC

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 15.02.2017

C = ivech(param(1:k*(k+1)/2),'lower');
C = C*C';
if isempty(type) || strcmpi(type,'scalar') % Parameter transformations acc. to different BEKK types.
    if  numel(param)~=k*(k+1)/2+p+q
        error('Number of elements in parameter input does not match BEKK specification')
    end
    A=NaN(k,k,p);
    for i=1:p
        A(:,:,i)=eye(k)*sqrt(param(k*(k+1)/2+i)); % scalar BEKK has scalar in front of the recursed matrices, so the quadratic matrix form has to have sqrt().
    end
    B=NaN(k,k,q);
    for i=1:q
        B(:,:,i)=eye(k)*sqrt(param(k*(k+1)/2+p+i));
    end
elseif strcmpi(type,'diagonal')
    if  numel(param)~=k*(k+1)/2+(p+q)*k
        error('Number of elements in parameter input does not match BEKK specification')
    end
    A=NaN(k,k,p);
    for i=1:p
        A(:,:,i)=diag(param(k*(k+1)/2+1+k*(i-1):k*(k+1)/2+k*i));
    end
    B=NaN(k,k,q);
    for i=1:q
        B(:,:,i)=diag(param(k*(k+1)/2+1+k*p+k*(i-1):k*(k+1)/2+k*p+k*i));
    end
elseif strcmpi(type,'full')
    if  numel(param)~=k*(k+1)/2+(p+q)*k^2
        error('Number of elements in parameter input does not match BEKK specification')
    end
    A = reshape(param(k*(k+1)/2 + 1 : k*(k+1)/2+p*k^2),k,k,p);
    B = reshape(param(k*(k+1)/2 + p*k^2+1 : k*(k+1)/2 + (p+q)*k^2),k,k,q);
else error('BEKK type has to either be "[]", "scalar", "diagonal" or "full".')
end

end