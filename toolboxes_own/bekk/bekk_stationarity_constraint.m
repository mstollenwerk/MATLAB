function [c,ceq] = bekk_stationarity_constraint(param,p,q,type,k)
% Non-linear constraint for estimation of BEKK(p,q) multivariate volatility models
%   These are stationarity conditions.
%
% USAGE:
%  [C,CEQ] = bekk_constraint(A,B,P,Q,TYPE,K) 
%
% INPUTS:
%   See bekkRec
%
% OUTPUTS:
%   C   - Vector of inequality constraints
%   CEQ - Empty
%
% COMMENTS:
%
%  EXAMPLES:
%
% See also BEKKREC
%
% REFERENCES:
%      [1] Lütkepol - New Introduction to Multiple Time Series Analysis,
%          p.565.
%
% based on:
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 3/27/2012

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 15.02.2017

[~,A,~,B] = bekk_parameter_transform(param,p,0,q,k,type);

ceq = [];

switch type
    case 'scalar'
        c = sum(A(1,1,:).^2,3) + sum(B(1,1,:).^2,3) - 1;
    case []
        c = sum(A(1,1,:).^2,3) + sum(B(1,1,:).^2,3) - 1;        
    case 'diagonal'
        c = diag(sum(A.^2,3) + sum(B.^2,3) - 1);
    case 'full'
        L = ELmatrix(k);
        D = Dmatrix(k);
        m = zeros(k*(k+1)/2,k*(k+1)/2);
        for i=1:p
            m = m + L*kron(A(:,:,i),A(:,:,i))*D;
        end
        for i=1:q
            m = m + L*kron(B(:,:,i),B(:,:,i))*D;
        end
        c = abs(eig(m)) - 1 + eps; % Biggest eigenvalue has to be smaller than 1 in absolute value.
end
