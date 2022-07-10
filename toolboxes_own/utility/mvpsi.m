function y = mvpsi(n,varargin)
%MVPSI Calculates my multivariate version of the polygamma function.
%   See Stollenwerk (2022)

p = length(n);
x = ((1:p)-1)/2;
x = reshape(x,size(n));

if nargin > 2
    error('Too many input arguments')
elseif nargin == 2
    k = varargin{1};
    y = psi(k, n-x);
else
    y = psi(n-x);
end

end

