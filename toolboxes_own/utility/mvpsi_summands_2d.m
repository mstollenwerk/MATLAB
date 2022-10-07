function y = mvpsi_summands_2d(x,varargin)
%MVPSI Calculates my multivariate version of the polygamma function.
%   See Stollenwerk (2022)

p = size(x,2);
d = ((1:p)-1)/2;

if nargin > 2
    error('Too many input arguments')
elseif nargin == 2
    k = varargin{1};
    y = psi(k, x-d);
else
    y = psi(x-d);
end

end

