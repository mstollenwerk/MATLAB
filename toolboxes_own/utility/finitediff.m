function df = finitediff(fun, x, d, varargin)
%FINITEDIFF estimates a gradient by finite-differencing method.
% (c) Stefan Harmeling, 2012-07-10.
sx = size(x);
nx = numel(x);
df = zeros(sx);
dx = zeros(sx);
for i = 1:nx
dx(i) = d;
df(i) = (fun(x+dx, varargin{:})-fun(x-dx, varargin{:}))/(2*d);
dx(i) = 0;
end