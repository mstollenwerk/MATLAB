function [g, H] = grad_hessian_2sided_ij(f,x,idx_i,idx_j,varargin)
% Computes 2-sided finite difference Hessian
%
% USAGE:
%   H=hessian_2sided(FUNC,X,VARARGIN)
%
% INPUTS:
%   FUNC         - Function name, fval = func(x,varargin)
%   X            - Vector of parameters (N x 1)
%   VARARGIN     - Optional arguments passed to the function
%
% OUTPUTS:
%   H            - Finite differnce, 2-sided hessian
%
% COMMENTS:
%   Modification of hessian() from the jpl toolbox


% Code originally from COMPECON toolbox [www4.ncsu.edu/~pfackler]
% documentation modified to fit the format of the Econometrics Toolbox
% by James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com
%
% Further modified (to do 2-sided numerical derivs, rather than 1) by:
% Modifications Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

n = size(x,1);

if size(x,2)~=1
    error('X must be a column vector.')
end

try
    feval(f,x,varargin{:});
catch FE
    errtxt=['There was an error evaluating the function.  Please check the arguements.  The specific error was:' FE.message];
    error(errtxt);
end


fx = feval(f,x,varargin{:});

% Compute the stepsize (h)
h = eps.^(1/3)*max(abs(x),1e-8);
xh = x+h;
h = xh-x;
ee = sparse(1:n,1:n,h,n,n);

% Compute forward and backward steps
gp = zeros(n,1);
gm = zeros(n,1);
parfor i=union(idx_i,idx_j)
    gp(i) = feval(f,x+ee(:,i),varargin{:});
    gm(i) = feval(f,x-ee(:,i),varargin{:});
end
g = (gp-gm)./(2*h);

hh=h*h';
Hm=NaN*ones(n);
Hp=NaN*ones(n);
% Compute "double" forward and backward steps
for i=idx_i
    parfor j=idx_j
        Hp(i,j) = feval(f,x+ee(:,i)+ee(:,j),varargin{:});
        Hm(i,j) = feval(f,x-ee(:,i)-ee(:,j),varargin{:});
    end
end
%Compute the hessian
H = zeros(n);
for i=idx_i
    for j=idx_j
        H(i,j) = (Hp(i,j)-gp(i)-gp(j)+fx+fx-gm(i)-gm(j)+Hm(i,j))/hh(i,j)/2;
    end
end
