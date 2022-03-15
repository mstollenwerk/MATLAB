function x=ivech3d(vechMatrix,type)
%% Transform a vector into 'sym' symmetric or 'lower' lower triangluar matrix
% 
% USAGE:
%   MATRIX = ivech(vechMatrix)
% 
% INPUTS:
%   vechMatrix   - A K(K+1)/2 by N vector of data to be transformed 
%   type         - 1 for symmetric, 2 for lower triangular
%
% OUTPUTS:
%   Matrix       - 'sym' K by K by N symmetric matrix of the form 
%                  [ data(1) data(2)    data(3)     ...               data(K)
%                    data(2) data(K+1)  data(K+2)   ...               ...
%                    data(3) data(K+2)  data(2K)    ...               ...
%                    ...     ....       ...         ...               data(K(K+1)/2-1)
%                    data(K) data(2K-1) ...         data(K(K+1)/2-1)  data(K(K+1)/2) ]
%
%                  OR 'lower' K by K by N lower triangular matrix of the form 
%                  [ data(1) 0          0           ...               0
%                    data(2) data(K+1)  0           ...               0
%                    data(3) data(K+2)  data(2K)    ...               ...
%                    ...     ....       ...         ...               0
%                    data(K) data(2K-1) ...         data(K(K+1)/2-1)  data(K(K+1)/2) ]

% Author: Kevin Sheppard
% Modified: Michael Stollenwerk
% Modification: Added option to create lower triangular matrix.
% 25.05.2020: Changed input format of vechMatrix to first dim being
% different matrices and second dim being elements of matrices.

if nargin==1
    type='sym';
end

%% Input Checking
[N,K_] = size(vechMatrix);
K = -.5 + sqrt(.25+2*K_);
x = NaN(K,K,N);
for tt = 1:N
    x(:,:,tt) = ivech(vechMatrix(tt,:),type);
end
end
