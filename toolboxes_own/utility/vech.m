function vechMatrix = vech(Matrix)
% Half-vectorizes a matrix
% 
%% Inputs/Outputs
%
% Inputs:
%   MATRIXDATA   - a K by K symmetric matrix
%
% Outputs:
%   STACKEDDATA - A K(K+1)/2 vector of stacked data
%
% Comments:
%   The data is stacked according to 
%     [ data(1) ...        ...         ...               ...
%       data(2) data(K+1)  ...         ...               ...
%       data(3) data(K+2)  data(2K)    ...               ...
%       ...     ....       ...         ...               ...
%       data(K) data(2K-1) ...         data(K(K+1)/2-1)  data(K(K+1)/2) ]

%% Input Checking
[k,l] = size(Matrix);
pl = ~tril(true(k));
if k~=l 
    error('MATRIXDATA must be a square matrix');
end
if ~issymmetric(Matrix) && any(Matrix(pl)~=0) && ~all(isnan(Matrix(:)))
    error('MATRIXDATA must be a lower triangular matrix or symmetric');
end

%% Core Code
sel = tril(true(k));
vechMatrix = Matrix(sel);
end
