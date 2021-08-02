function vecMatrix = vec3d(Matrix)
% Vectorizes a matrix
% 
%% Inputs/Outputs
%
% Inputs:
%   MATRIXDATA   - a K by K symmetric matrix
%
% Outputs:
%   STACKEDDATA - A K^2 vector of stacked data
%
% See also vech 

%% Input Checking
[p,~,n] = size(Matrix);
vecMatrix = NaN(n,p^2);
for ii = 1:n
    vecMatrix(ii,:) = vec(Matrix(:,:,ii));
end
end
