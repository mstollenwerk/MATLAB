function vechMatrix = vech_(Matrix)
% Half-vectorizes a matrix
% 
%% Inputs/Outputs
%
% Inputs:
%   MATRIXDATA   - A K by K by T symmetric matrix array
%
% Outputs:
%   STACKEDDATA  - A T by K(K+1)/2 vector of stacked data
%
% Comments:
%   The data in row t (element of 1:T) is stacked according to 
%     [ data(1) ...        ...         ...               ...
%       data(2) data(K+1)  ...         ...               ...
%       data(3) data(K+2)  data(2K)    ...               ...
%       ...     ....       ...         ...               ...
%       data(K) data(2K-1) ...         data(K(K+1)/2-1)  data(K(K+1)/2) ]

%% Input Checking
[k,~,T] = size(Matrix);
vechMatrix = NaN(T,k*(k+1)/2);
for tt = 1:T
    vechMatrix(tt,:) = vech(Matrix(:,:,tt));
end
end
