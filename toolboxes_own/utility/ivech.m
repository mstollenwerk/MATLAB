function Matrix=ivech(vechMatrix,type)
%% Transform a vector into 'sym' symmetric or 'lower' lower triangluar matrix
% 
% USAGE:
%   MATRIX = ivech(vechMatrix)
% 
% INPUTS:
%   vechMatrix   - A K(K+1)/2 vector of data to be transformed 
%   type         - 1 for symmetric, 2 for lower triangular
%
% OUTPUTS:
%   Matrix       - 'sym' K by K symmetric matrix of the form 
%                  [ data(1) data(2)    data(3)     ...               data(K)
%                    data(2) data(K+1)  data(K+2)   ...               ...
%                    data(3) data(K+2)  data(2K)    ...               ...
%                    ...     ....       ...         ...               data(K(K+1)/2-1)
%                    data(K) data(2K-1) ...         data(K(K+1)/2-1)  data(K(K+1)/2) ]
%
%                  OR 'lower' K by K lower triangular matrix of the form 
%                  [ data(1) 0          0           ...               0
%                    data(2) data(K+1)  0           ...               0
%                    data(3) data(K+2)  data(2K)    ...               ...
%                    ...     ....       ...         ...               0
%                    data(K) data(2K-1) ...         data(K(K+1)/2-1)  data(K(K+1)/2) ]

% Author: Kevin Sheppard
% Modified: Michael Stollenwerk
% Modification: Added option to create lower triangular matrix.
if nargin==1
    type='sym';
end
%% Input Checking
if size(vechMatrix,2)>size(vechMatrix,1)
    vechMatrix=vechMatrix';
end

if size(vechMatrix,2)~=1
    error('STACKED_DATA must be a column vector.')
end

K2=size(vechMatrix,1);%checking if numel in stackedData fits
K=(-1+sqrt(1+8*K2))/2; %checking if numel in stackedData fits

if floor(K)~=K %checking if numel in stackedData fits
    error(['The number of elements in STACKED_DATA must be conformable to' ...
    'the inverse chol2vec operation.'])
end

if ~strcmp(type,'sym') && ~strcmp(type,'lower')
    error('Please specify type: sym or lower')
end

%%
% Initialize the output data
Matrix=zeros(K);

% Use a logical trick to inverse the vech
pl=tril(true(K));
Matrix(pl)=vechMatrix;

if strcmp(type,'sym')
    diag_matrixData=diag(diag(Matrix));
    Matrix=Matrix+Matrix'-diag_matrixData;
end
end
