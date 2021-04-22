function Y = Cmatrix(m,n)
% Calculate commutation matrix as in Magnus Neuendecker 2005.
% Code take from stackexchange. 
% Input m = nrows, n = ncols.
I = reshape(1:m*n, [m, n]); % initialize a matrix of indices of size(A)
I = I'; % Transpose it
I = I(:); % vectorize the required indices
Y = eye(m*n); % Initialize an identity matrix
Y = Y(I,:); % Re-arrange the rows of the identity matrix
end