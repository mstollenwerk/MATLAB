% this function generates a random unitary matrix of order 'n' and verifies 
function [U,verify]= rndunitary(n) 
% generate a random complex matrix  
X = rand(n)/sqrt(2); 
% factorize the matrix 
[Q,R] = qr(X); 
R = diag(diag(R)./abs(diag(R))); 
% unitary matrix 
U = Q*R; 
% verification 
verify = U*U'; 
end 