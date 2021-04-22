function [symH] = symHessian(H)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here
p = size(H,1);

D = Dmatrix(sqrt(p));

symH = D*D'*H*(D*D'); % see https://tminka.github.io/papers/matrix/minka-matrix.pdf Formula (111)
end

