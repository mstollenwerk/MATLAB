function [symX] = symMatGradient(X)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here
p = size(X,1);
symX = 2*X - eye(p)*X;
end

