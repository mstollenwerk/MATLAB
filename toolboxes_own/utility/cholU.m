function [U] = cholU(X)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here
U = inv(chol(inv(X),'lower'))';
end

