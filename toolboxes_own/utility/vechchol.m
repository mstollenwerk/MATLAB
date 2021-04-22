function [y] = vechchol(sympdMatrix)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
y = vech(chol(sympdMatrix,'lower'));
end

