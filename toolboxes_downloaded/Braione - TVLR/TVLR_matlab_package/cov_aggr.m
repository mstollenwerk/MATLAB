function [mat]=cov_aggr(C)
%==========================================================================
%This function stacks the unique part of a (nxn) matrix into a (n(n+1)/2)
%row vector of parameters.
%
% INPUT
% C        =series of realized covariance matrices
%
% OUTPUT
% mat      = [T x (n(n+1)/2)] matrix with daily matrices in vech form per row
%==========================================================================
% Author : Manuela Braione, CORE-UCL 2016
%==========================================================================

[n, n, T]=size(C);

mat = zeros(T,n*(n+1)/2);

for t=1:T
    CC=vech(C(:,:,t));
    mat(t,:)=CC';
end





