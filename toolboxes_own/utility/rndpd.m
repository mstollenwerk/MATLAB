function [rndpd_] = rndpd(p,n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
narginchk(1,2)
if nargin == 1
    n = 1;
end

rndpd_ = NaN(p,p,n);
for ii = 1:n
    x = randn(p);
    x = x*x';
    rndpd_(:,:,ii) = x;
end

