function M = matviRiesz2expmat( n )
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 09.03.2021

if size(n,1) < size(n,2)
    n = n';
end

p = length(n);
d = NaN(p,1);

d(1) = 1/(n(1)-p-1);

for ii = 1:p
    d(ii) = 1/(n(ii)-p+ii-2)*(1+sum(d(1:ii-1)));
end
    
M = diag(d);

end
