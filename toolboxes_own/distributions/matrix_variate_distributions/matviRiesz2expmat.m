function M = matviRiesz2expmat( nu )
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 09.03.2021

if size(nu,1) < size(nu,2)
    nu = nu';
end

p = length(nu);
d = NaN(p,1);

d(1) = 1/(nu(1)-p-1);

for ii = 1:p
    d(ii) = 1/(nu(ii)-p+ii-2)*(1+sum(d(1:ii-1)));
end
    
M = spdiags(d,0,p,p);

end
