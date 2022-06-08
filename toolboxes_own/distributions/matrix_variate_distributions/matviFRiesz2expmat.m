function A = matviFRiesz2expmat( n, nu )


% Michael Stollenwerk
% michael.stollenwerk@live.com
% 28.04.2022

if length(n) ~= length(nu)
    error('n and nu must be same length.')
end

M = matviRiesz2expmat(nu);
m = diag(M);

p = length(nu);
a = NaN(p,1);

a(1) = n(1)*m(1);

for ii = 2:p
    a(ii) = sum(m(1:ii-1)) + (n(ii) - ii + 1)*m(ii);
end

A = diag(a);

end
