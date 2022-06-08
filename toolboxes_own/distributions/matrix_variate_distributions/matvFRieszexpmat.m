function M = matvFRieszexpmat( n, nu )


% Michael Stollenwerk
% michael.stollenwerk@live.com
% 16.02.2021

if length(n) ~= length(nu)
    error('df_1 and df_2 must be same length.')
end

k = length(n);
a = NaN(k,1);

a(1) = n(1)/(nu(1) - k - 1);

for ii = 2:k
    a(ii) = (n(ii) + sum(a(1:(ii-1))))/(nu(ii) - k + ii - 2);
end

M = diag(a);

end
