function M = matviRieszexpmat( n )
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 09.03.2021

if size(n,1) < size(n,2)
    n = n';
end

p = length(n);

d(p) = 1/(n(p) - (p+1));
for ii = fliplr(1:p-1)
    d(ii) = (n(ii+1) - (ii+1))/(n(ii) - (ii+1))*d(ii+1);
end
    
M = diag(d);

end
