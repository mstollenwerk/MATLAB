function M = matvtRieszexpmat( n, nu )
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 09.03.2021
p = length(n);
M = spdiags(n,0,p,p)*nu/(nu-2);

end
