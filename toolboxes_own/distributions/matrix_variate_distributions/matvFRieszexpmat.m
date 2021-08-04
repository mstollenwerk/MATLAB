function M = matvFRieszexpmat( df_1, df_2 )
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 16.02.2021

if length(df_1) ~= length(df_2)
    error('df_1 and df_2 must be same length.')
end

k = length(df_1);
a = NaN(k,1);

a(1) = df_1(1)/(df_2(1) - k - 1);

for ii = 2:k
    a(ii) = (df_1(ii) + sum(a(1:(ii-1))))/(df_2(ii) - k + ii - 2);
end

M = diag(a);

end
