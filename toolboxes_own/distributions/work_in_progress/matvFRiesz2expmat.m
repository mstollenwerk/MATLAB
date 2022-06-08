function M = matvFRiesz2expmat( df_1, df_2 )
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 06.09.2021

if length(df_1) ~= length(df_2)
    error('df_1 and df_2 must be same length.')
end

k = length(df_1);
a = NaN(k,1);

a(k) = df_1(k)/(df_2(k) - k - 1);

for ii = fliplr(1:k-1)
    a(ii) = (df_1(ii) + sum(a(ii+1:k)))/(df_2(ii) - ii - 1);
end

M = diag(a);

end
