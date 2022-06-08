function BU = BarlettU(df)
%BarlettU Random upper Barlett decomposition matrix.
%   Detailed explanation goes here

p = length(df);
p_ = p*(p+1)/2 - p;

chi2_u = NaN(p,1);
for pp = 1:p
    chi2_u(pp) = sqrt(chi2rnd( df(pp) + pp - p ));
end
BU = diag(chi2_u);
BU(triu(true(p),1)) = randn(p_,1);

end

