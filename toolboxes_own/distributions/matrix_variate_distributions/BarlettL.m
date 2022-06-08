function BL = BarlettL(df)
%BarlettL Random lower Barlett decomposition matrix.
%   Detailed explanation goes here

p = length(df);
p_ = p*(p+1)/2 - p;

chi2_l = NaN(p,1);
for pp = 1:p
    chi2_l(pp) = sqrt(chi2rnd( df(pp) - pp + 1 ));
end
BL = diag(chi2_l);
BL(tril(true(p),-1)) = randn(p_,1);

end

