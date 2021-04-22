function [covMat, zmax] = sym_covmat(covMat_almost_sym)

T = size(covMat_almost_sym,3);

zmax = NaN(T,1);
covMat = NaN(size(covMat_almost_sym));
for t = 1:T
    zmax(t) = max(max(abs(covMat_almost_sym(:,:,t)-covMat_almost_sym(:,:,t)')));
if  zmax(t) < 0.000001
    covMat(:,:,t)=(covMat_almost_sym(:,:,t)+covMat_almost_sym(:,:,t)')/2;
elseif ~isnan(zmax(t))
    warning(['Abs deviation of an element i,j from element j,i in covMat_almost_sym(:,:,' num2str(t) ') is' num2str(zmax(t))])
end
end