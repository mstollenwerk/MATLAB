function PDmatrix = makepd(nonPDmatrix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N = size(nonPDmatrix,3);

PDmatrix = NaN(size(nonPDmatrix));

for ii = 1:N
    [X,Y] = eig(nonPDmatrix(:,:,ii));
    if any(diag(Y) < 0)
        warning('pd adjustment made')
    end
    Y(Y<0) = eps^(1/3);
    y = X*Y*X';
    PDmatrix(:,:,ii) = .5*(y+y');
end
    
end

