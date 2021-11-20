function PDmatrix = makepd(nonPDmatrix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[X,Y] = eig(nonPDmatrix);
Y(Y<0) = eps^(1/3);
PDmatrix = X*Y*X';
PDmatrix = .5*(PDmatrix+PDmatrix');

end

