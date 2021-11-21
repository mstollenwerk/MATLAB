function [PDmatrix, varargout] = makepd(nonPDmatrix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N = size(nonPDmatrix,3);

PDmatrix = NaN(size(nonPDmatrix));

adjustments_made = zeros(N,1);
for ii = 1:N
    [X,Y] = eig(nonPDmatrix(:,:,ii));
    if any(diag(Y) < 0)
        Y(Y<0) = eps^(1/3);
        y = X*Y*X';
        PDmatrix(:,:,ii) = .5*(y+y');
        adjustments_made(ii) = sum(diag(Y) < 0);
    else
        PDmatrix(:,:,ii) = nonPDmatrix(:,:,ii);
    end
end

if nargout == 2
    varargout{1} = adjustments_made;
end
    
end

