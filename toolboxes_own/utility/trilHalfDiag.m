function Y = trilHalfDiag(X)
%TRILHALFDIAG Lower triangular part of matrix with its diagonals halfed.

Y = tril(X) - .5.*diag(diag(X));

end

