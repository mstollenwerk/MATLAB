function Y = trilSym(X)
%TRILHALFDIAG Lower triangular part of matrix with its diagonals halfed.

Y = tril(X) + tril(X,-1)';

end

