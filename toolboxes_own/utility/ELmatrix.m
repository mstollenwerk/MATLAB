function [ L ] = ELmatrix( n )
%Ematrix creates Elimination matrix as in Magnus and Neudecker (1988)

% % Follows: Mathematics for Econometrics, p.131, Phoebus J. Dhrymes (2013)
% storage = cell(n,1);
% for i=1:n
%     storage(i)= {[zeros(n+1-i,i-1), eye(n+1-i)]};
% end
% 
% L=blkdiag(storage{:});  

persistent C_ELmatrix
if isempty(C_ELmatrix)
   nCache = 1000;  % Set according to your needs
   C_ELmatrix      = cell(1, nCache);
end
if n <= numel(C_ELmatrix) && ~isempty(C_ELmatrix{n})
    L = C_ELmatrix{n};
else
    n_ = n*(n+1)/2;
    IndexMatrix = reshape(1:n^2,n,n);
    select = tril(true(n));
    v = IndexMatrix(select);

    L = sparse(1:n_, v, 1, n_, n^2);
    if n <= numel(C_ELmatrix)
      C_ELmatrix{n} = L;
    end
end

end