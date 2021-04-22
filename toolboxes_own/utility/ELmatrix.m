function [ L ] = ELmatrix( n )
%Ematrix creates Elimination matrix as in Magnus and Neudecker (1988)
% Follows: Mathematics for Econometrics, p.131, Phoebus J. Dhrymes (2013)

storage = cell(n,1);
for i=1:n
    storage(i)= {[zeros(n+1-i,i-1), eye(n+1-i)]};
end

L=blkdiag(storage{:});  
  
end