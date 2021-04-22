function [ B ] = kronroot( BkronB )
%kronroot calculates the "kronecker root"
%
%%Inputs/Outputs
%
%%%%%%%%
%Inputs%    BkronB - The Kronecker product of B with itself
%%%%%%%%
%
%%%%%%%%%
%Outputs%        B - The "root matrix" B
%%%%%%%%%
%
%%
n=sqrt(size(BkronB,1));

% Create logical wich picks out the squared elements of B in BkronB
pattern = [1 zeros(1,n)];
rowcolpattern = [repmat(pattern,1,n-1) 1];
pickout = logical(rowcolpattern'*rowcolpattern);

% Pick out the squared elements of B in BkronB and take their square roots
B = reshape(BkronB(pickout).^.5,n,n);
end