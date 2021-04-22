function [D, EL] = Dmatrix(n)
%DELmat Creates Duplication and Elimination Matrix as in Magnus and
%Neudecker (1988).
%
% USAGE:
%  [D,EL] = Dmatrix(n)
%
% INPUTS:
%   N            - Matrix Dimension
%
%
% OUTPUTS:
%   D            - Duplication Matrix
%   EL           - Elimination Matrix
%
% COMMENTS:
%
%  See also 
%
% REFERENCES:
%      [1] \url https://code.google.com/p/addi/source/browse/Addi/assets/m/linear-algebra/duplication_matrix.m?r=190
%      [2] 

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 25.06.2018

% DEPENDENCIES:

%% Input Checks

%% Duplication Matrix
% Not my code. See Reference [1]
D = zeros (n * n, n * (n + 1) / 2);
count = 0;
for j = 1 : n
    D ((j - 1) * n + j, count + j) = 1;
    for i = (j + 1) : n
      D ((j - 1) * n + i, count + i) = 1;
      D ((i - 1) * n + j, count + i) = 1;
   end
   count = count + n - j;
end
%% Elimination Matrix
% The Moore-Penrose inverse is an Elimination matrix, but not the canonical
% one! E = (D'*D)\D'

end