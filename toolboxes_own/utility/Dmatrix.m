function [D, iD] = Dmatrix(n)
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
%   iD           - Elimination Matrix based on pseudo-inverse of D.
%
% COMMENTS:
%
%  See also 
%
% REFERENCES:
%      [1] \url https://code.google.com/p/addi/source/browse/Addi/assets/m/linear-algebra/duplication_matrix.m?r=190
%      [2] \url https://de.mathworks.com/matlabcentral/answers/473737-efficient-algorithm-for-a-duplication-matrix?s_tid=srchtitle

% Michael Stollenwerk
% michael.stollenwerk@live.com
% 25.06.2018

% DEPENDENCIES:

%% Input Checks

%% Duplication Matrix
% % Not my code. See Reference [1], but WAAYY too slow compared to
% % Reference [2]!
% D = zeros (n * n, n * (n + 1) / 2);
% count = 0;
% for j = 1 : n
%     D ((j - 1) * n + j, count + j) = 1;
%     for i = (j + 1) : n
%       D ((j - 1) * n + i, count + i) = 1;
%       D ((i - 1) * n + j, count + i) = 1;
%    end
%    count = count + n - j;
% end

% Not my code. See Reference [2]
persistent C_Dmatrix
if isempty(C_Dmatrix)
   nCache = 100;  % Set according to your needs
   C_Dmatrix = cell(1, nCache);
end
if n <= numel(C_Dmatrix) && ~isempty(C_Dmatrix{n})
   D = C_Dmatrix{n};
else
   m   = n * (n + 1) / 2;
   nsq = n^2;
   r   = 1;
   a   = 1;
   v   = zeros(1, nsq);
   cn  = cumsum(n:-1:2);
   for i = 1:n
      v(r:r + i - 2) = i - n + cn(1:i - 1);
      r = r + i - 1;
      
      v(r:r + n - i) = a:a + n - i;
      r = r + n - i + 1;
      a = a + n - i + 1;
   end
   D = sparse(1:nsq, v, 1, nsq, m);
   if n <= numel(C_Dmatrix)
      C_Dmatrix{n} = D;
   end
end

%% Elimination Matrix
if nargout >= 2
    persistent C_iDmatrix
    if isempty(C_iDmatrix)
       C_iDmatrix = cell(1, nCache);
    end
    if n <= numel(C_iDmatrix) && ~isempty(C_iDmatrix{n})
       iD = C_iDmatrix{n};
    else
        % The Moore-Penrose inverse is an Elimination matrix, but not the canonical
        % one! E = (D'*D)\D'
        iD = (D'*D)\D';
        if n <= numel(C_iDmatrix)
          C_iDmatrix{n} = iD;
       end
    end
end