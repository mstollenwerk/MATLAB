function [A]=createA(n,intercept)
%
%This function creates a matrix containing the bounds to be imposed on the
% model parameters (used as fmincon input) as function of the number of
% assets n.
%--------------------------------------------------------------------------
% INPUT:
% n             - number of assets
% intercept     - intercept matrix in vech() form
%==========================================================================
% Author : Manuela Braione, CORE-UCL 2016
%==========================================================================

a=eye(n);
b=zeros(n,length(intercept)+2);
A=[b a a];
A(:,end+1)=0;
A(:,end+1)=0;
A(end+1,end-1:end)=1;

end