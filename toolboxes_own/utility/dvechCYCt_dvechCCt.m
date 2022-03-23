function y = dvechCYCt_dvechCCt(C,Y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = size(C,1);

[~,Gplus] = Dmatrix(p);
F = ELmatrix(p);

y = Gplus*kron(C*Y,eye(p))*F'/(Gplus*kron(C,eye(p))*F');

end

