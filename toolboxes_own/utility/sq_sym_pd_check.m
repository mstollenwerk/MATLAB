function [square,sym,pd] = sq_sym_pd_check(mat_)
square=1;
sym=1;
pd=1;
if mat_~=mat_'
    sym=0;
    warning('Matrix is not symmetric')
end
if size(mat_,1)~=size(mat_,2) || ~ismatrix(mat_)
    square=0;
    warning('Matrix is not square')
end
if min(eig(mat_))<0
    pd=0;
    warning('Matrix is not pd')
end
end
    