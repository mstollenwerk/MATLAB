function [ eparam, minval ] = ab_gridsearch( msh_sz, fun )
%By Walter Robinson on Matlab answers

a_ = 0:msh_sz:1;  %list of places to search for first parameter
b_ = 0:msh_sz:1;  %list of places to search for second parameter
[A_,B_] = ndgrid(a_, b_);
    ind = ((A_+B_)<1 & A_ <.2 & B_>.7) ;
    warning('Assuming a+b < 1 & a < .2 & b > .7')
    A_ = A_(ind);
    B_ = B_(ind);
fitresult = arrayfun(@(a,b) fun([a,b]), A_, B_); %run a fitting on every pair fun(F(J,K), S(J,K))
[minval, minidx] = min(fitresult);
a = A_(minidx);
b = B_(minidx);
eparam = [a,b];
end