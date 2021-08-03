function [eDetPower_h] = ncWish_Edet(df_, Sigma_, Omega_, h)
%NCWISH_EDET Calculates expected value of determinant of ncWish to the
%power of h.
%   Theorem 3.5.6 Gupta, Nagar (2000) - Matrix Variate Distributions.

p = size(Sigma_,1);
eDetPower_h = 2^(p*h)*mvgamma(.5*df_ + h, p) * det(Sigma_)^h /...
    mvgamma(.5*df_, p) * exp(trace(-.5*Omega_)) * ...
    mhg(25, 2, .5*df_ + h, .5*df_, eig(.5*Omega_));
end

