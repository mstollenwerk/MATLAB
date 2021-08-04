function y = laplacepdf( theta_, s, dta )
%LAPLACELIKE p.d.f. of univariate Laplace distribution
%
% REFERENCES:
%      [1] Kotz, Kozubowski and Podgorski (2001), p. 16, 2.1.1

narginchk(3,3);
nargoutchk(0,1);

if size(dta,1) < size(dta,2)
    dta = dta';
end

y = exp(-abs(dta - theta_)/s)/2/s;
end
