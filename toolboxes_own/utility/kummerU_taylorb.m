function y = kummerU_taylorb(a,b,z,tol)

y = ...
    gamma(1 - b)/gamma(a - b + 1)*gamma(b)*hypfun_M_taylorb(a,b,z,tol) ...
    + gamma(b - 1)/gamma( a )*z^(1 - b)*gamma(2-b)*hypfun_M_taylorb(a - b + 1, 2 - b, z ,tol);

end

