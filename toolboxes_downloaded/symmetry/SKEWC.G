

proc(1)=skewc(x,e,p,prewhite);
local xbar,m3,se,sig,n,a,m2,m4,m6,std,z,stat,omega;
a=1;
xbar=meanc(x);
z=x.*e;
n=rows(x);
omega=lrhac(e^3~z,prewhite);
a=a|-3*inv(x'x/n)*meanc(e.*e.*x);
/*omega=kernel(z^3~z,4);*/
se=sqrt(a'omega*a)/sqrt(n);
stat=meanc(e^3);
retp(stat/se);
endp;


