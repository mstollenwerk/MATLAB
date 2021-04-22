

proc(2)=skew(x,prewhite,kernel);
local xbar,m3,se,sig,n,a,m2,m4,m6,std,z,stat,omega,stat1;
a=ones(2,1);
xbar=meanc(x);
z=x-xbar;
n=rows(x);
m2=sumc(z^2)/(n-1);
m3=sumc(z^3)/(n-1);
m4=sumc(z^4)/(n-1);
m6=sumc(z^6)/(n-1);
omega=zeros(2,2);
/*omega[1,1]=m6;omega[2,1]=m4;omega[1,2]=m4;omega[2,2]=m2;*/
if kernel==1; omega=parzen(z^3~z,prewhite); endif;
if kernel==2; omega=nw(z^3~z,prewhite); endif;
if kernel==3; omega=qs(z^3~z,prewhite); endif;
a[2]=-3*m2;
/*omega=kernel(z^3~z,4);*/
/*se=sqrt(a'omega*a/m2^3)/sqrt(n);*/
stat1=m3/m2^(3/2);
se=sqrt(a'omega*a)/sqrt(n);
stat=m3;
retp(stat/se,stat1);
endp;


