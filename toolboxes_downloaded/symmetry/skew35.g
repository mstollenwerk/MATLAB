proc(1)=skew35(x,prewhite,kernel);
local xbar,m3,se,sig,n,a,m2,m4,m5,m6,std,z,stat,omega;
a=zeros(2,3);
a[1,2]=1.0;
a[2,1]=1.0;
xbar=meanc(x);
z=x-xbar;
n=rows(x);
m2=sumc(z^2)/(n-1);
m3=sumc(z^3)/(n-1);
m4=sumc(z^4)/(n-1);
m5=sumc(z^5)/(n-1);
m6=sumc(z^6)/(n-1);
a[1,3]=-3*m2;
a[2,3]=-5*m4;
omega=zeros(2,2);
if kernel==1; omega=parzen(z^5~z^3~z,prewhite); endif;
if kernel==2; omega=nw(z^5~z^3~z,prewhite); endif;
if kernel==3; omega=qs(z^5~z^3~z,prewhite); endif;
/*omega=kernel(z^5~z^3~z,4);*/
se=a*omega*a';
stat=(m3~m5);
retp(stat*inv(se)*stat'*n);
endp;

