proc(3)=arskew(dy,p);
local reg,i,b,e,fit;
local vnam,m,stb,vc,stderr,sigma,cx,rsq,dw;
reg=ones(rows(dy)-p,1);
i=1;
do while i<=p;
reg=reg~dy[p-i+1:rows(dy)-i];
i=i+1;
endo;
if p > 0;
__con=0;
_olsres=1;
output off;
{vnam,m,b,stb,vc,stderr,sigma,cx,rsq,e,dw}=
ols(0,dy[p+1:rows(dy)],reg);
output on;
{b,e,fit}=olsqr2(dy[p+1:rows(dy)],reg);
else;
{b,e,fit}=olsqr2(dy[p+1:rows(dy)],reg);
dw=sumc((e[2:rows(e)]-e[1:rows(e)-1])^2)/sumc(e[1:rows(e)-1]^2);
endif;
retp(e,dw,reg);
endp;

