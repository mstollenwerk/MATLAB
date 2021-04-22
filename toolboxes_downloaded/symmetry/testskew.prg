

new;
library pgraph;
format /ld 6,3;
outwidth 132;
prewhite=0;
kernel=2;

load gdpdef[154,2]=e:\gauss\unsymm\skew\gdpdef.asc;
load gdp[154,2]=e:\gauss\unsymm\skew\gdp92.asc;
load un[600,3]=e:\gauss\unsymm\skew\unrate.asc;
load ind[624,3]=e:\gauss\unsymm\skew\indpro.asc;
load ca[324,3]=e:\gauss\unsymm\skew\excaus.asc;
load ger[324,3]=e:\gauss\unsymm\skew\exgeus.asc;
load jp[324,3]=e:\gauss\unsymm\skew\exjpus.asc;
load cp30[759,2]=e:\gauss\unsymm\skew\cp30.asc;
load m2[899,2]=e:\gauss\unsymm\skew\m2.asc;
load cpi[612,3]=e:\gauss\unsymm\skew\cpi.asc;
load dur[150,2]=e:\gauss\unsymm\skew\dgc92.asc;
load ndur[150,2]=e:\gauss\unsymm\skew\ndgc92.asc;
load invest[154,3]=e:\gauss\unsymm\skew\gpdic92.asc;
load emp[587,2]=e:\gauss\unsymm\skew\em1664.asc;
load manemp[624,3]=e:\gauss\unsymm\skew\manemp.txt;
load nmfemp[611,2]=e:\gauss\unsymm\skew\nmfemp.txt;
load sales[154,3]=e:\gauss\unsymm\skew\finslc92.txt;
load invnres[154,3]=e:\gauss\unsymm\skew\pnfic92.txt;
load invres[154,3]=e:\gauss\unsymm\skew\prfic92.txt;
y1=ln(ca[.,3]);
y2=ln(ger[.,3]);
y3=ln(jp[.,3]);
y4=un[.,3];
y5=log(ind[.,3]);
y6=log(gdpdef[.,2]);
y7=log(gdp[.,2]);
y8=cp30[531:759,2];
y9=log(m2[583:811,2]);
y10=log(cpi[169:612.,3]);
y11=log(dur[.,2]);
y12=log(ndur[.,2]);
y13=log(emp[157:587,2]);
y14=log(invest[.,3]);
y15=log(manemp[181:624,3]);
y16=log(nmfemp[181:611,2]);
y17=log(sales[.,3]);
y18=log(invnres[.,3]);
y19=log(invres[.,3]);

stat1=zeros(21,1);
stat1a=zeros(21,1);
stat2=zeros(21,1);
stat3=zeros(21,3);
stat4=zeros(21,2);

output file=testskew.out reset;
y=y1;
p=0;
dy=y[2:rows(y)]-y[1:rows(y)-1];
{stat1[1,1],stat1a[1,1]}=skew(dy,prewhite,kernel);
stat2[1]=skew35(dy,prewhite,kernel);
stat3[1,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[1,1]=normal(dy,prewhite,kernel);
stat4[1,2]=stat1[1,1]^2+stat3[1,2]^2;


print "CA-US exchange rate";

y=y2;
p=0;
dy=y[2:rows(y)]-y[1:rows(y)-1];
{stat1[2,1],stat1a[2,1]}=skew(dy,prewhite,kernel);
stat2[2]=skew35(dy,prewhite,kernel);
stat3[2,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[2,1]=normal(dy,prewhite,kernel);
stat4[2,2]=stat1[2,1]^2+stat3[2,2]^2;
print "GER-US exchange rate";


y=y3;
p=0;
dy=y[2:rows(y)]-y[1:rows(y)-1];
{stat1[3,1],stat1a[3,1]}=skew(dy,prewhite,kernel);
stat2[3]=skew35(dy,prewhite,kernel);
stat3[3,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[3,1]=normal(dy,prewhite,kernel);
stat4[3,2]=stat1[3,1]^2+stat3[3,2]^2;
print "JP-US exchange rate";


y=y4;
p=2;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "Unemployment";
{stat1[4,1],stat1a[4,1]}=skew(dy,prewhite,kernel);
stat2[4]=skew35(dy,prewhite,kernel);
stat3[4,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[4,1]=normal(dy,prewhite,kernel);
stat4[4,2]=stat1[4,1]^2+stat3[4,2]^2;


y=y5;
p=2;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "Industrial Production";
{stat1[5,1],stat1a[5,1]}=skew(dy,prewhite,kernel);
stat2[5]=skew35(dy,prewhite,kernel);
stat3[5,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[5,1]=normal(dy,prewhite,kernel);
stat4[5,2]=stat1[5,1]^2+stat3[5,2]^2;


y=y6;
p=2;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "GDP deflator";
{stat1[6,1],stat1a[6,1]}=skew(dy,prewhite,kernel);
stat2[6]=skew35(dy,prewhite,kernel);
stat3[6,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[6,1]=normal(dy,prewhite,kernel);
stat4[6,2]=stat1[6,1]^2+stat3[6,2]^2;


y=y7;
p=2;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "GDP ";
{stat1[7,1],stat1a[7,1]}=skew(dy,prewhite,kernel);
stat2[7]=skew35(dy,prewhite,kernel);
stat3[7,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[7,1]=normal(dy,prewhite,kernel);
stat4[7,2]=stat1[7,1]^2+stat3[7,2]^2;


y=y10;
p=2;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "CPI";
{stat1[8,1],stat1a[8,1]}=skew(dy,prewhite,kernel);
stat2[8]=skew35(dy,prewhite,kernel);
stat3[8,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[8,1]=normal(dy,prewhite,kernel);
stat4[8,2]=stat1[8,1]^2+stat3[8,2]^2;


y=y8;
p=2;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "30 day interest rate";
{stat1[9,1],stat1a[9,1]}=skew(dy,prewhite,kernel);
stat2[9]=skew35(dy,prewhite,kernel);
stat3[9,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[9,1]=normal(dy,prewhite,kernel);
stat4[9,2]=stat1[9,1]^2+stat3[9,2]^2;


y=y9;
p=2;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "M2";
{stat1[10,1],stat1a[10,1]}=skew(dy,prewhite,kernel);
stat2[10]=skew35(dy,prewhite,kernel);
stat3[10,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[10,1]=normal(dy,prewhite,kernel);
stat4[10,2]=stat1[10,1]^2+stat3[10,2]^2;






y=y11;
p=2;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "Durables";
{stat1[11,1],stat1a[11,1]}=skew(dy,prewhite,kernel);
stat2[11]=skew35(dy,prewhite,kernel);
stat3[11,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[11,1]=normal(dy,prewhite,kernel);
stat4[11,2]=stat1[11,1]^2+stat3[11,2]^2;


y=y12;
p=2;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "Non-Durables";
{stat1[12,1],stat1a[12,1]}=skew(dy,prewhite,kernel);
stat2[12]=skew35(dy,prewhite,kernel);
stat3[12,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[12,1]=normal(dy,prewhite,kernel);
stat4[12,2]=stat1[12,1]^2+stat3[12,2]^2;


y=y13;
p=2;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "Employment";
{stat1[13,1],stat1a[13,1]}=skew(dy,prewhite,kernel);
stat2[13]=skew35(dy,prewhite,kernel);
stat3[13,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[13,1]=normal(dy,prewhite,kernel);
stat4[13,2]=stat1[13,1]^2+stat3[13,2]^2;


y=y14;
p=4;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "investment";
{stat1[14,1],stat1a[14,1]}=skew(dy,prewhite,kernel);
stat2[14]=skew35(dy,prewhite,kernel);
stat3[14,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[14,1]=normal(dy,prewhite,kernel);
stat4[14,2]=stat1[14,1]^2+stat3[14,2]^2;


y=y15;
p=2;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "Manu. emp.";
{stat1[15,1],stat1a[15,1]}=skew(dy,prewhite,kernel);
stat2[15]=skew35(dy,prewhite,kernel);
stat3[15,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[15,1]=normal(dy,prewhite,kernel);
stat4[15,2]=stat1[15,1]^2+stat3[15,2]^2;


y=y16;
p=2;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "Non-man. empl";
{stat1[16,1],stat1a[16,1]}=skew(dy,prewhite,kernel);
stat2[16]=skew35(dy,prewhite,kernel);
stat3[16,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[16,1]=normal(dy,prewhite,kernel);
stat4[16,2]=stat1[16,1]^2+stat3[16,2]^2;


y=y17;
p=2;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "final sales ";
{stat1[17,1],stat1a[17,1]}=skew(dy,prewhite,kernel);
stat2[17]=skew35(dy,prewhite,kernel);
stat3[17,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[17,1]=normal(dy,prewhite,kernel);
stat4[17,2]=stat1[17,1]^2+stat3[17,2]^2;

y=y18;
p=2;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "Non-res invest";
{stat1[18,1],stat1a[18,1]}=skew(dy,prewhite,kernel);
stat2[18]=skew35(dy,prewhite,kernel);
stat3[18,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[18,1]=normal(dy,prewhite,kernel);
stat4[18,2]=stat1[18,1]^2+stat3[18,2]^2;


y=y19;
p=2;
dy=y[2:rows(y)]-y[1:rows(y)-1];
print "Rest. invest";
{stat1[19,1],stat1a[19,1]}=skew(dy,prewhite,kernel);
stat2[19]=skew35(dy,prewhite,kernel);
stat3[19,1:3]=kurtosis(dy,5,prewhite,kernel);
stat4[19,1]=normal(dy,prewhite,kernel);
stat4[19,2]=stat1[19,1]^2+stat3[19,2]^2;


load stocks[8685,7]=e:\gauss\unsymm\skew\stockret.txt;
y20=stocks[6915:8685,2];
y20=stocks[7421:8685,2];
y=y20;
{stat1[20,1],stat1a[20,1]}=skew(y,prewhite,kernel);
stat2[20]=skew35(y,prewhite,kernel);
stat3[20,1:3]=kurtosis(y,5,prewhite,kernel);
stat4[20,1]=normal(y,prewhite,kernel);
stat4[20,2]=stat1[20,1]^2+stat3[20,2]^2;
print "Stock 1";

y21=stocks[6915:8665,4];
y21=stocks[7421:8665,4];
y=y21;
{stat1[21,1],stat1a[21,1]}=skew(y,prewhite,kernel);
stat2[21]=skew35(y,prewhite,kernel);
stat3[21,1:3]=kurtosis(y,5,prewhite,kernel);
stat4[21,1]=normal(y,prewhite,kernel);
stat4[21,2]=stat1[21,1]^2+stat3[21,2]^2;
y=y21;
p=0;
print "Stock 2";

output file=testskew.out reset;
print "skewness ";
stat1~stat1a~stat2;
print;
print "kurtosis";
stat3;
print;
print "normality";
stat4;
output off;

stop;
y22=stocks[6915:8685,7];
y22=stocks[7421:8685,7];
y=y22;
{stat1[22,1],stat1a[22,1]}=skew(y,1);
stat2[22]=skew35(y,1);
stat3[22,1:3]=kurtosis(y,5,1);
stat4[22,1]=normal(y,1);
stat4[22,2]=stat1[22,1]^2+stat3[22,2]^2;
p=2;
y=y22;
print "Stock 3";


