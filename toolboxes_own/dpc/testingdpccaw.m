function [eparam, logL, elogL]=testingdpccaw()

rep=1;
n=5;

param=[2.29880857373819,1.38662371087377,0.888061259380541,0.811659185461097,0.610045611773473,2.72757681110731,0.459090434626201,0.381994272087158,0.277217739028575,1.66613927533985,0.426449826679048,0.369486035781381,1.69572216069270,0.282785535017920,1.27957205996982,0.05,0.06,0.07,0.08,0.09,0.9,0.89,0.88,0.87,0.86,0.25,0.3,0.35,0.4,0.45,0.7,0.65,0.6,0.55,0.5,21];
nu=param(end);

logL=NaN(rep,1);
enu=NaN(rep,1);
elogL=NaN(rep,1);
eparam=NaN(rep,n*(n+1)/2+4*n);
for i=1:rep
    [ R, S, ~, ~ ] = dpccaw_sim( param, 1000, n, 'restrfull' );
    [ ~, LogLcontr_2ndpart ] = cawqlike( S, R );
    [ nLogL, ~ ] = cawlike_givn2ndpart( nu, R, LogLcontr_2ndpart );
    logL(i)=-nLogL;
    [ x0, ~, ~, ~, ~, ~ ] = dpc_ltarg( R, 1, 1, [] );
    [Lc,~]=dpceig(mean(R,3));
    Sc=Lc*diag(x0(1:n))*Lc';
    [ eparam(i,:), elogL(i), ~, ~, ~, enu(i) ] = dpc( R, 1, 1, [vech(chol(Sc,'lower'))' x0(n+1:end)] );
end
figure
plot(elogL-logL)