%
% This is called within MAIN.m for the purpose of estimating the model.
% Initial values and bouds are pre-imposed. 
%==========================================================================
% Author : Manuela Braione, CORE-UCL 2016
%==========================================================================


% Define parameters
imat=ones(n,n)*0.1;
for k=1:n
    imat(k,k)=0.2;
end
qq=vech(imat);
c=qq';

gamma =0.1*ones(1,n);
delta =0.8*ones(1,n);
alpha = 0.2;
beta  = 0.7;
teta  = 0.6;
omega = 9;
pars  =[c teta omega gamma delta alpha beta];

% Define bounds on parameters
g=zeros(1,length(c));
A=createA(n,g);
B=ones(size(A,1),1);
LB=-1*ones(1,length(pars));    %intercept
LB(length(c)+1)=0;             %theta
LB(length(c)+2)=1.001;         %omega
LB(length(c)+3:end)=0;         %garch parameters
UB=3*ones(1,length(pars));
UB(length(c)+1)=50;            %theta
UB(length(c)+2)=500;            %omega
UB(length(c)+3:end)=0.9999;    %garch parameters


% Reallocation in vech form of RC matrices to speed up estimation code
[mat]=cov_aggr(Vin);

tic
[estimates,loglik,exit]=...
    fmincon('TVLR_HEAVY_lik',pars,A,B,[],[],LB,UB,[],options,mat,Vin,DRtin,DRin,numlags,h,disp_lik);
TVLR_HEAVY.time=toc/60;


TVLR_HEAVY.stime    = estimates;
TVLR_HEAVY.int      = estimates(1:n*(n+1)/2); %use ivech to get it back
TVLR_HEAVY.longrun  = estimates((n*(n+1)/2)+1:(n*(n+1)/2)+2) ;
TVLR_HEAVY.GARCH    = [estimates((n*(n+1)/2)+3:(n*(n+1)/2)+2+n)' estimates(end-2-n+1:end-2)'];
TVLR_HEAVY.DCC      = [estimates(end-1) estimates(end)];
TVLR_HEAVY.lik      = -loglik;
TVLR_HEAVY.exit     = exit;

if exist('computeSTD','var') && computeSTD ==1
    [VCV]  =robustvcv('TVLR_HEAVY_lik',estimates,0,mat,Vin,DRtin,DRin,numlags,h);
    TVLR_HEAVY.std =diag(VCV).^(1/2);
    TVLR_HEAVY.tstat =estimates'./TVLR_HEAVY.std;
end


% Compute AIC BIC criterion for in-sample evaluation
npar =length(TVLR_HEAVY.stime);
[TVLR_HEAVY.aic ,TVLR_HEAVY.bic] = aic_bic(TVLR_HEAVY.lik,npar,tin-(numlags+1));

disp(' ')
disp('===================================================')
disp('         TIME VARYING LONG RUN HEAVY               ')
disp('===================================================')
disp(['Long run pars    : ',num2str(TVLR_HEAVY.longrun)])
disp(['GARCH parameters : ',num2str(TVLR_HEAVY.GARCH(1,:))])
disp(['                 : ',num2str(TVLR_HEAVY.GARCH(2,:))])
disp(['DCC parameters   : ',num2str(TVLR_HEAVY.DCC)])
disp(['Likelihood       : ',num2str(TVLR_HEAVY.lik)])
disp(['Estimation Time  : ',num2str(TVLR_HEAVY.time)])
disp('===================================================')

