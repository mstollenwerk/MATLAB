function [llik,llog]=TVLR_HEAVY_lik(pars,mat,Cin,DRtin,ret,numlags,h,disp_out)
%
% This function estimates IN ONE STEP the parameters of the TVLT-HEAVY-P (equation)
% assuming a time-varying long run specification.
%==========================================================================
%INPUT
%  pars    (numeric) vector of initial values for the model parameters
%  C       (numeric array) series of realized covariance matrices
%  DRt     (numeric array) series of daily covariance matrices
%  ret     (numeric matrix) matrix of daily returns
%  numlags (numeric scalar) aggregation period (usually 260 days)
%  h       (numeric scalar) time horizon for direct forecasting
%==========================================================================
%Manuela Braione, CORE 2016

if nargin<8
    disp_out=0;
end

% Initialization and variables declaration
[~, n, T]=size(Cin);
LLstore=zeros(n,n,T);
RC=zeros(n,n,T);
I=eye(n);
llog=zeros(T,1);
zt=zeros(T,n);
DD=zeros(n,n,T);


%Distributed lag component
intercept=pars(1:(n*(n+1)/2))';
intL=ivech(intercept)*ivech(intercept)';
teta=pars((n*(n+1)/2)+1);
k2=pars((n*(n+1)/2)+2);
k1=1;
Omega=repmat(beta_weights_cres(numlags,k1,k2),1,n*(n+1)/2);


%univariate garch
alpha=pars((n*(n+1)/2)+3:((n*(n+1)/2)+3)+n-1);
beta=pars((((n*(n+1)/2)+3)+n):end-2);
univ=[alpha; beta]';
const=1-sum(univ,2);

%correlation parameters
aa=pars(end-1);
bb=pars(end);



for t=(numlags+h):T
    COV=mat(t-numlags-h+1:t-h,:);
    Q=sum(Omega.*COV);
    Qt=vecTomat(Q);
    Mt= intL + teta*Qt;
    LL=chol(Mt)';
    LLstore(:,:,t)=LL;
    Linv(:,:,t)=LL^(-1);
    RC(:,:,t)=Linv(:,:,t)*Cin(:,:,t)*Linv(:,:,t)';
end


for i=1:n
    v(:,i)=reshape(RC(i,i,:),T,1);
end
S(numlags+h,:)=const;

for t=numlags+h+1:T
    S(t,:)=const' + (univ(:,1)'.*v(t-h,:)) + (univ(:,2)'.*S(t-h,:));
    DD(:,:,t)=(diag(S(t,:).^0.5));
    invD=DD(:,:,t)^(-1);
    invL=Linv(:,:,t);
    zt(t,:)=(invD*invL)*ret(t,:)';
end

Rdcc=I;

for t=(numlags+h+2):T
    DDR=DRtin(:,:,t);
    l=LLstore(:,:,t);
    Dt=DD(:,:,t);
    Rdcc=(1-aa-bb)*I +aa*(zt(t-h,:)'*zt(t-h,:))+(bb*Rdcc);
    Rdcc=diag(diag(Rdcc).^-.5)*Rdcc*diag(diag(Rdcc).^-.5);
    S=Dt*Rdcc*Dt;
    LSL=l*S*l';
    llog(t)=.5*(log(det(LSL))+ trace((LSL)^(-1)*DDR));
end

llik=sum(llog(numlags+h+2:end));

if disp_out==1
    fprintf(1,['\nSum log likelihood -> ', num2str(llik)]);
end
end

function [be_values]=beta_weights_cres(dayLag,k1,k2)
%Beta weights function: weights are INCREASING over time.

u=linspace(1-eps,eps,dayLag)';
be_values=u.^(k1-1).*(1-u).^(k2-1);
be_values=be_values/sum(be_values);
end