function [Ht,St]=predict_TVLR_HEAVY_h(pars,mat,Vt,ret,numlags,h)
% This function computes direct forecasting at h steps-ahead of the
% TVLR-HEAVY model given the estimated parameters at each point in time.
%==========================================================================
%INPUT
%  pars    = vector of initial parameter values
%  mat     = Txn*(n+1)/2 matrix of realized series in vech form
%  Vt      = series of realized covariance matrices
%  ret     = Txn matrix of daily returns
%  numlags = number of days backword in the long-term component (i.e. 260)
%  h       = time horizon for direct forecasting
%==========================================================================
% Author : Manuela Braione, CORE-UCL 2016
%==========================================================================

[~, n, T]=size(Vt);
N        =n*(n+1)/2;

%Long run parameters
intercept=pars(1:(n*(n+1)/2))';
Lambdabar=ivech(intercept)*ivech(intercept)';
teta     =pars((n*(n+1)/2)+1);
omega    =pars((n*(n+1)/2)+2);
Omega    =repmat(beta_weights_cres(numlags,1,omega),1,N);

%Short run univariate garch parameters
gamma =pars((n*(n+1)/2)+3:((n*(n+1)/2)+3)+n-1);
delta =pars((((n*(n+1)/2)+3)+n):end-2);
univ  =[gamma; delta]';
const =1-sum(univ,2);

%Short run correlation parameters
alpha=pars(end-1);
beta=pars(end);

% Initialize variables
epsilont=zeros(T,n);
DD=zeros(n,n,T);
St=zeros(n,n,T);
Ht=zeros(n,n,T);
Ginv=zeros(n,n,T);
Gstore=zeros(n,n,T);
RC=zeros(n,n,T);
I=eye(n);
v=zeros(T,n);
Hstar=zeros(T,n);
Qtstar=zeros(n,n,T);

% Reconstruct process up to time T and predict in T+h
for t=(numlags+h):T
    COV=mat(t-numlags-h+1:t-h,:);
    Q=sum(Omega.*COV);
    Vstore=vecTomat(Q);
    St(:,:,t)= Lambdabar + teta*Vstore;
    G=chol(St(:,:,t))';
    Gstore(:,:,t)=G;
    Ginv(:,:,t)=G^(-1);
    RC(:,:,t)=Ginv(:,:,t)*Vt(:,:,t)*Ginv(:,:,t)';
end

for i=1:n
    v(:,i)=reshape(RC(i,i,:),T,1);
end

Hstar(numlags+1,:)=const;
Qtstar(:,:,numlags+h+1)=I;

for t=numlags+h+1:T
    Hstar(t,:)=const' + (univ(:,1)'.*v(t-h,:)) + (univ(:,2)'.*Hstar(t-h,:));
    DD(:,:,t)=(diag(Hstar(t,:).^0.5));
    invD=DD(:,:,t)^(-1);
    epsilont(t,:)=(invD*Ginv(:,:,t))*ret(t,:)';
end

for t=numlags+h+1:T+h
    Hstar(t,:)=const' + (univ(:,1)'.*v(t-h,:)) + (univ(:,2)'.*Hstar(t-h,:));
    DD(:,:,t)=(diag(Hstar(t,:).^0.5));
    Dt=DD(:,:,t);
    Qtstar(:,:,t)=(1-alpha-beta)*I +alpha*(epsilont(t-h,:)'*epsilont(t-h,:))+(beta*Qtstar(:,:,t-h));
    Rdcc=diag(diag(Qtstar(:,:,t)).^-.5)*Qtstar(:,:,t)*diag(diag(Qtstar(:,:,t)).^-.5);
    Rtstar=Dt*Rdcc*Dt;
    
    if t>T
        COV=mat(t-numlags-h+1:t-h,:);
        Q=sum(Omega.*COV);
        Vstore=vecTomat(Q);
        St(:,:,t)= Lambdabar + teta*Vstore;
        G=chol(St(:,:,t))';
        Ht(:,:,t)=G*Rtstar*G';
    else
        Ht(:,:,t)=Gstore(:,:,t)*Rtstar*Gstore(:,:,t)';
        
    end
end

end

function [be_values]=beta_weights_cres(dayLag,k1,k2)

%Computes the weights in the Beta function. Weights are INCREASING to be
%compatible with the way in which "mat" is constructed.

u=linspace(1-eps,eps,dayLag)';
be_values=u.^(k1-1).*(1-u).^(k2-1);  
be_values=be_values/sum(be_values);
end

