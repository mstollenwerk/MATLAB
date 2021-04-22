function [EST,FORE]=est_dccgarch_fore(r,STARTER,fore_)
% fore_     : forecasting horizon
% r         : T X n -dimensional return vector
% STARTER   : starting values

warning off;

randn('seed',123);
rand('seed',123);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zz=size(r);
n=zz(2); % series
T=zz(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S=cov(r);


x0=STARTER;

ub = 1*ones(1,3*n+2);
lb = 0.000001*ones(1,3*n+2);


param_trans=zeros(3*n+2,1);
trans=sqrt((x0-lb)./(ub-lb));
param_trans=atan( trans./sqrt(1-(trans.^2)) );


% Optimierung
%------------

% BFGS
%-----

options = optimset('LargeScale','off','HessUpdate','bfgs','Display','iter','MaxFunEvals',5000000,'MaxIter',5000000,'TolFun',1e-4,...
                  'TolX',1e-4); % use medium scale algorithm

              
% First step of two-step estimation procedure              
%--------------------------------------------
llf_ = zeros(n,1);
estimate_ = zeros(n,3);
eps=zeros(T,n);
h_=zeros(T+fore_,n);
for i=1:n 
    param_trans_=[param_trans(i),param_trans(n+i),param_trans(2*n+i)];
    
    [x,fval] = fminunc(@negloglikeli_1,param_trans_,options,T,r,lb,ub,i);
    
    estimate(1)=lb(i)+(ub(i)-lb(i)).*(sin(x(1)).^2);
    estimate(2)=lb(n+i)+(ub(n+i)-lb(n+i)).*(sin(x(2)).^2);
    estimate(3)=lb(2*n+i)+(ub(2*n+i)-lb(2*n+i)).*(sin(x(3)).^2);
    estimate_(i,:)=estimate;
    
    llf_(i)=-fval;
    
    omega=estimate(1);
    kappa=estimate(2);
    lambda=estimate(3);

    h=zeros(1,T+fore_);
    h(1)=var(r(:,i)); 
    % do forecasting:
    for t=2:T+fore_
        if t<= T+1
            h(t)=omega+kappa*(r(t-1,i)^2)+lambda*h(t-1);        
        elseif t>T+1 % multiperiod forecasting
            h(t)=omega+kappa*h(t-1)+lambda*h(t-1);        
        end
    end    
    eps(:,i)=r(:,i)./sqrt(h(1:T)');
    h_(:,i)=h';
end


    
% Second step of two-step estimation procedure              
%--------------------------------------------

param_trans_=[param_trans(3*n+1),param_trans(3*n+2)];
[x,fval] = fminunc(@negloglikeli_2,param_trans_,options,T,r,eps,S,lb,ub);    
        
estimate(1)=lb(3*n+1)+(ub(3*n+1)-lb(3*n+1)).*(sin(x(1)).^2);
estimate(2)=lb(3*n+2)+(ub(3*n+2)-lb(3*n+2)).*(sin(x(2)).^2);

a=estimate(1);
b=estimate(2);
estimate_2(1,:)=estimate(1:2);

Q=zeros(n,n,T+fore_);
R=zeros(n,n,T+fore_);
Q(:,:,1)=S; % Set to unconditional MGARCH mean
for t=2:T+fore_
    if t<= T+1
        Q(:,:,t)=S.*(1-a-b)+a*(eps(t-1,:)'*eps(t-1,:))+b.*Q(:,:,t-1); % Covariance matrix of eps
        R(:,:,t)= diag((1./( sqrt( diag(Q(:,:,t)) ) ))) *Q(:,:,t)* diag((1./( sqrt( diag(Q(:,:,t)) ) ))); % Korrelationsmatrix! Einsen auf der Hauptdiagonalen        
    elseif t>T+1 % multiperiod forecasting
        Q(:,:,t)=S.*(1-a-b)+a*Q(:,:,t-1)+b.*Q(:,:,t-1); % Covariance matrix of eps
        R(:,:,t)= diag((1./( sqrt( diag(Q(:,:,t)) ) ))) *Q(:,:,t)* diag((1./( sqrt( diag(Q(:,:,t)) ) ))); % Korrelationsmatrix! Einsen auf der Hauptdiagonalen        
    end
end
llf_(3)=-fval;


% Total log-likelihood:
%---------------------
llf=sum(llf_);

% Cova-Matrix aufstellen
%-----------------------

R(:,:,1)=corrcoef(r);
for t=1:T+fore_
    D(:,:,t)=diag(sqrt(h_(t,:)));
    H(:,:,t)=D(:,:,t)*R(:,:,t)*D(:,:,t);
end

FORE=H(:,:,T+1:T+fore_);

EST=[estimate_(:,1)', estimate_(:,2)', estimate_(:,3)', estimate_2];

%save([mfilename,'_',int2str(number_)]);



%-------------------------------------------------------------------------


function [llf_1]=negloglikeli_1(para,T,r,lb,ub,i) % individual GARCH estimation; 'i' selects the time-series
 
pa=zeros(3,1);
pa(1)=lb(i)+(ub(i)-lb(i)).*(sin(para(1)).^2);
pa(2)=lb(5+i)+(ub(5+i)-lb(5+i)).*(sin(para(2)).^2);
pa(3)=lb(10+i)+(ub(10+i)-lb(10+i)).*(sin(para(3)).^2);

omega=pa(1);
kappa=pa(2);
lambda=pa(3);

%check for stationarity
if  kappa+lambda>=1
    llf_1=1000000;
    return
end

    
h=var(r(:,i));
for t=2:T
    h(t)=omega+kappa*(r(t-1,i)^2)+lambda*h(t-1);
    llf(t-1)=-0.5*log(2*pi)-0.5*log(h(t))-0.5*(r(t,i)^2)/h(t);
end

llf_1=-sum(llf);




function [llf_2]=negloglikeli_2(para,T,r,eps,S,lb,ub) % individual GARCH estimation; 'i' selects the time-series

pa=zeros(2,1);
pa(1)=lb(16)+(ub(16)-lb(16)).*(sin(para(1)).^2);
pa(2)=lb(17)+(ub(17)-lb(17)).*(sin(para(2)).^2);
pa';

a=pa(1);
b=pa(2);

%check for stationarity
if a+b>=1 
    llf_2=1000000;
    return
end


Q=S; % Set to unconditional MGARCH mean
for t=2:T
    Q(:,:,t)=S.*(1-a-b)+a*(eps(t-1,:)'*eps(t-1,:))+b.*Q(:,:,t-1); % Covariance matrix of eps
    R(:,:,t-1)= diag((1./( sqrt( diag(Q(:,:,t)) ) ))) *Q(:,:,t)* diag((1./( sqrt( diag(Q(:,:,t)) ) ))); % Korrelationsmatrix! Einsen auf der Hauptdiagonalen
    llf(t-1)=-0.5*(log(det(R(:,:,t-1)))+eps(t,:)*inv(R(:,:,t-1))*eps(t,:)'-eps(t,:)*eps(t,:)');
end

llf_2=-sum(llf);

