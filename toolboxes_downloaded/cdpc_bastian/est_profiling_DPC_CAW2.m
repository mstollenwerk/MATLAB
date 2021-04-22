function est_profiling_DPC_CAW2

% dbstop if error
% dbstop if warning
% warning off

load('Scalar_CAW_est','C','R_adj','R');

T = size(R_adj,3);

k = 5;

options = optimset('HessUpdate','bfgs','Display','off','MaxFunEvals',50000,'MaxIter',50000,'TolFun',1e-5,...
                    'TolX',1e-5); % use medium scale algorithm

%%%%%%   construct a grid for a and b (parameters of eigenvalue recursion)

scal = 0.001;
a_ = scal:scal:0.5;
b_ = 0.5:scal:1-scal;

le_ = length(a_); % entspicht length(b_)

I = eye(k);

x0 = [0.01 0.01];

%--------------------------------------------------------------------------

x1_  = zeros(le_,le_,k);
x2_  = zeros(le_,le_,k);

FVAL = zeros(le_,le_);

zz = 0;
for i=1:le_
     for j=1:le_
         zz = zz+1;
         disp(zz/(le_*le_))

         a = abs(a_(i));
         b = abs(b_(j));

         if a + b < 1
             Q=zeros(k,k,T);
             Q(:,:,1) = eye(k);

             for t=2:T
                 Q(:,:,t)= (1-a-b)*I + a*R_adj(:,:,t-1) + b*Q(:,:,t-1);
                 [L(:,:,t), D] = eig(Q(:,:,t));
                 g(:,t) = diag( L(:,:,t)'*R_adj(:,:,t)*L(:,:,t) );
             end

             fval_ = 0;
             for indi=1:k

                     [x,fval,exitflag] = fminunc(@eigenvalue,x0,options,g,T,indi);
                     x1_(i,j,indi) = x(1);
                     x2_(i,j,indi) = x(2);

                     fval_ = fval_ + fval;
             end
             FVAL(i,j) = fval_;
         else
             FVAL(i,j) = realmax;
         end

     end
end
%--------------------------------------------------------------------------

save('est_profiling_DPC_CAW2')
% load('est_profiling_DPC_CAW2')

%--------------------------------------------------------------------------
% find (i,j)-Kombination im Optimum (minimales FVAL:
%--------------------------------------------------------------------------

min_ = min(min(FVAL));
finder = find(FVAL==min_);

if mod(finder,le_) == 0
     Zeile = le_;
     Spalte = (finder/le_);
else
     Zeile  = mod(finder,le_);
     Spalte = floor(finder/le_)+1;
end

% berechne Parameter-Schätzer:

a  = a_(Zeile);
b  = b_(Spalte);

al = abs (squeeze(x1_(Zeile,Spalte,:))); be = abs (squeeze(x2_(Zeile,Spalte,:)));

%--------------------------------------------------------------------------
% Compute in-sample forecasts:
%--------------------------------------------------------------------------

Q=zeros(k,k,T);
Q(:,:,1) = eye(k);
d = zeros(k,T);
d(:,1) = ones(k,1);
g = zeros(k,T);
g(:,1) = ones(k,1);

for t=2:T
     Q(:,:,t)= (1-a-b)*I + a*R_adj(:,:,t-1) + b*Q(:,:,t-1);
     [L(:,:,t), D] = eig(Q(:,:,t));
     g(:,t) = diag( L(:,:,t)'*R_adj(:,:,t)*L(:,:,t) );

     for i=1:k
         d(i,t) = (1-al(i)-be(i)) + al(i)*g(i,t-1) + be(i)*d(i,t-1);
     end
     S(:,:,t) = L(:,:,t)*diag(d(:,t))*L(:,:,t)'; end

figure(1)
subplot(5,3,1), plot(squeeze(R_adj(1,1,:))), hold on,
plot(squeeze(S(1,1,:)),'r')
subplot(5,3,2), plot(squeeze(R_adj(2,1,:))), hold on,
plot(squeeze(S(2,1,:)),'r')
subplot(5,3,3), plot(squeeze(R_adj(3,1,:))), hold on,
plot(squeeze(S(3,1,:)),'r')
subplot(5,3,4), plot(squeeze(R_adj(4,1,:))), hold on,
plot(squeeze(S(4,1,:)),'r')
subplot(5,3,5), plot(squeeze(R_adj(5,1,:))), hold on,
plot(squeeze(S(5,1,:)),'r')
subplot(5,3,6), plot(squeeze(R_adj(2,2,:))), hold on,
plot(squeeze(S(2,2,:)),'r')
subplot(5,3,7), plot(squeeze(R_adj(3,2,:))), hold on,
plot(squeeze(S(3,2,:)),'r')
subplot(5,3,8), plot(squeeze(R_adj(4,2,:))), hold on,
plot(squeeze(S(4,2,:)),'r')
subplot(5,3,9), plot(squeeze(R_adj(5,2,:))), hold on,
plot(squeeze(S(5,2,:)),'r')
subplot(5,3,10), plot(squeeze(R_adj(3,3,:))), hold on,
plot(squeeze(S(3,3,:)),'r')
subplot(5,3,11), plot(squeeze(R_adj(4,3,:))), hold on,
plot(squeeze(S(4,3,:)),'r')
subplot(5,3,12), plot(squeeze(R_adj(5,3,:))), hold on,
plot(squeeze(S(5,3,:)),'r')
subplot(5,3,13), plot(squeeze(R_adj(4,4,:))), hold on,
plot(squeeze(S(4,4,:)),'r')
subplot(5,3,14), plot(squeeze(R_adj(5,4,:))), hold on,
plot(squeeze(S(5,4,:)),'r')
subplot(5,3,15), plot(squeeze(R_adj(5,5,:))), hold on,
plot(squeeze(S(5,5,:)),'r')



end



function llf = eigenvalue(para,g,T,indi)


al = abs(para(1));
be = abs(para(2));

d = 1;
llf = 0;
if al + be < 1

     for t=2:T

         d = (1-al-be) + al*g(indi,t-1) + be*d;

         llf = llf + ( -0.5*( log(d) + g(indi,t)/d ) );
     end

     llf = -sum(llf);

     if isnan(llf) || isinf(llf) || isreal(llf) == 0
         llf=1e6;
     end

else
     llf=1e6;
end


end
