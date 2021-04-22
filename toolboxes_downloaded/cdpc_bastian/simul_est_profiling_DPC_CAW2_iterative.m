%%
function simul_est_profiling_DPC_CAW2_iterative(R,i)

% dbstop if error
% dbstop if warning
% warning off

R_adj = R;

T = size(R_adj,3);

k = size(R_adj,2);

x0_eigenvalue = [0.01 0.01];

x0_eigenvector = [0.1 0.1];

options = optimset('HessUpdate','bfgs','MaxFunEvals',50000,'MaxIter',50000,'TolFun',1e-5,...
                   'TolX',1e-5); % use medium scale algorithm

%options_fminsearch = optimset('MaxIter',20,'Display','iter'); % use medium scale algorithm
options_patternsearch = optimset('MaxIter',20,'Display','iter');

abbruch = 0;

a = abs(x0_eigenvector(1));
b = abs(x0_eigenvector(2));

% This is initialization
Q=zeros(k,k,T);
Q(:,:,1) = eye(k);
I = eye(k);
for t=2:T
     Q(:,:,t)= (1-a-b)*I + a*R_adj(:,:,t-1) + b*Q(:,:,t-1);
     [L(:,:,t), D] = eig(Q(:,:,t));
     g(:,t) = diag( L(:,:,t)'*R_adj(:,:,t)*L(:,:,t) );
end

llf_end_ = 0;
checker  = 1;

while abbruch == 0

     for indi=1:k

         [x,fval,exitflag] = fminunc(@eigenvalue,x0_eigenvalue,options,g,T,indi);

         x0_eigenvalue = x;

         al = abs(x(1));
         be = abs(x(2));

         d(1,indi) = 1;

         for t=2:T
             d(t,indi) = (1-al-be) + al*g(indi,t-1) + be*d(t-1,indi);
         end

         al_(indi) = al;
         be_(indi) = be;
     end


     disp(checker)

%--------------------------------------------------------------------------


     %x = patternsearch(@(x) eigenvector(x,R_adj,d,T,k),x0_eigenvector,[],[],[],[],[],[],[],options_patternsearch);
     x = fminunc(@eigenvector,x0_eigenvector,options,R_adj,d,T,k);   
     
     x0_eigenvector = x;

     a = abs(x(1));
     b = abs(x(2));

     Q=zeros(k,k,T);
     Q(:,:,1) = eye(k);
     I = eye(k);
     for t=2:T
         Q(:,:,t)= (1-a-b)*I + a*R_adj(:,:,t-1) + b*Q(:,:,t-1);
         [L(:,:,t), D] = eig(Q(:,:,t));
         g(:,t) = diag( L(:,:,t)'*R_adj(:,:,t)*L(:,:,t) );
     end


%--------------------------------------------------------------------------
     % compute overall likelihood

     llf = 0;
     for t=2:T

         llf = llf -0.5*(  log(det(diag(d(t,:))))  +  trace( diag(1./d(t,:))*L(:,:,t)'*R_adj(:,:,t)*L(:,:,t) )   );

     end

     llf_end = llf;

     checker = abs(  (llf_end - llf_end_)/llf_end_  );

     llf_end_ = llf;


     if checker < 0.000001
        abbruch = 1;
     end

end


%--------------------------------------------------------------------------

save(strcat('simul_est_profiling_DPC_CAW2_iterative_',num2str(i)))
% load('simul_est_profiling_DPC_CAW2_iterative')

end


%%
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


%%
function llf = eigenvector(para,R_adj,d,T,k)


a = abs(para(1));
b = abs(para(2));q

llf = 0;

Q = eye(k);
I = eye(k);

if a + b < 1

     for t=2:T
         Q = (1-a-b)*I + a*R_adj(:,:,t-1) + b*Q;
         [L(:,:,t), D] = eig(Q);

         S = L(:,:,t)*diag(d(t,:))*L(:,:,t)';

         llf = llf + (   -0.5*( log(det(S)) + trace(inv(S)*R_adj(:,:,t)) )   );

     end

     llf = -sum(llf);

     if isnan(llf) || isinf(llf) || isreal(llf) == 0
         llf=1e6;
     end

else
     llf=1e6;
end


end
