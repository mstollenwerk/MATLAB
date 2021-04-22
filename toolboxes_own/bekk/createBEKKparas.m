function [ paras ] = createBEKKparas( n, p, q )
%createBEKKparas creates parameters for the BEKK(p,q) recursion with 
%existing finite mean
%
%% Inputs/Outputs
%
%%%%%%%%
%Inputs%
%%%%%%%%
%
%          n - Dimension of recursion
%          p - Order of autoregressive components 
%          q - Number of data lags
%
%%%%%%%%%
%Outputs%
%%%%%%%%%
%
%      paras - parameters for BEKK(p,q) recursion:
%              [vech(C); vec(B_1); ...; vec(B_p); vec(A_1); ...;vec(A_q)]
%
%% 
%%

M=CAWmean([2*ones(n*(n+1)/2+(p+q)*n^2,1);21],p,q);

while M==0
    A=-eye(n);
    B=-eye(n);

    while A(1,1)<0||B(1,1)<0
        A=randn(n)/2+.25;
        B=randn(n)/2+.25;
    end
    %%
    ARpara_p = (rand(1)+rand(1))/2;
    ARpara_q = (rand(1)+rand(1))/2;

    wp = ARpara_p.^(1:p);
    wp = wp/sum(wp);
    wq = ARpara_q.^(1:q);
    wq = wq/sum(wq);

    matrices=zeros(n,n,p+q);
    for i=1:p
        matrices(:,:,i) = wp(i)*A;
    end
    for j=1:q
        matrices(:,:,p+j) = wq(j)*B;
    end

    % Arbitrary constant matrix C
    C = mean(matrices,3)*mean(matrices,3)';

    % Parameters [vech(C); vec(B_1); ...; vec(B_p); vec(A_1); ...;vec(A_q)]
    paras = [vech(C);reshape(matrices,[],1)];
    
    M=CAWmean([paras;21],p,q);
end
    end