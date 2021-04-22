function [ Y ] = nw_kernel( X, bandwidth )
%NW_KERNEL extracts long-term trend w Nadaraya-Watson kernel estimator
%   Detailed explanation goes here

[n,~,T] = size(X);

Xvech = zeros(n*(n+1)/2,T);
for i=1:T
    Xvech(:,i) = vech(X(:,:,i));
end

kern = normpdf((repmat(1:T,T,1)-repmat(1:T,T,1)')/T/bandwidth);
kernsum = sum(kern,2);

Y = zeros(n,n,T);
for i=1:T
    Y(:,:,i) = ivech(Xvech*kern(i,:)'./kernsum(i),'sym');
end

end

