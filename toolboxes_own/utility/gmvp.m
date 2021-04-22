function [ weights, var ] = gmvp( covmat )
%GMVP weights+variance of global min variance portf given covariance matrix
%   DEMEANED RETURNS ONLY!
%% Initializations
[n,~,T]=size(covmat);                                                      
weights=NaN(T,n); 
var=NaN(T,1);                   
%%
for t=1:T
    inv_=covmat(:,:,t)\ones(n,1);
    weights(t,:)=inv_/(ones(1,n)*inv_);
    var(t)=weights(t,:)*covmat(:,:,t)*weights(t,:)';
end

end
