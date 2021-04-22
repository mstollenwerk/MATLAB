function [ var ] = portfvar( weights, realCovMat )
%PORTFVAR returns portfolio variance given portf. weights and cov matrices

if size(weights,1)>size(weights,2)
    weights = weights';
    warning('weights input has been transposed');
end
if size(weights,1)~=size(realCovMat,1)
    error('Same cross-sectional size of weights and realCovMat is required')
elseif size(weights,2)~=size(realCovMat,3)
    error('Same number of portfolios in weights and realCovMat is required')
end
[~,T]=size(weights);                                                      
%%
var=NaN(T,1);
for t=1:T
    var(t)=weights(:,t)'*realCovMat(:,:,t)*weights(:,t);
end

end

