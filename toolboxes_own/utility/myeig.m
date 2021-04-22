function [ V, D, I ] = myeig( X, varargin )
%Find eigenvector/value and enforce a sign convention 
%   by making the largest eigenvector element non-negative. Also return
%   the corresponding indices.

optargs = {'Biggest' 'Vector'};
optargs(1:length(varargin)) = varargin;
[type, option] = optargs{:};

[V, D]=eig(X,option);

if strcmp(type,'First')
    V = bsxfun(@times,V,sign(V(1,:)));
elseif strcmp(type,'Biggest')
    [~,I] = max(abs(V),[],1);
    I = sub2ind(size(V),I,1:size(V,2));
    V = bsxfun(@times,V,sign(V(I)));
end

end