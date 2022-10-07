function y = diag3d(X)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% if numel(size(X))==3
    [p,~,N] = size(X);

    y = NaN(N,p);
    for ii = 1:N
        y(ii,:) = diag(X(:,:,ii));
    end
% elseif numel(size(X))==2
%     [N,p] = size(X);
% 
%     y = NaN(p,p,N);
%     for ii = 1:N
%         y(:,:,ii) = diag(X(ii,:));
%     end    
% end

end

