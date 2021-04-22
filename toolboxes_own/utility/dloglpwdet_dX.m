function [y] = dloglpwdet_dX(X,powers)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

p = size(X,1);

if size(powers,1) < size(powers,2)
    powers = powers';
end

%%
L = chol(X,'lower');
y = L'\diag(powers)/L;
% %%
% % This version is slower but equivalent.
% % y = powers(p)*inv(X);
% % for ii = 1:(p-1)
% %     y = y + (powers(ii) - powers(ii+1)) * [ inv(X(1:ii,1:ii))   zeros(ii,p-ii);
% %                                             zeros(p-ii,ii)      zeros(p-ii,p-ii) ];
% % end

end

%% DEPRECATED:
% p_ = p*(p+1)/2;
% L = chol(X,'lower');
% y = zeros(p_,1);
% for ii = 1:p
%     dcholXii_dvechSigma = NaN(p_,1);
%     jj = 1;
%     for l = 1:p 
%         for k = l:p
%             dcholXii_dvechSigma(jj) = dcholXij_dXkl(X,ii,ii,k,l);
%             jj = jj + 1;
%         end
%     end
%     y = y + 2*powers(ii)/L(ii,ii)*dcholXii_dvechSigma;
% end
%%
% p_ = p*(p+1)/2;
% L = chol(X,'lower');
% dcholX_dX_ = dcholX_dX(X);
% y = dcholX_dX_'*vech(diag(2*powers./diag(L)));