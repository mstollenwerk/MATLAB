function [y] = loglpwdet3d(X,powers,varargin)
%LPWDET Lower Power Weighted Determinant
%   Detailed explanation goes here

if nargin == 3
    dgC = varargin{:};
end

T = max(size(powers,1),size(dgC,1));
y = NaN(T,1);
for ii = 1:T
    if nargin == 3
        if size(powers,1) > 1 && size(dgC,1) > 1
            y(ii) = loglpwdet([],powers(ii,:),dgC(ii,:)');
        elseif size(powers,1) == 1 && size(dgC,1) > 1
            y(ii) = loglpwdet([],powers,dgC(ii,:)');   
        elseif size(powers,1) > 1 && size(dgC,1) == 1
            y(ii) = loglpwdet([],powers(ii,:),dgC');
        else
            error('Incompatible input sizes.')
        end
    else
        if size(powers,1) > 1 && size(X,1) > 1
            y(ii) = loglpwdet(X(:,:,ii),powers(ii,:));
        elseif size(powers,1) == 1 && size(X,1) > 1
            y(ii) = loglpwdet(X(:,:,ii),powers);   
        elseif size(powers,1) > 1 && size(X,1) == 1
            y(ii) = loglpwdet(X,powers(ii,:));
        else
            error('Incompatible input sizes.')
        end        
    end
end


end

