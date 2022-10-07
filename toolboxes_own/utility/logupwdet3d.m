function [y] = logupwdet3d(X,powers,varargin)
%LPWDET Lower Power Weighted Determinant
%   Detailed explanation goes here

if nargin == 3
    dgCholU = varargin{:};
end

T = max(size(powers,1),size(dgCholU,1));
y = NaN(T,1);
for ii = 1:T
    if nargin == 3
        if size(powers,1) > 1 && size(dgCholU,1) > 1
            y(ii) = logupwdet([],powers(ii,:),dgCholU(ii,:)');
        elseif size(powers,1) == 1 && size(dgCholU,1) > 1
            y(ii) = logupwdet([],powers,dgCholU(ii,:)');   
        elseif size(powers,1) > 1 && size(dgCholU,1) == 1
            y(ii) = logupwdet([],powers(ii,:),dgCholU');
        else
            error('Incompatible input sizes.')
        end
    else
        if size(powers,1) > 1 && size(X,1) > 1
            y(ii) = logupwdet(X(:,:,ii),powers(ii,:));
        elseif size(powers,1) == 1 && size(X,1) > 1
            y(ii) = logupwdet(X(:,:,ii),powers);   
        elseif size(powers,1) > 1 && size(X,1) == 1
            y(ii) = logupwdet(X,powers(ii,:));
        else
            error('Incompatible input sizes.')
        end        
    end
end


end

