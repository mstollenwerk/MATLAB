function [y, varargout] = logpdf( dist, dta, varargin )

if any(strcmp(dist, {'Wish', 'iWish', 'Riesz', 'iRiesz'})) && numel(varargin)==3
    C = varargin{3};
elseif any(strcmp(dist, {'tWish', 'itWish', 'F', 'tRiesz', 'itRiesz', 'FRiesz'})) && numel(varargin)==4
    C = varargin{4};
end

if exist('C', 'var')
    if nargout == 1
        if strcmp( dist, 'Wish' )
            y = matvWishlike(varargin{1:2}, dta, [], C);
        elseif strcmp( dist, 'iWish' )
            y = matviWishlike(varargin{1:2}, dta, [], C);
        elseif strcmp( dist, 'tWish' )
            y = matvtWishlike(varargin{1:3}, dta, [], C);
        elseif strcmp( dist, 'itWish' )
            y = matvitWishlike(varargin{1:3}, dta, [], C);
        elseif strcmp( dist, 'F' )
            y = matvFlike(varargin{1:3}, dta, [], C);
        elseif strcmp( dist, 'Riesz' )
            y = matvRieszlike(varargin{1:2}, dta, [], C);
        elseif strcmp( dist, 'iRiesz' )
            y = matviRieszlike(varargin{1:2}, dta, [], C);
        elseif strcmp( dist, 'tRiesz' )
            y = matvtRieszlike(varargin{1:3}, dta, [], C);
        elseif strcmp( dist, 'itRiesz' )
            y = matvitRieszlike(varargin{1:3}, dta, [], C);
        elseif strcmp( dist, 'FRiesz' )
            y = matvFRieszlike(varargin{1:3}, dta, [], C);
        end
    y = -y;
    elseif nargout>=2
        if strcmp( dist, 'Wish' )
            [y, ~, score, hessian, param, fisherinfo] = matvWishlike(varargin{1:2}, dta, [], C);
        elseif strcmp( dist, 'iWish' )
            [y, ~, score, hessian, param, fisherinfo] = matviWishlike(varargin{1:2}, dta, [], C);
        elseif strcmp( dist, 'tWish' )
            [y, ~, score, hessian, param, fisherinfo] = matvtWishlike(varargin{1:3}, dta, [], C);
        elseif strcmp( dist, 'itWish' )
            [y, ~, score, hessian, param, fisherinfo] = matvitWishlike(varargin{1:3}, dta, [], C);
        elseif strcmp( dist, 'F' )
            [y, ~, score, hessian, param, fisherinfo] = matvFlike(varargin{1:3}, dta, [], C);
        elseif strcmp( dist, 'Riesz' )
            [y, ~, score, hessian, param, fisherinfo] = matvRieszlike(varargin{1:2}, dta, [], C);
        elseif strcmp( dist, 'iRiesz' )
            [y, ~, score, hessian, param, fisherinfo] = matviRieszlike(varargin{1:2}, dta, [], C);
        elseif strcmp( dist, 'tRiesz' )
            [y, ~, score, hessian, param, fisherinfo] = matvtRieszlike(varargin{1:3}, dta, [], C);
        elseif strcmp( dist, 'itRiesz' )
            [y, ~, score, hessian, param, fisherinfo] = matvitRieszlike(varargin{1:3}, dta, [], C);
        elseif strcmp( dist, 'FRiesz' )
            [y, ~, score, hessian, param, fisherinfo] = matvFRieszlike(varargin{1:3}, dta, [], C);
        end
        y = -y;
        varargout{1} = score;
        varargout{2} = hessian;
        varargout{3} = param;
        varargout{4} = fisherinfo;
    end    
else
    if nargout == 1
        if strcmp( dist, 'Wish' )
            y = matvWishlike(varargin{1:2}, dta);
        elseif strcmp( dist, 'iWish' )
            y = matviWishlike(varargin{1:2}, dta);
        elseif strcmp( dist, 'tWish' )
            y = matvtWishlike(varargin{1:3}, dta);
        elseif strcmp( dist, 'itWish' )
            y = matvitWishlike(varargin{1:3}, dta);
        elseif strcmp( dist, 'F' )
            y = matvFlike(varargin{1:3}, dta);
        elseif strcmp( dist, 'Riesz' )
            y = matvRieszlike(varargin{1:2}, dta);
        elseif strcmp( dist, 'iRiesz' )
            y = matviRieszlike(varargin{1:2}, dta);
        elseif strcmp( dist, 'tRiesz' )
            y = matvtRieszlike(varargin{1:3}, dta);
        elseif strcmp( dist, 'itRiesz' )
            y = matvitRieszlike(varargin{1:3}, dta);
        elseif strcmp( dist, 'FRiesz' )
            y = matvFRieszlike(varargin{1:3}, dta);
        end
    y = -y;
    elseif nargout>=2
        if strcmp( dist, 'Wish' )
            [y, ~, score, hessian, param, fisherinfo] = matvWishlike(varargin{1:2}, dta);
        elseif strcmp( dist, 'iWish' )
            [y, ~, score, hessian, param, fisherinfo] = matviWishlike(varargin{1:2}, dta);
        elseif strcmp( dist, 'tWish' )
            [y, ~, score, hessian, param, fisherinfo] = matvtWishlike(varargin{1:3}, dta);
        elseif strcmp( dist, 'itWish' )
            [y, ~, score, hessian, param, fisherinfo] = matvitWishlike(varargin{1:3}, dta);
        elseif strcmp( dist, 'F' )
            [y, ~, score, hessian, param, fisherinfo] = matvFlike(varargin{1:3}, dta);
        elseif strcmp( dist, 'Riesz' )
            [y, ~, score, hessian, param, fisherinfo] = matvRieszlike(varargin{1:2}, dta);
        elseif strcmp( dist, 'iRiesz' )
            [y, ~, score, hessian, param, fisherinfo] = matviRieszlike(varargin{1:2}, dta);
        elseif strcmp( dist, 'tRiesz' )
            [y, ~, score, hessian, param, fisherinfo] = matvtRieszlike(varargin{1:3}, dta);
        elseif strcmp( dist, 'itRiesz' )
            [y, ~, score, hessian, param, fisherinfo] = matvitRieszlike(varargin{1:3}, dta);
        elseif strcmp( dist, 'FRiesz' )
            [y, ~, score, hessian, param, fisherinfo] = matvFRieszlike(varargin{1:3}, dta);
        end
        y = -y;
        varargout{1} = score;
        varargout{2} = hessian;
        varargout{3} = param;
        varargout{4} = fisherinfo;
    end
end

end

