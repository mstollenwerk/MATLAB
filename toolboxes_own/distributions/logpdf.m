function [y, varargout] = logpdf( dist, dta, varargin )

if nargout == 1
    if strcmp( dist, 'Wish' )
        y = matvWishlike(varargin{:}, dta, []);      
    elseif strcmp( dist, 'iWish' )
        y = matviWishlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'tWish' )
        y = matvtWishlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'itWish' )
        y = matvitWishlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'F' )
        y = matvFlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'Riesz' )
        y = matvRieszlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'iRiesz' )
        y = matviRieszlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'tRiesz' )
        y = matvtRieszlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'itRiesz' )
        y = matvitRieszlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'FRiesz' )
        y = matvFRieszlike(varargin{:}, dta, []);
    end
y = -y;
elseif nargout>=2
    if strcmp( dist, 'Wish' )
        [y, ~, score, hessian, param, fisherinfo] = matvWishlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'iWish' )
        [y, ~, score, hessian, param, fisherinfo] = matviWishlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'tWish' )
        [y, ~, score, hessian, param, fisherinfo] = matvtWishlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'itWish' )
        [y, ~, score, hessian, param, fisherinfo] = matvitWishlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'F' )
        [y, ~, score, hessian, param, fisherinfo] = matvFlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'Riesz' )
        [y, ~, score, hessian, param, fisherinfo] = matvRieszlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'iRiesz' )
        [y, ~, score, hessian, param, fisherinfo] = matviRieszlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'tRiesz' )
        [y, ~, score, hessian, param, fisherinfo] = matvtRieszlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'itRiesz' )
        [y, ~, score, hessian, param, fisherinfo] = matvitRieszlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'FRiesz' )
        [y, ~, score, hessian, param, fisherinfo] = matvFRieszlike(varargin{:}, dta, []);
    end
    y = -y;
    varargout{1} = score;
    varargout{2} = hessian;
    varargout{3} = param;
    varargout{4} = fisherinfo;
end

end

