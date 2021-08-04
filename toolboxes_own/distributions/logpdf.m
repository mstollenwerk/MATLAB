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
elseif nargout==2
    if strcmp( dist, 'Wish' )
        [y, ~, score] = matvWishlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'iWish' )
        [y, ~, score] = matviWishlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'tWish' )
        [y, ~, score] = matvtWishlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'itWish' )
        [y, ~, score] = matvitWishlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'F' )
        [y, ~, score] = matvFlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'Riesz' )
        [y, ~, score] = matvRieszlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'iRiesz' )
        [y, ~, score] = matviRieszlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'tRiesz' )
        [y, ~, score] = matvtRieszlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'itRiesz' )
        [y, ~, score] = matvitRieszlike(varargin{:}, dta, []);
    elseif strcmp( dist, 'FRiesz' )
        [y, ~, score] = matvFRieszlike(varargin{:}, dta, []);
    end
y = -y;
varargout{1} = score;
else
    error('Number of output arguments can only be one or two.') 
end

end

