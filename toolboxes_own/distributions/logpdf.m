function [y, varargout] = logpdf( dist, scaling, dta, varargin )

if nargout == 1
    if strcmp( dist, 'Wish' )
        y = matvWishlike(varargin{:}, dta, [], scaling);      
    elseif strcmp( dist, 'iWish' )
        y = matviWishlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'tWish' )
        y = matvtWishlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'itWish' )
        y = matvitWishlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'F' )
        y = matvFlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'Riesz' )
        y = matvRieszlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'iRiesz' )
        y = matviRieszlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'tRiesz' )
        y = matvtRieszlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'itRiesz' )
        y = matvitRieszlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'FRiesz' )
        y = matvFRieszlike(varargin{:}, dta, [], scaling);
    end
y = -y;
elseif nargout==2
    if strcmp( dist, 'Wish' )
        [y, ~, score] = matvWishlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'iWish' )
        [y, ~, score] = matviWishlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'tWish' )
        [y, ~, score] = matvtWishlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'itWish' )
        [y, ~, score] = matvitWishlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'F' )
        [y, ~, score] = matvFlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'Riesz' )
        [y, ~, score] = matvRieszlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'iRiesz' )
        [y, ~, score] = matviRieszlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'tRiesz' )
        [y, ~, score] = matvtRieszlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'itRiesz' )
        [y, ~, score] = matvitRieszlike(varargin{:}, dta, [], scaling);
    elseif strcmp( dist, 'FRiesz' )
        [y, ~, score] = matvFRieszlike(varargin{:}, dta, [], scaling);
    end
y = -y;
varargout{1} = score;
else
    error('Number of output arguments can only be one or two.') 
end

end

