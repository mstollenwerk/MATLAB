function [ eparam, tstats, logL, optimoutput ] = matvest(X,dist,x0,varargin)
%MATVEST Simple wrapper around my matrix-variate distribution estimations

if strcmp( dist, 'Wish' )
    [ eparam, tstats, logL, optimoutput ] = matvWishest(X,x0,varargin{:});
elseif strcmp( dist, 'iWish' )
    [ eparam, tstats, logL, optimoutput ] = matviWishest(X,x0,varargin{:});
elseif strcmp( dist, 'tWish' )
    [ eparam, tstats, logL, optimoutput ] = matvtWishest(X,x0,varargin{:});
elseif strcmp( dist, 'itWish' )
    [ eparam, tstats, logL, optimoutput ] = matvitWishest(X,x0,varargin{:});
elseif strcmp( dist, 'F' )
    [ eparam, tstats, logL, optimoutput ] = matvFest(X,x0,varargin{:});   
elseif strcmp( dist, 'Riesz' )
    [ eparam, tstats, logL, optimoutput ] = matvRieszest(X,x0,varargin{:});
elseif strcmp( dist, 'iRiesz' )
    [ eparam, tstats, logL, optimoutput ] = matviRieszest(X,x0,varargin{:});
elseif strcmp( dist, 'tRiesz' )
    [ eparam, tstats, logL, optimoutput ] = matvtRieszest(X,x0,varargin{:});
elseif strcmp( dist, 'itRiesz' )
    [ eparam, tstats, logL, optimoutput ] = matvitRieszest(X,x0,varargin{:});
elseif strcmp( dist, 'FRiesz' )
    [ eparam, tstats, logL, optimoutput ] = matvFRieszest(X,x0,varargin{:});
end

end

