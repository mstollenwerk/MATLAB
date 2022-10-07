function x0 = matvX0df(dist,p)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
x0_df_barlett = 2*p;
x0_df_barlettL = 2*p*ones(1,p);
x0_df_barlettU = 2*p*ones(1,p);
x0_df_chi = 5;

if strcmp( dist, 'Wish' )
    x0 = x0_df_barlett;
elseif strcmp( dist, 'iWish' )
    x0 = x0_df_barlett;
elseif strcmp( dist, 'tWish' )
    x0 = [ x0_df_barlett, x0_df_chi ];
elseif strcmp( dist, 'itWish' )
    x0 = [ x0_df_chi, x0_df_barlett ];
elseif strcmp( dist, 'F' )
    x0 = [ x0_df_barlett, x0_df_barlett ];
elseif strcmp( dist, 'Riesz' )
    x0 = x0_df_barlettL;
elseif strcmp( dist, 'iRiesz2' )
    x0 = x0_df_barlettU;
elseif strcmp( dist, 'tRiesz' )
    x0 = [ x0_df_barlettL, x0_df_chi ];
elseif strcmp( dist, 'itRiesz2' )
    x0 = [ x0_df_chi, x0_df_barlettU ];
elseif strcmp( dist, 'FRiesz' )
    x0 = [ x0_df_barlettL, x0_df_barlettU ];
elseif strcmp( dist, 'iFRiesz2' )
    x0 = [ x0_df_barlettL, x0_df_barlettU ];
end

end

