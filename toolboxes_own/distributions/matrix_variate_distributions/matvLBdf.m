function lb_df = matvLBdf(dist,k)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
lb_df_barlett = [(k-1)];
lb_df_barlettL = [(0:k-1)];
lb_df_barlettU = [flip(0:k-1)];
lb_df_chi = [0];

if strcmp( dist, 'Wish' )
    lb_df = lb_df_barlett;
%     Aeq = [ zeros(4,6) eye(4) ];
%     beq = zeros(4,1);       
elseif strcmp( dist, 'iWish' )
    lb_df = lb_df_barlett;
%     Aeq = [ zeros(4,5*p+1) eye(4) ];
%     beq = zeros(4,1);    
elseif strcmp( dist, 'tWish' )
    lb_df = [ lb_df_barlett, lb_df_chi ];
%     Aeq = [ zeros(4,6) eye(4) zeros(4,5) ;
%             zeros(4,11) eye(4)];
%     beq = zeros(8,1);    
elseif strcmp( dist, 'itWish' )
    lb_df = [ lb_df_chi, lb_df_barlett ];
%     Aeq = [ zeros(4,6) eye(4) zeros(4,5) ;
%             zeros(4,11) eye(4)];
%     beq = zeros(8,1);      
elseif strcmp( dist, 'F' )
    lb_df = [ lb_df_barlett, lb_df_barlett ]; 
%     Aeq = [ zeros(4,6) eye(4) zeros(4,5) ;
%             zeros(4,11) eye(4)];
%     beq = zeros(8,1);        
elseif strcmp( dist, 'Riesz' )
    lb_df = lb_df_barlettL;
elseif strcmp( dist, 'iRiesz2' )
    lb_df = lb_df_barlettU;
elseif strcmp( dist, 'tRiesz' )
    lb_df = [ lb_df_barlettL, lb_df_chi ];
elseif strcmp( dist, 'itRiesz2' )
    lb_df = [ lb_df_chi, lb_df_barlettU];
elseif strcmp( dist, 'FRiesz' )
    lb_df = [ lb_df_barlettL, lb_df_barlettU];
elseif strcmp( dist, 'iFRiesz2' )
    lb_df = [ lb_df_barlettL, lb_df_barlettU ];
end

end