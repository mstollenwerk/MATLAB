function [ nLogL, logLcontr, weights, scores_weights ] = ...
    gas_weights_mix_2_likeRec( param, L1, L2 )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 05.04.2022

T = length(L1);
t_ahead = 220;
%% Parameters
%% Data Storage
weights = NaN(T+t_ahead,1);
scores_weights = NaN(T+t_ahead,1);
logLcontr = NaN(T,1);
%% Recursion
weights(1) = param(1);
for tt=1:T
    
    if tt > 1
        weights(tt) = (1-param(3))*param(2) ...
                       + param(3)*weights(tt-1) ...
                       + param(4)*scores_weights(tt-1);          
    end
    
    if weights(tt) > 2 || weights(tt) < -1
        nLogL = inf;
        logLcontr = NaN;
        weights = NaN;
        scores_weights = NaN;
        return
    end
    
    Like_tt = weights(tt)*L1(tt) + (1 - weights(tt))*L2(tt);
    scores_weights(tt) = (weights(tt)-weights(tt)^2)*(L1(tt) - L2(tt))/Like_tt;
    
    logLcontr(tt) = log(Like_tt);
    
end

%% Fcst
for tt=T+1:T+t_ahead
    
    if tt == T+1
        weights(tt) = (1-param(3))*param(2) ...
                       + param(3)*weights(tt-1) ...
                       + param(4)*scores_weights(tt-1);  
                   
    else
        
        weights(tt) = (1-param(3))*param(2) ...
                       + param(3)*weights(tt-1);        
                   
    end

end
%% Log-Likelihood
nLogL = -sum(logLcontr);
