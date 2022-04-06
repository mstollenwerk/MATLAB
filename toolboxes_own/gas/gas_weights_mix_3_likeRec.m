function [ nLogL, logLcontr, weights, scores_weights ] = ...
    gas_weights_mix_3_likeRec( param, logL1, logL2, logL3 )
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 05.04.2022

T = length(logL1);
t_ahead = 220;
%% Parameters
param_1 = param(1:4);
param_2 = param(5:8);
%% Data Storage
weights{1} = NaN(T+t_ahead,1);
weights{2} = NaN(T+t_ahead,1);
scores_weights{1} = NaN(T+t_ahead,1);
scores_weights{2} = NaN(T+t_ahead,1);
logLcontr = NaN(T,1);
%% Recursion
weights{1}(1) = param_1(1);
weights{2}(1) = param_2(1);
for tt=1:T
    
    if tt > 1
        weights{1}(tt) = (1-param_1(3))*param_1(2) ...
                       + param_1(3)*weights{1}(tt-1) ...
                       + param_1(4)*scores_weights{1}(tt-1);

        weights{2}(tt) = (1-param_2(3))*param_2(2) ...
                       + param_2(3)*weights{2}(tt-1) ...
                       + param_2(4)*scores_weights{2}(tt-1);               
    end
    
    if weights{1}(tt) > 2 || weights{1}(tt) < -1 || weights{2}(tt) > 2 || weights{2}(tt) < -1
        nLogL = inf;
        logLcontr = NaN;
        weights = NaN;
        scores_weights = NaN;
        return
    end
    
    Like_tt = weights{1}(tt)*logL1(tt) + weights{2}(tt)*logL2(tt) + (1 - weights{1}(tt) - weights{2}(tt))*logL3(tt);
    wertebereich_1_grenzen = [-1; 1; ...
                               logL1(tt)/(weights{1}(tt)*logL1(tt) + weights{2}(tt)*logL2(tt)); ...
                               -logL3(tt)/(weights{2}(tt)*logL2(tt) + (1-weights{1}(tt)-weights{2}(tt))*logL3(tt)); ...
                               (logL1(tt) - logL3(tt))/(weights{1}(tt)*logL1(tt) + (1-weights{1}(tt)-weights{2}(tt))*logL3(tt))];
    scores_weights{1}(tt) = (logL1(tt) - logL3(tt))/Like_tt/(max(wertebereich_1_grenzen)-min(wertebereich_1_grenzen)); %(weights{1}(tt)-weights{1}(tt)^2)*

    wertebereich_2_grenzen = [-1; 1; ...
                               logL2(tt)/(weights{2}(tt)*logL2(tt) + weights{1}(tt)*logL1(tt)); ...
                               -logL3(tt)/(weights{1}(tt)*logL1(tt) + (1-weights{1}(tt)-weights{2}(tt))*logL3(tt)); ...
                               (logL2(tt) - logL3(tt))/(weights{2}(tt)*logL2(tt) + (1-weights{1}(tt)-weights{2}(tt))*logL3(tt))];
    
    scores_weights{2}(tt) = (logL2(tt) - logL3(tt))/Like_tt/(max(wertebereich_2_grenzen)-min(wertebereich_2_grenzen)); %(weights{2}(tt)-weights{2}(tt)^2)*
    logLcontr(tt) = log(Like_tt);
    
end

%% Fcst
for tt=T+1:T+t_ahead
    
    if tt == T+1
        weights{1}(tt) = (1-param_1(3))*param_1(2) ...
                       + param_1(3)*weights{1}(tt-1) ...
                       + param_1(4)*scores_weights{1}(tt-1);

        weights{2}(tt) = (1-param_2(3))*param_2(2) ...
                       + param_2(3)*weights{2}(tt-1) ...
                       + param_2(4)*scores_weights{2}(tt-1);                
    else
        
        weights{1}(tt) = (1-param_1(3))*param_1(2) ...
                       + param_1(3)*weights{1}(tt-1);

        weights{2}(tt) = (1-param_2(3))*param_2(2) ...
                       + param_2(3)*weights{2}(tt-1);                     
    end

end
%% Log-Likelihood
nLogL = -sum(logLcontr);
