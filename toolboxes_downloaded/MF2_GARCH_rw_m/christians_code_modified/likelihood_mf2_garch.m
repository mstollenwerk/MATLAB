function [ll]=likelihood_mf2_garch(param, rv, m)
    [esq, h, tau, V_m ] = mf2_garch_core(param, rv, m);
    
    try
        lls = -0.5*(log(2*pi) + log(h.*tau) + esq);
    catch
        ll = inf;
    end
    
    % Use these to comput the LL
    ll = -sum(lls);
    
end