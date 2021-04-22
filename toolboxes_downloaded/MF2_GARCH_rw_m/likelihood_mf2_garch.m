function [ll]=likelihood_mf2_garch(param, y, m)
    [e, h, tau, V_m ] = mf2_garch_core(param, y, m);

    lls = -0.5*(log(2*pi) + log(h.*tau) + e.^2);
    
    % Use these to comput the LL
    ll = -sum(lls);
    
end