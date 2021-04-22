function [ e, h, tau, V_m ] = mf2_garch_core(parameters, y, m)


    mu = parameters(1);
    alpha = parameters(2);
    gamma = parameters(3);
    beta = parameters(4);
    
    lambda_0 = parameters(5);
    lambda_1 = parameters(6);
    lambda_2 = parameters(7);

    h = ones(size(y));
    
    tau = ones(size(y))*mean(y.^2);
    V = zeros(size(y));
    V_m = zeros(size(y));
    
    
    
    for t = 2:m
        
       if ((y(t-1) - mu) < 0) 
       h(t) = (1-alpha-gamma/2-beta) + (alpha + gamma).*(y(t-1) - mu).^2 ./ tau(t-1) + beta .* h(t-1);
       else
       h(t) = (1-alpha-gamma/2-beta) + alpha .* (y(t-1) - mu).^2 ./ tau(t-1) + beta .* h(t-1); 
       end
       
    end
    
    
    for t = (m+1):length(y)
        
       if ((y(t-1) - mu) < 0) 
       h(t) = (1-alpha-gamma/2-beta) + (alpha + gamma).*(y(t-1) - mu).^2 ./ tau(t-1) + beta .* h(t-1);
       else
       h(t) = (1-alpha-gamma/2-beta) + alpha .* (y(t-1) - mu).^2 ./ tau(t-1) + beta .* h(t-1); 
       end
       
       V(t) = (y(t) - mu).^2 ./ h(t);
       V_m(t) = sum(V(t-(m-1):t))./m;
       
       tau(t) = lambda_0 + lambda_1 * V_m(t-1) + lambda_2 * tau(t-1);
       
    end    
    
    e = (y-mu) ./ sqrt(h.*tau);
    
    sample = (2*252+1):length(y);
    h=h(sample);
    tau=tau(sample);
    e=e(sample);
    
    end

