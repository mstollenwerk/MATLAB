classdef LogLik_I % in LogLik.m
    methods (Static = true)
        
        
        %% multi GAMMA in logs
        function [mg] = lnmultigamma(k, a)
            % return multivariate gamma
            mg= log(pi^(k*(k-1)/4));
            
            for i=1:k
                mg= mg + gammaln(a + (1-i)/2);
            end
            
        end
        
        
        %% multi GAMMA
        function [mg] = multigamma(k, a)
            % return multivariate gamma
            mg= pi^(k*(k-1)/4);
            
            for i=1:k
                mg= mg*gamma(a + (1-i)/2);
            end
            
        end
        
        function [LLF,LLF_r2,LLF_RC] = LogLik_GAS_tF_HAR_CT(k,T,params,r2,RC);
            
            alpha = params(1);
            beta_vec  = params(2:4);
            nu_t  = params(end-2);
            nu_f1 = params(end-1);
            nu_f2 = params(end);
            
            constant = nu_f1/ (nu_f2-k-1);
            [V,~,Flag]  = Filter_New_I.GAS_HAR_CT( k, T, [alpha beta_vec' nu_t nu_f1 nu_f2], r2, RC );
            
            if Flag==1
                LLF = 1e12;
                LLF_RC = 1e12;
                LLF_r2 = 1e12;
            else
                % get loglike
                loglik_r2= zeros(1,T);
                loglik_RC= zeros(1,T);
                
                for j=1:T
                    j_invV       = inv(V(:,:,j));
                    j_logdetV    = log(det(V(:,:,j)));
                    j_RC         = RC(:,:,j);
                    j_logdetRC   = log(det(RC(:,:,j)));
                    
                    % Student T
                    loglik_r2(j) = - 1/2 * j_logdetV  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(j_invV*r2(:,:,j)));
                    % matrix F
                    loglik_RC(j) = 0.5*nu_f1 * log(det(constant*j_invV)) + 0.5*(nu_f1-k-1)*j_logdetRC - 0.5*(nu_f1+nu_f2)*log(det(eye(k)+ constant*j_invV*j_RC));
                end
                
                loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);
                loglik_RC = loglik_RC + LogLik_I.lnmultigamma(k, 0.5*(nu_f1+nu_f2)) - LogLik_I.lnmultigamma(k, 0.5*nu_f1) - LogLik_I.lnmultigamma(k, 0.5*nu_f2);
                
                LLF_r2 = sum(loglik_r2);
                LLF_RC = sum(loglik_RC);
                
                LLF= -sum(loglik_r2 + loglik_RC);
                
                if  isnan(LLF)
                    LLF=1e6;
                end
                
            end
            
        end
                        
        
        function [LLF,LLF_r2,LLF_RC,V] = LogLik_GAS_tF_CT(k,T,params,r2,RC);
            
            
            alpha = params(1);
            beta  = params(2);
            nu_t  = params(3);
            nu_f1 = params(4);
            nu_f2 = params(5);
            
            constant = nu_f1/ (nu_f2-k-1);
            
            [V,~,Flag]  = Filter_New_I.GAS_CT( k, T, [alpha beta nu_t nu_f1 nu_f2], r2, RC );
            
            if Flag==1
                LLF = 1e12;
                LLF_r2 = 1e12;
                LLF_RC = 1e12;
                
            else
                % get loglike
                loglik_r2= zeros(1,T);
                loglik_RC= zeros(1,T);
                
                for j=1:T
                    j_invV       = inv(V(:,:,j));
                    j_logdetV    = log(det(V(:,:,j)));
                    j_RC         = RC(:,:,j);
                    j_logdetRC   = log(det(RC(:,:,j)));
                    
                    %loglik_r2(j) = - 1/2 * j_logdetV  - 1/2 * trace( j_invV *r2(:,:,j) ); normal
                    
                    % Student T
                    loglik_r2(j) = - 1/2 * j_logdetV  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(j_invV*r2(:,:,j)));
                    % matrix F
                    loglik_RC(j) = 0.5*nu_f1 * log(det(constant*j_invV)) + 0.5*(nu_f1-k-1)*j_logdetRC - 0.5*(nu_f1+nu_f2)*log(det(eye(k)+ constant*j_invV*j_RC));
                end
                
                %loglik_r2= loglik_r2 - k/2 * log(2*pi); % normal
                
                loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);
                loglik_RC = loglik_RC + LogLik_I.lnmultigamma(k, 0.5*(nu_f1+nu_f2)) - LogLik_I.lnmultigamma(k, 0.5*nu_f1) - LogLik_I.lnmultigamma(k, 0.5*nu_f2);
                
                LLF_r2 = sum(loglik_r2);
                LLF_RC = sum(loglik_RC);
                
                LLF= -sum(loglik_r2 + loglik_RC);
                
                if  isnan(LLF)
                    LLF=1e6;
                end
                
            end
            
        end
        
        
        
        
        
        
    end
end