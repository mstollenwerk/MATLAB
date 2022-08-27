classdef LogLik_MV % in LogLik.m
    methods (Static = true)
        
        
        %% multi GAMMA in logs
        function [mg] = lnmultigamma(k, a)
            % retfurn multivariate gamma
            mg= (k*(k-1)/4)*log(pi);
            
            for i=1:k
                mg= mg + gammaln(a + (1-i)/2);
            end
            
        end
        
        function [mg] = lnmultigamma_nu_vec(k, nu_vec)
            % return multivariate gamma
            mg= (k*(k-1)/4)*log(pi);
            
            for i=1:k
                mg= mg + gammaln(nu_vec(i) + (1-i)/2);
            end
            
        end
        
        function [mg] = lnmultigamma_U_nu_vec(k, nu_vec)
            % return upper multivariate gamma
            mg= (k*(k-1)/4)*log(pi);
            
            for i=1:k
                mg= mg + gammaln(nu_vec(i) + (i-k)/2);
            end
            
        end
        
        
        function [ln_PWD_V] = ln_PWD(V, nu_vec)
            % the lower power weighted determinant
            
            L = chol(V)';
            ln_PWD_V = 2 * nu_vec' * log(diag(L));
            
        end
        
        
        function [ln_PWD_V] = ln_PWD_U(V, nu_vec)
            % the UPPER power weighted determinant
            
            ln_PWD_V = NaN;
            k = size(V,1);
            [chol_V_inv,ind_pd] = chol(V\eye(k));
            if ind_pd==0
                U = chol_V_inv\eye(k);
                ln_PWD_V = 2 * nu_vec' * log(diag(U));
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
        
        function [LLF,loglik_RC] = LogLik_HAR_IW(k,T,params,RC);
            
            %omega = params(1);
            beta_1 = params(end-3);
            beta_2  = params(end-2);
            beta_3  = params(end-1);
            nu = params(end);
            
            [V,Flag]  = Filter_New.IW_HAR(k, T, [beta_1 beta_2 beta_3],RC );
            
            loglik_RC= zeros(1,T);
            constant = (nu-k-1);
            
            if Flag ==1
                loglik_RC = [];
                LLF = 1e12;
                
            else
                
                for j =1 :T
                    j_RC        = RC(:,:,j);
                    j_detRC     = det(j_RC);
                    j_inv_RC    = inv(j_RC);
                    V_j         = V(:,:,j);
                    j_detV      = det(V_j);
                    
                    loglik_RC(j) = (nu/2) * log(j_detV) -0.5*(nu+k+1)*log(j_detRC)-constant/2 * trace(j_inv_RC*V_j);
                    
                end
                
                loglik_RC = loglik_RC + 0.5*nu*k*log(constant/2) - LogLik.lnmultigamma(k, 0.5*nu);
                LLF = -sum(loglik_RC);
                
                if isnan(LLF)
                    LLF = 1e12;
                    
                end
            end
        end
        
        function [LLF,LLF_r2,LLF_RC] = LogLik_GAS_HAR_CT_tF(k,T,params,r2,RC);
            
            alpha = params(1);
            beta_vec  = params(2:4);
            nu_t  = params(end-2);
            nu_f1 = params(end-1);
            nu_f2 = params(end);
            
            constant = nu_f1/ (nu_f2-k-1);
            [V,~,Flag]  = Filter_New.GAS_HAR_CT( k, T, [alpha beta_vec' nu_t nu_f1 nu_f2], r2, RC );
            
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
                loglik_RC = loglik_RC + LogLik.lnmultigamma(k, 0.5*(nu_f1+nu_f2)) - LogLik.lnmultigamma(k, 0.5*nu_f1) - LogLik.lnmultigamma(k, 0.5*nu_f2);
                
                LLF_r2 = sum(loglik_r2);
                LLF_RC = sum(loglik_RC);
                
                LLF= -sum(loglik_r2 + loglik_RC);
                
                if  isnan(LLF)
                    LLF=1e6;
                end
                
            end
            
        end
        
        function [LLF,loglik_RC,V] = LogLik_HEAVY_matrix_F(k,T,T_w,params,RC);
            
            alpha = params(1);
            beta  = params(2);
            nu_f1 = params(end-1);
            nu_f2 = params(end);
            
            constant = nu_f1/ (nu_f2-k-1);
            [V,Flag]  = Filter_New.CAW(k,1,T,[alpha beta]', RC );
            
            if Flag==1
                LLF = 1e12;
            else
                % get loglike
                loglik_RC= zeros(1,T);
                
                for j=1:T
                    j_invV       = inv(V(:,:,j));
                    j_RC         = RC(:,:,j);
                    j_logdetRC   = log(det(RC(:,:,j)));
                    
                    % matrix F
                    loglik_RC(j) = 0.5*nu_f1 * log(det(constant*j_invV)) + 0.5*(nu_f1-k-1)*j_logdetRC - 0.5*(nu_f1+nu_f2)*log(det(eye(k)+ constant*j_invV*j_RC));
                end
                
                loglik_RC = loglik_RC + LogLik.lnmultigamma(k, 0.5*(nu_f1+nu_f2)) - LogLik.lnmultigamma(k, 0.5*nu_f1) - LogLik.lnmultigamma(k, 0.5*nu_f2);
                
                LLF = -sum(loglik_RC(T_w+1:end));
                
                if  isnan(LLF)
                    LLF=1e6;
                end
                
            end
            
        end
        
        
        function [LLF,loglik_RC,H_TOT_mat,H_ST_mat] = LogLik_HEAVY_Matrix_F_LT_spec(k,T,T_w,params,RC_mat);
            
            
            para_HEAVY_ST = params(1:2);
            nu_f1 = params(end-1);
            nu_f2 = params(end);
            
            constant = nu_f1/ (nu_f2-k-1);
            I_k = eye(k);
            
            H_bar = mean(RC_mat,3);
            chol_H_LT_bar = chol(H_bar)';
            inv_chol_H_LT_bar = chol_H_LT_bar\I_k;
            
            
            H_ST_mat = zeros(k,k,T);
            H_TOT_mat = zeros(k,k,T);
            H_ST_mat(:,:,1) = eye(k);
            Omega_ST = (1-sum(para_HEAVY_ST))*I_k;
            loglik_RC= zeros(T,1);
            Flag = 0;
            
            for j =1 :T
                
                RC_j = RC_mat(:,:,j);
                
                % compute H_tot
                H_TOT_mat(:,:,j) = chol_H_LT_bar  * H_ST_mat(:,:,j) * chol_H_LT_bar';
                
                % compute V_t
                
                [~,ind_pd_H_ST] = chol(H_ST_mat(:,:,j));
                if ind_pd_H_ST>0
                    Flag = 1;
                    loglik_RC = NaN;
                    break
                end
                
                if j<T
                    
                    % update H_ST
                    H_ST_mat(:,:,j+1) = Omega_ST + para_HEAVY_ST(1)*(inv_chol_H_LT_bar*RC_j)*inv_chol_H_LT_bar' + para_HEAVY_ST(2) * H_ST_mat(:,:,j);
                    
                    
                    j_invV = H_TOT_mat(:,:,j)\I_k;
                    j_logdetRC = log(det(RC_j));
                    
                    loglik_RC(j,1) = 0.5*nu_f1 * log(det(constant*j_invV)) + 0.5*(nu_f1-k-1)*j_logdetRC - 0.5*(nu_f1+nu_f2)*log(det(eye(k)+ constant*j_invV*RC_j));
                end
                
            end
            
            if Flag==1
                LLF = 1e12;
            else
                
                loglik_RC = loglik_RC + LogLik_MV.lnmultigamma(k, 0.5*(nu_f1+nu_f2)) - LogLik_MV.lnmultigamma(k, 0.5*nu_f1) - LogLik_MV.lnmultigamma(k, 0.5*nu_f2);
                LLF = -sum(loglik_RC(T_w+1:end));
                
                if isnan(LLF)
                    LLF = 1e12;
                end
            end
            
            
        end
        
        
        function [LLF,loglik_RC,H_ST_mat,Z_t_mat,H_TOT_mat] = LogLik_MF2_HEAVY_matrix_F(k,T,T_w,params,RC);
            
            alpha_ST = params(1);
            beta_ST  = params(2);
            
            nu_f1 = params(3);
            nu_f2 = params(4);
            
            alpha_LT = params(end-1);
            beta_LT  = params(end);
            
            I_k = eye(k);
            constant = nu_f1/ (nu_f2-k-1);
            
            H_LT_bar = mean(RC,3);
            
            Omega_LT = H_LT_bar * (1-alpha_LT - beta_LT);
            chol_H_LT_bar = chol(H_LT_bar)';
            inv_chol_H_LT_bar = chol_H_LT_bar\I_k;
            
            chol_Z_new = chol_H_LT_bar;
            inv_chol_Z_new = inv_chol_H_LT_bar;
            
            Z_t_mat = zeros(k,k,T);
            Z_t_mat(:,:,1) = H_LT_bar;
            V_t_mat = zeros(k,k,T);
            
            H_ST_mat = zeros(k,k,T);
            H_TOT_mat = zeros(k,k,T);
            H_ST_mat(:,:,1) = eye(k);
            Omega_ST = (1-alpha_ST-beta_ST)*I_k;
            loglik_RC= zeros(T,1);
            Flag = 0;
            
            for j =1 :T
                
                j_RC         = RC(:,:,j);
                
                % compute H_tot
                H_TOT_mat_j= chol_Z_new * H_ST_mat(:,:,j) * chol_Z_new';
                
                if rcond(H_TOT_mat_j)<1e-12
                    Flag = 1;
                    loglik_RC = NaN;
                    break
                end
                
                H_TOT_mat(:,:,j) = H_TOT_mat_j;
                
                % compute V_t
                
                [chol_H_ST_mat_j,ind_pd_H_ST] = chol(H_ST_mat(:,:,j));
                if ind_pd_H_ST>0
                    Flag = 1;
                    loglik_RC = NaN;
                    break
                end
                
                chol_H_ST_j = chol_H_ST_mat_j';
                inv_chol_H_ST_j = chol_H_ST_j\eye(k);
                V_t_mat(:,:,j) = inv_chol_H_ST_j * j_RC* inv_chol_H_ST_j';
                
                if j<T
                    
                    % update H_ST
                    H_ST_mat(:,:,j+1) = Omega_ST +  alpha_ST*(inv_chol_Z_new*j_RC)*inv_chol_Z_new' +  beta_ST * H_ST_mat(:,:,j);
                    
                    % update Z_t
                    if j>=T_w
                        %Z_t_mat(:,:,j+1) = Omega_LT + A_mat_LT * mean(V_t_mat(:,:,j-T_w+1:j),3)*A_mat_LT' + para_HEAVY_LT(3)*Z_t_mat(:,:,j);
                        Z_t_mat(:,:,j+1) = Omega_LT + alpha_LT*mean(V_t_mat(:,:,j-T_w+1:j),3) + beta_LT*Z_t_mat(:,:,j);
                        Z_t_mat_new = Z_t_mat(:,:,j+1);
                    else
                        Z_t_mat(:,:,j+1) = H_LT_bar;
                        Z_t_mat_new = H_LT_bar;
                    end
                    
                    
                    [chol_Z_new_j,ind_pd_Z_j] = chol(Z_t_mat_new);
                    if ind_pd_Z_j>0
                        Flag = 1;
                        loglik_RC = NaN;
                        break
                    end
                    
                    
                    chol_Z_new = chol_Z_new_j';
                    inv_chol_Z_new = chol_Z_new\eye(k);
                    
                    
                    inv_Sigma_t = H_TOT_mat(:,:,j)\eye(k);
                    j_logdetRC  = log(det(j_RC));
                    
                    % matrix F
                    loglik_RC(j) = 0.5*nu_f1 * log(det(constant*inv_Sigma_t)) + 0.5*(nu_f1-k-1)*j_logdetRC - 0.5*(nu_f1+nu_f2)*log(det(eye(k)+ constant*inv_Sigma_t*j_RC));
                end
                
            end
            
            loglik_RC = loglik_RC + LogLik_MV.lnmultigamma(k, 0.5*(nu_f1+nu_f2)) - LogLik_MV.lnmultigamma(k, 0.5*nu_f1) - LogLik_MV.lnmultigamma(k, 0.5*nu_f2);
            
            if Flag==1
                LLF = 1e12;
            else
                
                LLF = -sum(loglik_RC(T_w+1:end));
                
                if isnan(LLF)
                    LLF = 1e12;
                end
            end
            
        end
        
        
        
        
        function [LLF,V] = LogLik_CAW(k,p,T, params,RC)
            % for now: diagonal CAW with targeting
            % p = order of AR and MA, usually 1
            
            nu = params(end);
            %[V,Flag] = Filter_New.CAW_diag(k,p,T, params,RC);
            [V,Flag] = Filter_New.CAW(k,p,T, params,RC);
            
            if Flag==1
                LLF = 1e6;
                
            else
                
                % get loglike
                loglik_RC= zeros(1,T);
                
                for j=1:T
                    j_RC         = RC(:,:,j);
                    j_invV       = inv(V(:,:,j));
                    j_logdetV    = log(det(V(:,:,j)));
                    j_logdetRC   = log(det(RC(:,:,j)));
                    
                    loglik_RC(j) = - nu/2 * j_logdetV   +  (nu-k-1)/2 * j_logdetRC  - 0.5 * nu * trace(j_invV * j_RC);
                end
                
                loglik_RC= loglik_RC -(k*nu)/2 * log(2)  + (k*nu)/2*log(nu) -  LogLik.lnmultigamma(k, nu/2 );
                
                LLF= -sum(loglik_RC);
                
                if  isnan(LLF) || isreal(LLF)==0
                    LLF=1e6;
                end
            end
        end
        
        
        %% log lik Wishart
        function [LLF,LLF_r2,LLF_RC,V] = LogLik_GAS_NW_CT(k, T, params, r2, RC)
            % PURPOSE:
            %     Gaussian loglikelihood for the Time-varying CAPM model: _TVCAPM
            
            % Gamma parameters
            
            %alpha= exp(params(1));
            %beta = exp(params(2)) / (1+exp(params(2))) - 0.0001;
            %nu   = exp(params(3)) + k + 0.01 - 1;
            alpha= params(1);
            beta = params(2);
            nu   = params(3);
            
            
            [ V, ~, ~,Flag]  = Filter.GAS_CT( k, T, [alpha beta nu], r2, RC );
            
            if Flag==1
                LLF = 1e6;
                
            else
                
                % get loglike
                loglik_r2= zeros(1,T);
                loglik_RC= zeros(1,T);
                %data_mat = zeros(T,6);
                
                for j=1:T
                    j_RC         = RC(:,:,j);
                    j_invV       = inv(V(:,:,j));
                    j_logdetV    = log(det(V(:,:,j)));
                    j_logdetRC   = log(det(RC(:,:,j)));
                    
                    %data_mat(j,:) = [rcond(j_RC) rcond(r2(:,:,j)) j_logdetRC  j_logdetV trace(j_invV * j_RC)  trace(j_invV * r2(:,:,j))];
                    
                    loglik_r2(j) = - 1/2 * j_logdetV  - 1/2 * trace( j_invV * r2(:,:,j) );
                    loglik_RC(j) = - nu/2 * j_logdetV   +  (nu-k-1)/2 * j_logdetRC  - 0.5 * nu * trace( j_invV * j_RC);
                end
                
                loglik_r2= loglik_r2 - k/2 * log(2*pi);
                loglik_RC= loglik_RC -(k*nu)/2 * log(2)  + (k*nu)/2*log(nu) -  LogLik.lnmultigamma(k, nu/2 );
                
                LLF= -sum(loglik_r2 + loglik_RC);
                
                LLF_r2 = sum(loglik_r2);
                LLF_RC = sum(loglik_RC);
                
                
                if  isnan(LLF) || isreal(LLF)==0
                    LLF=1e6;
                end
            end
        end
        
        
        
        
        function [LLF,LLF_r2,LLF_RC,V] = LogLik_GAS_CT_tF(k,T,params,return_mat,RC);
            
            
            alpha = params(1);
            beta  = params(2);
            nu_t  = params(3);
            nu_f1 = params(4);
            nu_f2 = params(5);
            
            constant = nu_f1/ (nu_f2-k-1);
            
            [V,~,Flag]  = Filter_New.GAS_tF_CT( k, T, [alpha beta nu_t nu_f1 nu_f2], return_mat, RC );
            
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
                    r2_j         = return_mat(j,:)'*return_mat(j,:);
                    
                    %loglik_r2(j) = - 1/2 * j_logdetV  - 1/2 * trace( j_invV *r2(:,:,j) ); normal
                    
                    % Student T
                    loglik_r2(j) = - 1/2 * j_logdetV  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(j_invV*r2_j));
                    % matrix F
                    loglik_RC(j) = 0.5*nu_f1 * log(det(constant*j_invV)) + 0.5*(nu_f1-k-1)*j_logdetRC - 0.5*(nu_f1+nu_f2)*log(det(eye(k)+ constant*j_invV*j_RC));
                end
                
                %loglik_r2= loglik_r2 - k/2 * log(2*pi); % normal
                
                loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);
                loglik_RC = loglik_RC + LogLik.lnmultigamma(k, 0.5*(nu_f1+nu_f2)) - LogLik.lnmultigamma(k, 0.5*nu_f1) - LogLik.lnmultigamma(k, 0.5*nu_f2);
                
                LLF_r2 = sum(loglik_r2);
                LLF_RC = sum(loglik_RC);
                
                LLF= -sum(loglik_r2 + loglik_RC);
                
                if  isnan(LLF)
                    LLF=1e6;
                end
                
            end
            
        end
        
        
        
        function [LLF,V] = LogLik_GAS_CT_t(k,T,params,return_mat);
            
            
            alpha = params(1);
            beta  = params(2);
            nu_t  = params(3);
            
            
            [V,~,Flag]  = Filter_New.GAS_t_CT( k, T, [alpha beta nu_t], return_mat);
            
            if Flag==1
                LLF = 1e12;
                
            else
                % get loglike
                loglik_r2= zeros(1,T);
                
                
                for j=1:T
                    j_invV       = inv(V(:,:,j));
                    j_logdetV    = log(det(V(:,:,j)));
                    r2_j         = return_mat(j,:)'*return_mat(j,:);
                    
                    % Student T
                    loglik_r2(j) = - 1/2 * j_logdetV  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(j_invV*r2_j));
                end
                
                %loglik_r2= loglik_r2 - k/2 * log(2*pi); % normal
                
                loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);
                
                LLF= -sum(loglik_r2);
                
                if  isnan(LLF)
                    LLF=1e6;
                end
                
            end
            
        end
        
        
        function [LLF,V,nu_vec,nu_vec_score] = LogLik_GAS_CT_t_tv_nu(k,T,params,return_mat);
            
            
            %alpha = params(1);
            %beta  = params(2);
            %nu_bar  = params(3);
            
            
            [V,~,Flag,nu_vec,nu_vec_score]  = Filter_New.GAS_t_CT_tv_nu( k, T, params, return_mat);
            
            if Flag==1
                LLF = 1e12;
                
            else
                % get loglike
                loglik_r2= zeros(T,1);
                
                
                for j=1:T
                    nu_t         = nu_vec(j,1);
                    j_invV       = inv(V(:,:,j));
                    j_logdetV    = log(det(V(:,:,j)));
                    r2_j         = return_mat(j,:)'*return_mat(j,:);
                    
                    % Student T
                    loglik_r2(j,1) = - 1/2 * j_logdetV  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(j_invV*r2_j));
                end
                
                %loglik_r2= loglik_r2 - k/2 * log(2*pi); % normal
                
                loglik_r2 = loglik_r2  + gammaln(0.5*(nu_vec + k)) - gammaln(nu_vec./2) - 0.5*k*log((nu_vec-2).*pi);
                
                LLF= -sum(loglik_r2);
                
                if  isnan(LLF)
                    LLF=1e6;
                end
                
            end
            
        end
        
        
        
        function [LLF,loglik_RC] = LogLik_Matrix_F(k,T,params,RC);
            
            nu1 = params(1);
            nu2 = params(2);
            
            %vechDelta = params(3:end);
            %Delta = Admin.Vech2Sym(k,vechDelta);
            Delta = mean(RC,3);
            
            [test1,test2] = chol(Delta);
            if test2>0
                LLF = 1e12;
            else
                
                inv_Delta = inv(Delta);
                loglik_RC= zeros(T,1);
                constant = nu1/ (nu2-k-1);
                for j =1 :T
                    j_RC = RC(:,:,j);
                    j_detRC = det(j_RC);
                    loglik_RC(j) = 0.5*(nu1-k-1)*log(j_detRC) - 0.5*(nu1+nu2)*log(det(eye(k)+ constant*inv_Delta*j_RC));
                    %loglik_RC(j) = 0.5*(nu1-k-1)*log(j_detRC) - 0.5*(nu1+nu2)*log(det(eye(k)+ inv_Delta*j_RC));
                    
                end
                
                loglik_RC = loglik_RC + LogLik.lnmultigamma(k, 0.5*(nu1+nu2)) - LogLik.lnmultigamma(k, 0.5*nu1) - LogLik.lnmultigamma(k, 0.5*nu2);
                loglik_RC = loglik_RC + 0.5*nu1 * log(det(constant*inv_Delta));
                
                LLF = -sum(loglik_RC);
                
                if isnan(LLF)
                    LLF = 1e12;
                    
                end
            end
        end
        
        
        function [LLF,loglik_RC] = LogLik_Wishart(k,T,params,RC);
            
            nu = params(1);
            
            Delta = mean(RC,3);
            log_det_Delta = log(det(Delta));
            inv_Delta = inv(Delta);
            
            loglik_RC= zeros(T,1);
            for j =1 :T
                j_RC         = RC(:,:,j);
                j_logdetRC   = log(det(RC(:,:,j)));
                loglik_RC(j) = - nu/2 * log_det_Delta   +  (nu-k-1)/2 * j_logdetRC  - 0.5 * nu * trace(inv_Delta * j_RC);
            end
            
            loglik_RC = loglik_RC -(k*nu)/2 * log(2)  + (k*nu)/2*log(nu) -  LogLik.lnmultigamma(k, nu/2 );
            LLF       = -sum(loglik_RC);
            
            if isnan(LLF)
                LLF = 1e12;
            end
            
        end
        
        
        function [LLF,loglik_RC] = LogLik_Riesz_I(k,T,params,RC);
            
            k_star = k*(k+1)/2;
            nu_vec = params(k_star+1:end);
            L_tilde = [params(1) 0; params(2) params(3)];
            
            L = L_tilde * diag(1./sqrt(nu_vec));
            Sigma = L * L';
            
            [chol_Sigma,ind_pd] = chol(Sigma);
            
            if ind_pd>0
                LLF = NaN;
            else
                
                ln_PWD_Sigma = LogLik.ln_PWD(Sigma,0.5*nu_vec);
                inv_Sigma = inv(Sigma);
                
                loglik_RC= zeros(T,1);
                for j =1 :T
                    j_RC       = RC(:,:,j);
                    ln_PWD_j_RC  = LogLik.ln_PWD(j_RC,0.5*(nu_vec-k-1));
                    
                    loglik_RC(j) = -ln_PWD_Sigma  + ln_PWD_j_RC - 0.5 * trace(inv_Sigma * j_RC);
                end
                
                loglik_RC = loglik_RC - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
                LLF       = -sum(loglik_RC);
                
            end
            if isnan(LLF)
                LLF = 1e12;
            end
            
        end
        
        function [LLF,loglik_RC] = LogLik_Riesz_II(k,T,params,RC);
            
            k_star = k*(k+1)/2;
            nu_vec = params(k_star+1:end);
            U_tilde = [params(1) params(2); 0 params(3)];
            
            U = U_tilde * diag(1./sqrt(nu_vec));
            Sigma = U * U';
            
            [chol_Sigma,ind_pd] = chol(Sigma);
            
            if ind_pd>0
                LLF = NaN;
            else
                
                ln_PWD_Sigma_U = LogLik.ln_PWD_U(Sigma,0.5*nu_vec);
                
                if isnan(ln_PWD_Sigma_U)
                    LLF = NaN;
                    
                else
                    inv_Sigma = inv(Sigma);
                    
                    loglik_RC= zeros(T,1);
                    for j =1 :T
                        j_RC       = RC(:,:,j);
                        ln_PWD_j_RC_U  = LogLik.ln_PWD_U(j_RC,0.5*(nu_vec-k-1));
                        
                        loglik_RC(j) = -ln_PWD_Sigma_U  + ln_PWD_j_RC_U - 0.5 * trace(inv_Sigma * j_RC);
                    end
                    
                    loglik_RC = loglik_RC - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_U_nu_vec(k, nu_vec/2 );
                    LLF       = -sum(loglik_RC);
                end
            end
            if isnan(LLF)
                LLF = 1e12;
            end
            
        end
        
        
        function [LLF,loglik_RC] = LogLik_Riesz_II_CT(k,T,params,RC);
            
            % CT_ means with covariance targeting
            
            nu_vec = params;
            Sigma_0= mean(RC,3);
            
            U_tilde = chol(Sigma_0\eye(k))\eye(k);
            U = U_tilde * diag(1./sqrt(nu_vec));
            Sigma = U * U';
            
            [chol_Sigma,ind_pd] = chol(Sigma);
            
            if ind_pd>0
                LLF = NaN;
            else
                
                ln_PWD_Sigma_U = LogLik.ln_PWD_U(Sigma,0.5*nu_vec);
                
                if isnan(ln_PWD_Sigma_U)
                    LLF = NaN;
                    
                else
                    inv_Sigma = inv(Sigma);
                    
                    loglik_RC= zeros(T,1);
                    for j =1 :T
                        j_RC       = RC(:,:,j);
                        ln_PWD_j_RC_U  = LogLik.ln_PWD_U(j_RC,0.5*(nu_vec-k-1));
                        
                        loglik_RC(j) = -ln_PWD_Sigma_U  + ln_PWD_j_RC_U - 0.5 * trace(inv_Sigma * j_RC);
                    end
                    
                    loglik_RC = loglik_RC - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_U_nu_vec(k, nu_vec/2 );
                    LLF       = -sum(loglik_RC);
                end
            end
            if isnan(LLF)
                LLF = 1e12;
            end
            
        end
        
        
        
        function [LLF,loglik_RC] = LogLik_Riesz_I_CT(k,T,params,RC);
            
            % CT_ means with covariance targeting
            
            nu_vec = params;
            Sigma_0= mean(RC,3);
            [chol_Sigma_0,ind_pd] = chol(Sigma_0);
            
            L_tilde = chol_Sigma_0';
            L = L_tilde * diag(1./sqrt(nu_vec));
            Sigma = L * L';
            
            [chol_Sigma,ind_pd] = chol(Sigma);
            
            if ind_pd>0
                LLF = NaN;
            else
                
                ln_PWD_Sigma = LogLik.ln_PWD(Sigma,0.5*nu_vec);
                inv_Sigma = inv(Sigma);
                
                loglik_RC= zeros(T,1);
                for j =1 :T
                    j_RC       = RC(:,:,j);
                    ln_PWD_j_RC  = LogLik.ln_PWD(j_RC,0.5*(nu_vec-k-1));
                    
                    loglik_RC(j) = -ln_PWD_Sigma  + ln_PWD_j_RC - 0.5 * trace(inv_Sigma * j_RC);
                end
                
                loglik_RC = loglik_RC - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
                LLF       = -sum(loglik_RC);
                
            end
            if isnan(LLF)
                LLF = 1e12;
            end
            
        end
        
        
        function [LLF,loglik_RC] = LogLik_InvWishart(k,T,params,RC);
            
            %omega = params(1);
            %nr_l = k*(k+1)/2;
            %lower_vech_Sigma = params(1:nr_l);
            nu = params(end);
            Sigma = mean(RC,3);
            
            loglik_RC= zeros(1,T);
            constant = (nu-k-1);
            j_detSigma      = det(Sigma);
            
            for j =1 :T
                j_RC        = RC(:,:,j);
                j_detRC     = det(j_RC);
                j_inv_RC    = inv(j_RC);
                
                loglik_RC(j) = (nu/2) * log(j_detSigma) -0.5*(nu+k+1)*log(j_detRC)-constant/2 * trace(j_inv_RC*Sigma);
                
            end
            
            loglik_RC = loglik_RC + 0.5*nu*k*log(constant/2) - LogLik.lnmultigamma(k, 0.5*nu);
            LLF = -sum(loglik_RC);
            
            if isnan(LLF)
                LLF = 1e12;
                
            end
        end
        
        
        function [LLF,loglik_RC] = LogLik_Inv_Riesz_I(k,T,params,RC);
            
            %nr_l = k*(k+1)/2;
            %lower_vech_Sigma = params(1:nr_l);
            nu_vec = params(4:end);
            
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            L_00_trans = [params(1) 0; params(2) params(3)];
            Sigma_0 = L_00_trans * L_00_trans';
            inv_Sigma_0 = inv(Sigma_0);
            
            
            [L_0_trans,ind_pd] = chol(inv_Sigma_0);
            
            if ind_pd>0
                LLF = 1e14;
            else
                
                L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                inv_Sigma = L_trans' * L_trans;
                [~,ind_chol] = chol(inv_Sigma);
                
                loglik_RC= zeros(T,1);
                
                if ind_chol>0
                    LLF = 1e14;
                else
                    
                    Sigma = inv(inv_Sigma);
                    ln_PWD_Sigma = LogLik.ln_PWD(inv_Sigma,-0.5*nu_vec);
                    
                    for j =1 :T
                        j_RC       = RC(:,:,j);
                        j_inv_RC   = inv(j_RC);
                        
                        ln_PWD_j_RC  = LogLik.ln_PWD(j_inv_RC,0.5*(nu_vec+k+1));
                        loglik_RC(j) = ln_PWD_Sigma  + ln_PWD_j_RC - 0.5 * trace(Sigma * j_inv_RC);
                    end
                    
                    loglik_RC = loglik_RC - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
                    LLF       = -sum(loglik_RC);
                    
                    if isnan(LLF)
                        LLF = 1e12;
                    end
                    
                end
            end
        end
        
        function [LLF,loglik_RC] = LogLik_Inv_Riesz_I_CT(k,T,params,RC);
            
            % CT_ means with covariance targeting
            
            nu_vec = params;
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            Sigma_0 = mean(RC,3);
            inv_Sigma_0 = inv(Sigma_0);
            
            
            [L_0_trans,ind_pd] = chol(inv_Sigma_0);
            
            if ind_pd>0
                LLF = 1e14;
            else
                
                L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                inv_Sigma = L_trans' * L_trans;
                Sigma = inv(inv_Sigma);
                
                loglik_RC= zeros(T,1);
                
                [~,ind_chol] = chol(Sigma);
                if ind_chol>0
                    LLF = 1e14;
                else
                    
                    Sigma_inv = inv(Sigma);
                    ln_PWD_Sigma = LogLik.ln_PWD(Sigma_inv,-0.5*nu_vec);
                    
                    
                    for j =1 :T
                        j_RC       = RC(:,:,j);
                        j_inv_RC   = inv(j_RC);
                        %chol_j_RC_0  = chol(j_RC_0)';
                        %j_RC_01      = chol_j_RC_0 * diag(sqrt(nu_vec));
                        %j_RC         = j_RC_01 *j_RC_01';
                        ln_PWD_j_RC  = LogLik.ln_PWD(j_inv_RC,0.5*(nu_vec+k+1));
                        
                        loglik_RC(j) = ln_PWD_Sigma  + ln_PWD_j_RC - 0.5 * trace(Sigma * j_inv_RC);
                    end
                    
                    loglik_RC = loglik_RC - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
                    LLF       = -sum(loglik_RC);
                    
                    if isnan(LLF)
                        LLF = 1e12;
                    end
                    
                end
            end
        end
        
        
        function [LLF,loglik_RC] = LogLik_TRiesz_I_CT(k,T,params,y_mat);
            
            % CT_ means with covariance targeting
            
            nu_vec = params;
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            Sigma_0 = cov(y_mat);
            inv_Sigma_0 = inv(Sigma_0);
            
            
            [L_0_trans,ind_pd] = chol(inv_Sigma_0);
            
            if ind_pd>0
                LLF = 1e14;
            else
                
                L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                inv_Sigma = L_trans' * L_trans;
                Sigma = inv(inv_Sigma);
                
                loglik_RC= zeros(T,1);
                
                [~,ind_chol] = chol(Sigma);
                if ind_chol>0
                    LLF = 1e14;
                else
                    
                    %Sigma_inv = inv(Sigma);
                    %ln_PWD_Sigma = LogLik.ln_PWD(Sigma,0.5*nu_vec);
                    ln_PWD_Sigma = LogLik.ln_PWD(inv_Sigma,0.5*nu_vec);
                    
                    for j =1 :T
                        j_yy       = y_mat(j,:)'*y_mat(j,:);
                        j_yy_Sigma = j_yy + Sigma;
                        
                        %ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma,-0.5*(nu_vec+1));
                        %loglik_RC(j) = ln_PWD_Sigma  + ln_PWD_j_yy_Sigma;
                        
                        j_yy_Sigma_inv = j_yy_Sigma\eye(k);
                        
                        %ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma,-0.5*(nu_vec+1));
                        ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma_inv,0.5*(nu_vec+1));
                        
                        loglik_RC(j) = -ln_PWD_Sigma  + ln_PWD_j_yy_Sigma;
                        
                        
                    end
                    
                    loglik_RC = loglik_RC - (k/2)*log(pi) + LogLik.lnmultigamma_nu_vec(k, (nu_vec+1)/2 )   -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
                    LLF       = -sum(loglik_RC);
                    
                    if isnan(LLF)
                        LLF = 1e12;
                    end
                    
                end
            end
        end
        
        
        function [LLF,loglik_RC] = LogLik_TRiesz_I(k,T,params,y_mat);
            
            
            
            nu_vec = params(4:end);
            
            %mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            L_00_trans = [params(1) 0; params(2) params(3)];
            Sigma_0 = L_00_trans * L_00_trans';
            inv_Sigma_0 = inv(Sigma_0);
            
            %Sigma_0 = cov(y_mat);
            %inv_Sigma_0 = inv(Sigma_0);
            
            
            [L_0_trans,ind_pd] = chol(inv_Sigma_0);
            
            if ind_pd>0
                LLF = 1e14;
            else
                
                %L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                %inv_Sigma = L_trans' * L_trans;
                
                %Sigma = inv(inv_Sigma);
                
                Sigma = Sigma_0;
                inv_Sigma = inv_Sigma_0;
                
                loglik_RC= zeros(T,1);
                
                [~,ind_chol] = chol(Sigma);
                if ind_chol>0
                    LLF = 1e14;
                else
                    
                    %Sigma_inv = inv(Sigma);
                    %ln_PWD_Sigma = LogLik.ln_PWD(Sigma,0.5*nu_vec);
                    
                    ln_PWD_Sigma = LogLik.ln_PWD(inv_Sigma,0.5*nu_vec);
                    
                    for j =1 :T
                        j_yy       = y_mat(j,:)'*y_mat(j,:);
                        j_yy_Sigma = j_yy + Sigma;
                        
                        j_yy_Sigma_inv = j_yy_Sigma\eye(k);
                        
                        %ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma,-0.5*(nu_vec+1));
                        ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma_inv,0.5*(nu_vec+1));
                        
                        loglik_RC(j) = -ln_PWD_Sigma  + ln_PWD_j_yy_Sigma;
                    end
                    
                    loglik_RC = loglik_RC - (k/2)*log(pi) + LogLik.lnmultigamma_nu_vec(k, (nu_vec+1)/2 )   -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
                    LLF       = -sum(loglik_RC);
                    
                    if isnan(LLF)
                        LLF = 1e12;
                    end
                    
                end
            end
        end
        
        function [LLF,loglik_r2] = LogLik_T(k,T,params,y_mat);
            
            nu_t = params(1);
            
            Sigma = cov(y_mat);
            
            loglik_r2= zeros(T,1);
            inv_Sigma = inv(Sigma);
            logdet_Sigma = log(det(Sigma));
            
            
            for j =1 :T
                %loglik_r2(j) = - 1/2 * logdet_Delta  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(inv_Delta*r2(:,:,j)));
                loglik_r2(j) = - 1/2 * logdet_Sigma  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*y_mat(j,:)*inv_Sigma*y_mat(j,:)');
                
            end
            
            loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);
            LLF = -sum(loglik_r2);
            
            if isnan(LLF)
                LLF = 1e12;
                
            end
            % end
            
            
            
        end
        
        
        
        function [LLF,loglik_r2,V] = LogLik_BEKK_T(k,T,params,y_mat,ind_BEKK);
            
            
            
            
            if ind_BEKK==1
                para_BEKK = params(1:2);
                [V,Flag] = Filter_New.BEKK_CT(k,T,para_BEKK, y_mat);
            elseif ind_BEKK==2
                para_BEKK = params(1:k+1);
                [V,Flag] = Filter_New.BEKK_DIAG_CT(k,T,para_BEKK, y_mat);
            end
            
            
            
            
            nu_t = params(end);
            
            if Flag ==1
                LLF = 1e14;
                loglik_r2= NaN;
            else
                
                loglik_r2= zeros(T,1);
                
                
                for j =1 :T
                    
                    inv_Sigma_t = V(:,:,j)\eye(k);
                    logdet_Sigma_t = log(det(V(:,:,j)));
                    %loglik_r2(j) = - 1/2 * logdet_Delta  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(inv_Delta*r2(:,:,j)));
                    loglik_r2(j) = - 1/2 * logdet_Sigma_t  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*y_mat(j,:)*inv_Sigma_t*y_mat(j,:)');
                    
                end
                
                loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);
                LLF = -sum(loglik_r2);
                
                if isnan(LLF)
                    LLF = 1e12;
                    
                end
            end
            
            
            
        end
        
        
        function [LLF,loglik_r2] = LogLik_Mult_T_c_given_V_DAY(k,T,params,y_mat,V_DAY_vech,ind_c);
            
            
            if ind_c ==1
                c_hat = params(1);
            elseif ind_c==2
                c_vec = params(1:k);
                c_hat = c_vec*c_vec';
            end
            
            nu_t = params(end);
            loglik_r2= zeros(T,1);
            
            for j =1 :T
                V_DAY_j = Admin.Vech2Sym(k,V_DAY_vech(j,:));
                V_CTC_t  = c_hat.*V_DAY_j;
                inv_Sigma_t = V_CTC_t\eye(k);
                logdet_Sigma_t = log(det(V_CTC_t));
                %loglik_r2(j) = - 1/2 * logdet_Delta  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(inv_Delta*r2(:,:,j)));
                loglik_r2(j) = - 1/2 * logdet_Sigma_t  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*y_mat(j,:)*inv_Sigma_t*y_mat(j,:)');
                
            end
            
            loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);
            LLF = -sum(loglik_r2);
            
            if isnan(LLF)
                LLF = 1e12;
                
            end
            
            
            
            
        end
        
        
        
        function [LLF,loglik_r2] = LogLik_Mult_TRiesz_c_given_V_DAY(k,T,params,y_mat,V_DAY_vech);
            
            % let op: we doen c * V_t, wat uiteindeljk ook c * Sigma_t is
            % (schaling triesz)
            
            c = params(1);
            nu_vec = params(2:end);
            loglik_r2= zeros(T,1);
            
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            for j = 1:T
                V_DAY_t = Admin.Vech2Sym(k,V_DAY_vech(j,:));
                
                V_CTC_t  = c* V_DAY_t;
                Sigma_inv_j_0    = V_CTC_t\eye(k);
                
                [L_0_trans,ind_chol] = chol(Sigma_inv_j_0);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                inv_Sigma_j = L_trans' * L_trans;
                
                [~,ind_chol] = chol(inv_Sigma_j);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                
                j_Sigma = inv(inv_Sigma_j);
                j_yy       = y_mat(j,:)'*y_mat(j,:);
                j_yy_Sigma = j_yy + j_Sigma;
                
                
                if rcond(j_yy_Sigma)<1e-15
                    loglik_r2= NaN;
                    break
                end
                
                ln_PWD_Sigma = LogLik.ln_PWD(inv_Sigma_j,0.5*nu_vec);
                j_yy_Sigma_inv = j_yy_Sigma\eye(k);
                
                [~,ind_chol] = chol(j_yy_Sigma_inv);
                if ind_chol>0
                    loglik_r2= NaN;
                    break
                end
                
                ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma_inv,0.5*(nu_vec+1));
                loglik_r2(j) = -ln_PWD_Sigma  + ln_PWD_j_yy_Sigma;
                
                
            end
            
            loglik_r2 = loglik_r2 - (k/2)*log(pi) + LogLik.lnmultigamma_nu_vec(k, (nu_vec+1)/2 )   -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
            LLF       = -sum(loglik_r2);
            
            if isnan(LLF)
                LLF = 1e12;
            end
            
        end
        
        
        function [LLF,loglik_r2,c_vec, V_CTC_mat,s_vec] = LogLik_Mult_TRiesz_tv_c_given_V_DAY(k,T,params,y_mat,V_DAY_vech,c_1);
            
            % let op: we doen c * V_t, wat uiteindeljk ook c * Sigma_t is
            % (schaling triesz)
            
            c_vec = zeros(T,1);
            c_vec(1) = c_1;
            
            omega = params(1);
            alfa = params(2);
            beta = params(3);
            
            %c_vec(1) = omega/(1-beta);
            
            nu_vec = params(4:end);
            V_CTC_mat = zeros(k,k,T);
            loglik_r2= zeros(T,1);
            s_vec = zeros(T,1);
            
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            K_k = Admin.commutation(k,k);
            D_k = Admin.duplication(k);
            D_tilde_k = Admin.duplication_gen(k);
            I_k = eye(k);
            I_k_vech = eye(k*(k+1)/2);
            
            inv_D_prime_D_D_prime    = inv(D_k'*D_k) *D_k';
            %temp_score_part = temp_score_D * (eye(k^2) + K_k);
            
            K_k_D_tilde = K_k*D_tilde_k;
            S_D_tilde = Admin.Compute_S_d_matrix(1:k)';
            dau_diag_U_dau_vech_U_prime = S_D_tilde*K_k_D_tilde;
            
            for j = 1:T
                V_DAY_t = Admin.Vech2Sym(k,V_DAY_vech(j,:));
                %Sigma_inv_j_0    = V_DAY_t\eye(k);
                
                V_CTC_t  = c_vec(j)* V_DAY_t;
                V_CTC_mat(:,:,j) =  V_CTC_t ;
                Sigma_inv_j_0    = V_CTC_t\eye(k);
                
                [L_0_trans,ind_chol] = chol(Sigma_inv_j_0);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                inv_Sigma_j = L_trans' * L_trans;
                
                [~,ind_chol] = chol(inv_Sigma_j);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                %j_Sigma_DAY = inv(inv_Sigma_j);
                %j_Sigma_CTC = c_vec(j)* j_Sigma_DAY;
                %inv_Sigma_j = inv(j_Sigma_CTC);
                
                %[~,ind_chol] = chol(inv_Sigma_j);
                %if ind_chol>0
                %    loglik_r2 = NaN;
                %    break;
                %end
                
                
                j_Sigma = inv(inv_Sigma_j);
                j_yy       = y_mat(j,:)'*y_mat(j,:);
                j_yy_Sigma = j_yy + j_Sigma;
                
                %j_yy_Sigma = j_yy + j_Sigma_CTC;
                
                if rcond(j_yy_Sigma)<1e-15
                    loglik_r2= NaN;
                    break
                end
                
                ln_PWD_Sigma = LogLik.ln_PWD(inv_Sigma_j,0.5*nu_vec);
                j_yy_Sigma_inv = j_yy_Sigma\eye(k);
                
                [~,ind_chol] = chol(j_yy_Sigma_inv);
                if ind_chol>0
                    loglik_r2= NaN;
                    break
                end
                
                
                if j<T
                    L_inv_S  = chol(j_yy_Sigma_inv);
                    U_t      = inv(L_inv_S);
                    diag_U_t_vec = diag(U_t);
                    
                    s_t_1 = 0.5*sum(nu_vec)/c_vec(j);
                    
                    dau_ln_PWD_dau_diag_U = -(nu_vec'+1)./diag_U_t_vec';
                    
                    %s_t_2_2 =  S_D_tilde;
                    %s_t_2_3 =  K_k_D_tilde;
                    %s_t_2_6 = j_Sigma_DAY(:);
                    
                    Z = kron(I_k,U_t) + kron(U_t,I_k)*K_k;
                    temp_score_2                 = inv_D_prime_D_D_prime * Z*D_tilde_k;
                    dau_vech_U_prime_dau_vech_S  = temp_score_2\I_k_vech;
                    dau_vech_S_dau_vec_S         = inv_D_prime_D_D_prime;
                    dau_vec_S_dau_c              = j_Sigma(:)/c_vec(j);
                    
                    s_t = c_vec(j)^2 * (s_t_1 + dau_ln_PWD_dau_diag_U*dau_diag_U_dau_vech_U_prime*dau_vech_U_prime_dau_vech_S...
                        *dau_vech_S_dau_vec_S*dau_vec_S_dau_c);
                    
                    s_vec(j) = s_t;
                    
                    c_vec(j+1) = omega + alfa * s_t + beta * c_vec(j);
                end
                
                
                ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma_inv,0.5*(nu_vec+1));
                loglik_r2(j) = -ln_PWD_Sigma  + ln_PWD_j_yy_Sigma;
                
                
            end
            
            loglik_r2 = loglik_r2 - (k/2)*log(pi) + LogLik.lnmultigamma_nu_vec(k, (nu_vec+1)/2 )   -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
            LLF       = -sum(loglik_r2);
            
            if isnan(LLF)
                LLF = 1e12;
            end
            
        end
        
        
        
        function [LLF,loglik_r2,c_mat, V_CTC_mat,s_mat] = LogLik_Mult_TRiesz_tv_c_vec_given_V_DAY(k,T,params,y_mat,V_DAY_vech,c_1_vec);
            
            % let op: we doen c * V_t, wat uiteindeljk ook c * Sigma_t is
            % (schaling triesz)
            
            c_mat = zeros(T,k);
            c_mat(1,:) = c_1_vec';
            
            omega = params(1:k)';
            %omega = params(1);
            alfa = params(k+1)';
            beta = params(k+2);
            
            %c_vec(1) = omega/(1-beta);
            
            nu_vec = params(k+3:end);
            V_CTC_mat = zeros(k,k,T);
            loglik_r2= zeros(T,1);
            s_mat = zeros(T,k);
            
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            K_k = Admin.commutation(k,k);
            D_k = Admin.duplication(k);
            D_tilde_k = Admin.duplication_gen(k);
            I_k = eye(k);
            I_k_vech = eye(k*(k+1)/2);
            
            inv_D_prime_D_D_prime    = inv(D_k'*D_k) *D_k';
            %temp_score_part = temp_score_D * (eye(k^2) + K_k);
            
            K_k_D_tilde = K_k*D_tilde_k;
            S_D_tilde = Admin.Compute_S_d_matrix(1:k)';
            dau_diag_U_dau_vech_U_prime = S_D_tilde*K_k_D_tilde;
            dau_vec_C_dau_diag_C = S_D_tilde';
            
            for j = 1:T
                V_DAY_t = Admin.Vech2Sym(k,V_DAY_vech(j,:));
                Sigma_inv_j_0    = V_DAY_t\I_k;
                
                C_mat = diag(c_mat(j,:));
                
                %V_CTC_t  = c_mat(j,:)* V_DAY_t;
                %V_CTC_mat(:,:,j) =  V_CTC_t ;
                %Sigma_inv_j_0    = V_CTC_t\eye(k);
                
                [L_0_trans,ind_chol] = chol(Sigma_inv_j_0);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                inv_Sigma_j = L_trans' * L_trans;
                
                [~,ind_chol] = chol(inv_Sigma_j);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                j_Sigma_DAY = inv_Sigma_j\I_k;
                j_Sigma_CTC = C_mat* j_Sigma_DAY * C_mat;
                
                if rcond(j_Sigma_CTC)<1e-15
                    loglik_r2 = NaN;
                    break;
                end
                
                
                inv_Sigma_j = inv(j_Sigma_CTC);
                
                [~,ind_chol] = chol(inv_Sigma_j);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                
                % j_Sigma = inv(inv_Sigma_j);
                j_yy       = y_mat(j,:)'*y_mat(j,:);
                %j_yy_Sigma = j_yy + j_Sigma;
                j_yy_Sigma = j_yy + j_Sigma_CTC;
                
                if rcond(j_yy_Sigma)<1e-15
                    loglik_r2 = NaN;
                    break;
                end
                
                ln_PWD_Sigma = LogLik.ln_PWD(inv_Sigma_j,0.5*nu_vec);
                
                j_yy_Sigma_inv = j_yy_Sigma\eye(k);
                
                [~,ind_chol] = chol(j_yy_Sigma_inv);
                if ind_chol>0
                    loglik_r2= NaN;
                    break
                end
                
                
                if j<T
                    L_inv_S  = chol(j_yy_Sigma_inv);
                    U_t      = inv(L_inv_S);
                    diag_U_t_vec = diag(U_t);
                    
                    s_t_1 = nu_vec'./c_mat(j,:);
                    
                    dau_ln_PWD_dau_diag_U = -(nu_vec'+1)./diag_U_t_vec';
                    
                    
                    Z = kron(I_k,U_t) + kron(U_t,I_k)*K_k;
                    temp_score_2                 = inv_D_prime_D_D_prime * Z*D_tilde_k;
                    dau_vech_U_prime_dau_vech_S  = temp_score_2\I_k_vech;
                    dau_vech_S_dau_vec_S         = inv_D_prime_D_D_prime;
                    
                    
                    C_Sigma_j = C_mat * j_Sigma_DAY;
                    dau_vec_S_dau_vec_C = kron(C_Sigma_j,I_k) + kron(I_k,C_Sigma_j);
                    
                    s_t = c_mat(j,:).^2.* (s_t_1 + dau_ln_PWD_dau_diag_U*dau_diag_U_dau_vech_U_prime*dau_vech_U_prime_dau_vech_S...
                        *dau_vech_S_dau_vec_S*dau_vec_S_dau_vec_C*dau_vec_C_dau_diag_C);
                    
                    s_mat(j,:) = s_t;
                    
                    c_mat(j+1,:) = omega + alfa.* s_t + beta * c_mat(j,:);
                end
                
                
                ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma_inv,0.5*(nu_vec+1));
                loglik_r2(j,1) = -ln_PWD_Sigma  + ln_PWD_j_yy_Sigma;
                
                
            end
            
            loglik_r2 = loglik_r2 - (k/2)*log(pi) + LogLik.lnmultigamma_nu_vec(k, (nu_vec+1)/2 )   -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
            LLF       = -sum(loglik_r2);
            
            if isnan(LLF)
                LLF = 1e12;
            end
            
        end
        
        
        function [LLF,loglik_r2,c_vec, V_CTC_mat,s_vec] = LogLik_Mult_T_tv_c_given_V_DAY(k,T,params,y_mat,V_DAY_vech,c_1);
            
            
            c_vec = zeros(T,1);
            c_vec(1) = c_1;
            
            omega = params(1);
            alfa = params(2);
            beta = params(3);
            
            nu_t = params(end);
            V_CTC_mat = zeros(k,k,T);
            loglik_r2= zeros(T,1);
            s_vec = zeros(T,1);
            
            for j =1 :T
                
                V_DAY_t = Admin.Vech2Sym(k,V_DAY_vech(j,:));
                V_CTC_t  = c_vec(j)*V_DAY_t;
                
                V_CTC_mat(:,:,j) =  V_CTC_t ;
                
                if isnan(rcond(V_CTC_t))
                    loglik_r2 = NaN;
                    break
                end
                
                inv_Sigma_t = V_CTC_t\eye(k);
                y_V_inv_y = y_mat(j,:)*inv_Sigma_t*y_mat(j,:)';
                
                if j<T
                    w_t = (nu_t+k)/(nu_t-2 + y_V_inv_y);
                    nabla_t = 0.5 * ((w_t/c_vec(j)) *  y_V_inv_y - k/c_vec(j));
                    s_t = c_vec(j)^2 * nabla_t;
                    s_vec(j) = s_t;
                    c_vec(j+1) = omega + alfa * s_t + beta * c_vec(j);
                end
                
                logdet_Sigma_t = log(det(V_CTC_t));
                %loglik_r2(j) = - 1/2 * logdet_Delta  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(inv_Delta*r2(:,:,j)));
                loglik_r2(j) = - 1/2 * logdet_Sigma_t  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*y_V_inv_y);
                
            end
            
            loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);
            LLF = -sum(loglik_r2);
            
            if isnan(LLF)
                LLF = 1e12;
            end
            
            
            
            
        end
        
        
        
        function [LLF,loglik_r2,c_mat, V_CTC_mat] = LogLik_Mult_T_tv_c_vec_given_V_DAY(k,T,params,y_mat,V_DAY_vech,c_1_vec);
            
            
            c_mat = zeros(T,k);
            c_mat(1,:) = c_1_vec;
            
            omega_vec = params(1:k)';
            alfa = params(k+1);
            beta = params(k+2);
            
            nu_t = params(end);
            V_CTC_mat = zeros(k,k,T);
            loglik_r2= zeros(T,1);
            
            %S_d_mat = Admin.Compute_S_d_matrix(1:k);
            
            for j =1 :T
                
                V_DAY_j = Admin.Vech2Sym(k,V_DAY_vech(j,:));
                
                c_vec_j = c_mat(j,:);
                C_t = diag(c_vec_j);
                V_CTC_t  = C_t *V_DAY_j* C_t;
                
                V_CTC_mat(:,:,j) =  V_CTC_t ;
                
                if isnan(rcond(V_CTC_t)) || rcond(V_CTC_t)<1e-15
                    loglik_r2 = NaN;
                    break
                end
                
                inv_Sigma_t = V_CTC_t\eye(k);
                y_V_inv_y = y_mat(j,:)*inv_Sigma_t*y_mat(j,:)';
                
                
                
                if j<T
                    w_t = (nu_t+k)/(nu_t-2 + y_V_inv_y);
                    
                    
                    inv_C_t = diag(1./ c_vec_j);
                    A_t = inv_C_t *  (y_mat(j,:)'*y_mat(j,:))*inv_Sigma_t;
                    A_t_trans = A_t';
                    
                    % wegschalen met C_t ervoor en erna
                    % geeft zelfde effect als nable * c_vec.^2
                    
                    %A_t = y_mat(j,:)'*y_mat(j,:)*inv_Sigma_t*C_t;
                    %nabla_t = w_t *(A_t(:)' + A_t_trans(:)') *S_d_mat - c_vec_j;
                    
                    %A_t_trans = A_t';
                    
                    %nabla_t = w_t *(A_t(:)' + A_t_trans(:)') *S_d_mat - 1./c_vec_j;
                    % truc lucas: als je weet wat S_d is, dan geldt
                    % onderstaande:
                    
                    nabla_t   = 0.5 * w_t * diag(A_t + A_t_trans)' - 1./c_vec_j;
                    
                    s_t = c_mat(j,:).^2 .* nabla_t;
                    %s_t =  nabla_t;
                    c_mat(j+1,:) = omega_vec + alfa * s_t + beta * c_mat(j,:);
                end
                
                logdet_Sigma_t = log(det(V_CTC_t));
                %loglik_r2(j) = - 1/2 * logdet_Delta  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(inv_Delta*r2(:,:,j)));
                loglik_r2(j) = - 1/2 * logdet_Sigma_t  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*y_V_inv_y);
                
            end
            
            loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);
            LLF = -sum(loglik_r2);
            
            if isnan(LLF)
                LLF = 1e12;
            end
            
            
            
            
        end
        
        
        function [LLF,loglik_r2,V] = LogLik_MV_HEAVY_T(k,T,T_w,params,y_mat,RC_mat,ind_HEAVY);
            
            if ind_HEAVY==1
                para_HEAVY = params(1:2);
                [V,Flag] = Filter_New_MV.HEAVY_CT(k,T,para_HEAVY, y_mat,RC_mat);
            elseif ind_HEAVY==2
                para_HEAVY = params(1:k+1);
                [V,Flag] = Filter_New_MV.HEAVY_DIAG_CT(k,T,para_HEAVY, y_mat,RC_mat);
            end
            
            nu_t = params(end);
            
            if Flag ==1
                LLF = 1e14;
                loglik_r2= NaN;
            else
                
                loglik_r2= zeros(T,1);
                
                
                for j =1 :T
                    
                    inv_Sigma_t = V(:,:,j)\eye(k);
                    logdet_Sigma_t = log(det(V(:,:,j)));
                    %loglik_r2(j) = - 1/2 * logdet_Delta  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(inv_Delta*r2(:,:,j)));
                    loglik_r2(j) = - 1/2 * logdet_Sigma_t  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*y_mat(j,:)*inv_Sigma_t*y_mat(j,:)');
                    
                end
                
                loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);
                LLF = -sum(loglik_r2(T_w+1:end));
                
                if isnan(LLF)
                    LLF = 1e12;
                    
                end
            end
            
            
            
        end
        
        
        
        
        function [LLF,loglik_r2,V] = LogLik_BEKK_N(k,T,params,y_mat,ind_BEKK);
            
            
            
            if ind_BEKK==1
                para_BEKK = params(1:2);
                [V,Flag] = Filter_New_MV.BEKK_CT(k,T,para_BEKK, y_mat);
            elseif ind_BEKK==2
                para_BEKK = params(1:k+1);
                [V,Flag] = Filter_New_MV.BEKK_DIAG_CT(k,T,para_BEKK, y_mat);
            end
            
            
            
            
            if Flag ==1
                LLF = 1e14;
                loglik_r2= NaN;
            else
                
                loglik_r2= zeros(T,1);
                for j =1 :T
                    
                    inv_Sigma_t = V(:,:,j)\eye(k);
                    logdet_Sigma_t = log(det(V(:,:,j)));
                    %loglik_r2(j) = - 1/2 * logdet_Delta  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(inv_Delta*r2(:,:,j)));
                    loglik_r2(j) = - 1/2 * logdet_Sigma_t  - 0.5* y_mat(j,:)*inv_Sigma_t*y_mat(j,:)';
                    
                end
                
                loglik_r2 = loglik_r2  - 0.5*k*log(2*pi);
                LLF = -sum(loglik_r2);
                
                if isnan(LLF)
                    LLF = 1e12;
                    
                end
            end
            
            
            
        end
        
        
        function [LLF,loglik_r2] = LogLik_Mult_N_c_given_V_DAY(k,T,params,y_mat,V_DAY_vech,ind_c);
            
            if ind_c==1
                c_hat = params;
            elseif ind_c==2
                c_vec = params(1:k);
                c_hat =  c_vec* c_vec';
            end
            
            loglik_r2= zeros(T,1);
            for j =1 :T
                V_DAY_t = Admin.Vech2Sym(k,V_DAY_vech(j,:));
                V_CTC_t = c_hat.*V_DAY_t;
                inv_Sigma_t = V_CTC_t\eye(k);
                logdet_Sigma_t = log(det(V_CTC_t));
                %loglik_r2(j) = - 1/2 * logdet_Delta  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(inv_Delta*r2(:,:,j)));
                loglik_r2(j) = - 1/2 * logdet_Sigma_t  - 0.5* y_mat(j,:)*inv_Sigma_t*y_mat(j,:)';
                
            end
            
            loglik_r2 = loglik_r2  - 0.5*k*log(2*pi);
            LLF = -sum(loglik_r2);
            
            if isnan(LLF)
                LLF = 1e12;
                
            end
            
            
            
        end
        
        
        
        function [LLF,loglik_r2,V] = LogLik_MV_HEAVY_N(k,T,T_w,params,y_mat,RC_mat,ind_HEAVY);
            
            
            
            if ind_HEAVY==1
                para_HEAVY = params(1:2);
                [V,Flag] = Filter_New_MV.HEAVY_CT(k,T,para_HEAVY, y_mat,RC_mat);
            elseif ind_HEAVY==2
                para_HEAVY = params(1:k+1);
                [V,Flag] = Filter_New_MV.HEAVY_DIAG_CT(k,T,para_HEAVY,y_mat,RC_mat);
            end
            
            if Flag ==1
                LLF = 1e14;
                loglik_r2= NaN;
            else
                
                loglik_r2= zeros(T,1);
                for j =1 :T
                    
                    inv_Sigma_t = V(:,:,j)\eye(k);
                    logdet_Sigma_t = log(det(V(:,:,j)));
                    loglik_r2(j) = - 1/2 * logdet_Sigma_t  - 0.5* y_mat(j,:)*inv_Sigma_t*y_mat(j,:)';
                    
                end
                
                loglik_r2 = loglik_r2  - 0.5*k*log(2*pi);
                LLF = -sum(loglik_r2(T_w+1:end));
                
                if isnan(LLF)
                    LLF = 1e12;
                    
                end
            end
            
            
            
        end
        
        
        function [LLF,loglik_r2,H_TOT_mat,H_ST_mat] = LogLik_MV_HEAVY_LT_spec_N(k,T,T_w,params,return_mat,RC_mat);
            
            
            para_HEAVY_ST = params(1:2);
            I_k = eye(k);
            
            %H_bar = mean(RC_mat,3);
            H_bar = cov(return_mat);
            chol_H_LT_bar = chol(H_bar)';
            inv_chol_H_LT_bar = chol_H_LT_bar\I_k;
            
            
            H_ST_mat = zeros(k,k,T);
            H_TOT_mat = zeros(k,k,T);
            H_ST_mat(:,:,1) = eye(k);
            Omega_ST = (1-sum(para_HEAVY_ST))*I_k;
            loglik_r2= zeros(T,1);
            Flag = 0;
            
            for j =1 :T
                
                % compute H_tot
                H_TOT_mat(:,:,j) = chol_H_LT_bar  * H_ST_mat(:,:,j) * chol_H_LT_bar';
                
                % compute V_t
                
                [~,ind_pd_H_ST] = chol(H_ST_mat(:,:,j));
                if ind_pd_H_ST>0
                    Flag = 1;
                    loglik_r2 = NaN;
                    break
                end
                
                if j<T
                    
                    % update H_ST
                    H_ST_mat(:,:,j+1) = Omega_ST + para_HEAVY_ST(1)*(inv_chol_H_LT_bar*RC_mat(:,:,j))*inv_chol_H_LT_bar' + para_HEAVY_ST(2) * H_ST_mat(:,:,j);
                    
                    inv_Sigma_t = H_TOT_mat(:,:,j)\eye(k);
                    logdet_Sigma_t = log(det(H_TOT_mat(:,:,j)));
                    loglik_r2(j) = - 1/2 * logdet_Sigma_t  - 0.5* return_mat(j,:)*inv_Sigma_t*return_mat(j,:)';
                    
                end
                
                
            end
            
            if Flag==1
                LLF = 1e12;
            else
                
                loglik_r2 = loglik_r2  - 0.5*k*log(2*pi);
                LLF = -sum(loglik_r2(T_w+1:end));
                
                if isnan(LLF)
                    LLF = 1e12;
                end
            end
            
            
        end
        
        function [LLF,loglik_r2,H_TOT_mat,H_ST_mat] = LogLik_MV_HEAVY_LT_spec_T(k,T,T_w,params,return_mat,RC_mat);
            
            
            para_HEAVY_ST = params(1:2);
            nu_t = params(3);
            
            I_k = eye(k);
            
            H_bar = mean(RC_mat,3);
            chol_H_LT_bar = chol(H_bar)';
            inv_chol_H_LT_bar = chol_H_LT_bar\I_k;
            
            
            H_ST_mat = zeros(k,k,T);
            H_TOT_mat = zeros(k,k,T);
            H_ST_mat(:,:,1) = eye(k);
            Omega_ST = (1-sum(para_HEAVY_ST))*I_k;
            loglik_r2= zeros(T,1);
            Flag = 0;
            
            for j =1 :T
                
                % compute H_tot
                H_TOT_mat(:,:,j) = chol_H_LT_bar  * H_ST_mat(:,:,j) * chol_H_LT_bar';
                
                % compute V_t
                
                [~,ind_pd_H_ST] = chol(H_ST_mat(:,:,j));
                if ind_pd_H_ST>0
                    Flag = 1;
                    loglik_r2 = NaN;
                    break
                end
                
                if j<T
                    
                    % update H_ST
                    H_ST_mat(:,:,j+1) = Omega_ST + para_HEAVY_ST(1)*(inv_chol_H_LT_bar*RC_mat(:,:,j))*inv_chol_H_LT_bar' + para_HEAVY_ST(2) * H_ST_mat(:,:,j);
                    
                    inv_Sigma_t = H_TOT_mat(:,:,j)\eye(k);
                    logdet_Sigma_t = log(det(H_TOT_mat(:,:,j)));
                    
                    loglik_r2(j,1) = - 1/2 * logdet_Sigma_t  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*return_mat(j,:)*inv_Sigma_t*return_mat(j,:)');
                    
                end
                
            end
            
            
            if Flag==1
                LLF = 1e12;
            else
                
                loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);
                LLF = -sum(loglik_r2(T_w+1:end));
                
                if isnan(LLF)
                    LLF = 1e12;
                end
            end
            
            
        end
        
        
        function [LLF,loglik_r2,H_TOT_mat] = LogLik_BEKK_LT_spec_N(k,T,T_w,params,return_mat);
            
            
            para_BEKK_ST = params(1:2);
            I_k = eye(k);
            
            H_bar = cov(return_mat);
            chol_H_LT_bar = chol(H_bar)';
            inv_chol_H_LT_bar = chol_H_LT_bar\I_k;
            
            
            H_ST_mat = zeros(k,k,T);
            H_TOT_mat = zeros(k,k,T);
            H_ST_mat(:,:,1) = eye(k);
            Omega_ST = (1-sum(para_BEKK_ST))*I_k;
            loglik_r2= zeros(T,1);
            Flag = 0;
            
            for j =1 :T
                
                r_t_r_t_prime = return_mat(j,:)'*return_mat(j,:);
                
                % compute H_tot
                H_TOT_mat(:,:,j) = chol_H_LT_bar  * H_ST_mat(:,:,j) * chol_H_LT_bar';
                
                % compute V_t
                
                [~,ind_pd_H_ST] = chol(H_ST_mat(:,:,j));
                if ind_pd_H_ST>0
                    Flag = 1;
                    loglik_r2 = NaN;
                    break
                end
                
                if j<T
                    
                    % update H_ST
                    H_ST_mat(:,:,j+1) = Omega_ST + para_BEKK_ST(1)*(inv_chol_H_LT_bar*r_t_r_t_prime)*inv_chol_H_LT_bar' + para_BEKK_ST(2) * H_ST_mat(:,:,j);
                    
                    inv_Sigma_t = H_TOT_mat(:,:,j)\eye(k);
                    logdet_Sigma_t = log(det(H_TOT_mat(:,:,j)));
                    loglik_r2(j) = - 1/2 * logdet_Sigma_t  - 0.5* return_mat(j,:)*inv_Sigma_t*return_mat(j,:)';
                    
                end
                
                
            end
            
            if Flag==1
                LLF = 1e12;
            else
                
                loglik_r2 = loglik_r2  - 0.5*k*log(2*pi);
                LLF = -sum(loglik_r2(T_w+1:end));
                
                if isnan(LLF)
                    LLF = 1e12;
                end
            end
            
            
        end
        
        
        
        
        
        function [LLF,loglik_r2,H_TOT_mat,H_ST_mat, Z_t_mat,V_t_mat] = LogLik_MV_MF2_HEAVY_N(k,T,T_w,params,return_mat,RC_mat);
            
            Flag = 0;
            
            chol_omega_LT_vech = [params(1) 0; params(2) params(3)];
            Omega_LT = chol_omega_LT_vech*chol_omega_LT_vech';
            
            [chol_test,ind_pd_H_ST] = chol( Omega_LT);
            if ind_pd_H_ST>0
                Flag = 1;
                loglik_r2 = NaN;
                
            else
                
                
                alpha_ST  = params(end-3);
                beta_ST  = params(end-2);
                alpha_LT  = params(end-1);
                beta_LT  = params(end);
                
                % [i_LT] = Admin.UnitLT(k);
                % A_mat_LT = para_HEAVY_LT(1)*eye(k) + para_HEAVY_LT(2)* i_LT;
                
                I_k = eye(k);
                
                H_LT_bar = mean(RC_mat,3);
                
                %Omega_LT = H_LT_bar  - A_mat_LT * H_LT_bar  * A_mat_LT'  - para_HEAVY_LT(end)*H_LT_bar;
                %Omega_LT = H_LT_bar * (1-alpha_LT - beta_LT);
                chol_H_LT_bar = chol(H_LT_bar)';
                inv_chol_H_LT_bar = chol_H_LT_bar\I_k;
                
                chol_Z_new = chol_H_LT_bar;
                inv_chol_Z_new = inv_chol_H_LT_bar;
                
                Z_t_mat = zeros(k,k,T);
                Z_t_mat(:,:,1) = H_LT_bar;
                V_t_mat = zeros(k,k,T);
                
                H_ST_mat = zeros(k,k,T);
                H_TOT_mat = zeros(k,k,T);
                H_ST_mat(:,:,1) = eye(k);
                Omega_ST = (1-alpha_ST-beta_ST)*I_k;
                loglik_r2= zeros(T,1);
                
                
                for j =1 :T
                    
                    % compute H_tot
                    H_TOT_mat_j= chol_Z_new * H_ST_mat(:,:,j) * chol_Z_new';
                    
                    if rcond(H_TOT_mat_j)<1e-12
                        Flag = 1;
                        loglik_r2 = NaN;
                        break
                    end
                    
                    H_TOT_mat(:,:,j) = H_TOT_mat_j;
                    
                    % compute V_t
                    
                    [chol_H_ST_mat_j,ind_pd_H_ST] = chol(H_ST_mat(:,:,j));
                    if ind_pd_H_ST>0
                        Flag = 1;
                        loglik_r2 = NaN;
                        break
                    end
                    
                    chol_H_ST_j = chol_H_ST_mat_j';
                    inv_chol_H_ST_j = chol_H_ST_j\eye(k);
                    V_t_mat(:,:,j) = inv_chol_H_ST_j * RC_mat(:,:,j) * inv_chol_H_ST_j';
                    
                    if j<T
                        
                        % update H_ST
                        H_ST_mat(:,:,j+1) = Omega_ST + alpha_ST*(inv_chol_Z_new*RC_mat(:,:,j))*inv_chol_Z_new' + beta_ST * H_ST_mat(:,:,j);
                        
                        % update Z_t
                        if j>=T_w
                            %Z_t_mat(:,:,j+1) = Omega_LT + A_mat_LT * mean(V_t_mat(:,:,j-T_w+1:j),3)*A_mat_LT' + para_HEAVY_LT(3)*Z_t_mat(:,:,j);
                            Z_t_mat(:,:,j+1) = Omega_LT + alpha_LT*mean(V_t_mat(:,:,j-T_w+1:j),3) + beta_LT*Z_t_mat(:,:,j);
                            Z_t_mat_new = Z_t_mat(:,:,j+1);
                        else
                            Z_t_mat(:,:,j+1) = H_LT_bar;
                            Z_t_mat_new = H_LT_bar;
                        end
                        
                        
                        [chol_Z_new_j,ind_pd_Z_j] = chol(Z_t_mat_new);
                        if ind_pd_Z_j>0
                            Flag = 1;
                            loglik_r2 = NaN;
                            break
                        end
                        
                        
                        chol_Z_new = chol_Z_new_j';
                        inv_chol_Z_new = chol_Z_new\eye(k);
                        
                        inv_Sigma_t = H_TOT_mat(:,:,j)\eye(k);
                        logdet_Sigma_t = log(det(H_TOT_mat(:,:,j)));
                        loglik_r2(j) = - 1/2 * logdet_Sigma_t  - 0.5* return_mat(j,:)*inv_Sigma_t*return_mat(j,:)';
                        
                    end
                    
                    
                end
                
            end
            
            if Flag==1
                LLF = 1e12;
            else
                
                loglik_r2 = loglik_r2  - 0.5*k*log(2*pi);
                LLF = -sum(loglik_r2(T_w+1:end));
                
                if isnan(LLF)
                    LLF = 1e12;
                end
            end
            
        end
        
        
        
        function [LLF,loglik_r2,H_TOT_mat,H_ST_mat, Z_t_mat,V_t_mat] = LogLik_MV_MF2_HEAVY_diag_CCC_N(k,T,T_w,params,return_mat,RC_mat);
            
            Flag = 0;
            
            
            
            omega_vec_LT = params(1:k)';            
                        
            alpha_ST  = params(end-3);
            beta_ST  = params(end-2);
            alpha_LT  = params(end-1);
            beta_LT  = params(end);
                               
                        
            RCOV_bar = mean(RC_mat,3);            
            RCORR_bar = corr(return_mat);
            
            Z_t_mat = zeros(T,k);
            Z_t_mat(1,:) = diag(RCOV_bar);
            V_t_mat = zeros(T,k);
            
            H_ST_mat = zeros(T,k);
            H_TOT_mat = zeros(k,k,T);
            H_ST_mat(1,:) = ones(1,k);
            Omega_vec_ST = (1-alpha_ST-beta_ST)*ones(1,k);
            loglik_r2= zeros(T,1);
            
            
            for j =1 :T
                
                Z_t_vec      = Z_t_mat(j,:);
                D_tau_mat_j  = diag(sqrt(Z_t_vec));
                %inv_D_tau_mat_j = diag(1./sqrt(Z_t_mat(j,:)));
                D_H_ST_mat_j = diag(sqrt(H_ST_mat(j,:)));
                
                % compute H_tot
                H_TOT_mat_j=  D_tau_mat_j *  D_H_ST_mat_j * RCORR_bar * D_H_ST_mat_j * D_tau_mat_j;
               
                if rcond(H_TOT_mat_j)<1e-12 || isnan(rcond(H_TOT_mat_j))
                    Flag = 1;
                    loglik_r2 = NaN;
                    break
                end
                
                H_TOT_mat(:,:,j) = H_TOT_mat_j;
                
                % compute V_t
                RV_vec_j = diag(RC_mat(:,:,j))';                             
                V_t_mat(j,:) = RV_vec_j./(H_ST_mat(j,:));
                
                if j<T
                    
                    % update H_ST
                    H_ST_mat(j+1,:) = Omega_vec_ST + alpha_ST*RV_vec_j./Z_t_vec + beta_ST * H_ST_mat(j,:);
                    
                    % update Z_t
                    if j>=T_w
                        %Z_t_mat(:,:,j+1) = Omega_LT + A_mat_LT * mean(V_t_mat(:,:,j-T_w+1:j),3)*A_mat_LT' + para_HEAVY_LT(3)*Z_t_mat(:,:,j);
                        Z_t_mat(j+1,:) = omega_vec_LT + alpha_LT*mean(V_t_mat(j-T_w+1:j,:)) + beta_LT*Z_t_mat(j,:);                    
                    else
                        Z_t_mat(j+1,:) = Z_t_mat(j,:);
                       
                    end
                                                                                                  
                    
                    inv_Sigma_t = H_TOT_mat(:,:,j)\eye(k);
                    logdet_Sigma_t = log(det(H_TOT_mat(:,:,j)));
                    loglik_r2(j) = - 1/2 * logdet_Sigma_t  - 0.5* return_mat(j,:)*inv_Sigma_t*return_mat(j,:)';
                    
                end
                
                
            end
            
            
            
            if Flag==1
                LLF = 1e12;
            else
                
                loglik_r2 = loglik_r2  - 0.5*k*log(2*pi);
                LLF = -sum(loglik_r2(T_w+1:end));
                
                if isnan(LLF)
                    LLF = 1e12;
                end
            end
            
        end
        
        
        function [LLF,h_TOT_mat,tau_mat,h_mat] = LogLik_MV_MF2_HEAVY_VOL_PART(k,T,T_w,params,return_mat,RV_mat)
            
            omega_LT = params(1:k)';
            alpha_ST = params(end-3);           
            beta_ST  = params(end-2);                        
            alpha_LT= params(end-1);
            beta_LT = params(end);
            
            
            eps_mat = return_mat;
            eps2_mat = eps_mat.^2;            
                                               
            h_mat = zeros(T,k);
            h_mat(1,:) = ones(1,k);
            
            tau_mat = zeros(T,k);
            v_mat = zeros(T,k);
            tau_mat(1:T_w,:) = ones(T_w,1)*var(eps_mat);
            
            for j = 1:T
                
                if j<T                                                      
                        v_mat(j,:) = RV_mat(j,:)./h_mat(j,:);                    
                        h_mat(j+1,:) =  (1-alpha_ST-beta_ST) + beta_ST * h_mat(j,:) + alpha_ST* RV_mat(j,:)./tau_mat(j,:);
                        
                        if j>=T_w
                            tau_mat(j+1,:) = omega_LT + alpha_LT*mean(v_mat(j-T_w+1:j,:)) + beta_LT * tau_mat(j,:);
                        end
                        
                end
                
            end
            
            h_TOT_mat= h_mat.*tau_mat;
            loglik_mat_r2 = -0.5*log(2*pi)*ones(T,k) -0.5 * log(h_TOT_mat) - 0.5*eps2_mat./h_TOT_mat;              
            LLF_vec= -sum(loglik_mat_r2(T_w+1:end,:));
            LLF = sum(LLF_vec);
            
            if isnan(LLF) || ~isreal(LLF) || isinf(abs(LLF))
                LLF = 1e14;
            end
            
        end                                       
        
        
        function [LLF,loglik_r2,H_TOT_mat,H_ST_mat,Z_t_mat,Rt_mat,V_t_mat] = LogLik_MV_MF2_HEAVY_diag_DCC_N(k,T,T_w,params,return_mat,RC_mat);
            
            Flag = 0;
                       
            omega_vec_LT = params(1:k)';            
                                    
            alpha_ST  = params(end-5);
            beta_ST  = params(end-4);
            alpha_LT  = params(end-3);
            beta_LT  = params(end-2);
            
            alpha_DCC = params(end-1);
            beta_DCC = params(end);
            
            RCOV_bar = mean(RC_mat,3);          
            Qt_current = corr(return_mat);
            Qbar_cDCC = Qt_current;
                        
            Z_t_mat = zeros(T,k);
            Z_t_mat(1,:) = diag(RCOV_bar);
            V_t_mat = zeros(T,k);
            
            H_ST_mat = zeros(T,k);
            H_TOT_mat = zeros(k,k,T);
            Rt_mat = zeros(k,k,T);
            Rt_mat(:,:,1) = Qt_current;
            H_ST_mat(1,:) = ones(1,k);
            Omega_vec_ST = (1-alpha_ST-beta_ST)*ones(1,k);
            loglik_r2= zeros(T,1);
            
            
            for j =1 :T
                
                Z_t_vec      = Z_t_mat(j,:);
                D_tau_mat_j  = diag(sqrt(Z_t_vec));
                %inv_D_tau_mat_j = diag(1./sqrt(Z_t_mat(j,:)));
                D_H_ST_mat_j = diag(sqrt(H_ST_mat(j,:)));
                
                u_vec_j = return_mat(j,:)./sqrt((Z_t_vec.*H_ST_mat(j,:)));
                
                Rt_DCC_j=  Qt_current./(sqrt(diag(Qt_current))*sqrt(diag(Qt_current))');
                Rt_mat(:,:,j) = Rt_DCC_j;
                
                
                % compute H_tot
                H_TOT_mat_j=  D_tau_mat_j *  D_H_ST_mat_j * Rt_DCC_j * D_H_ST_mat_j * D_tau_mat_j;
               
                if rcond(H_TOT_mat_j)<1e-12 || isnan(rcond(H_TOT_mat_j))
                    Flag = 1;
                    loglik_r2 = NaN;
                    break
                end
                
                H_TOT_mat(:,:,j) = H_TOT_mat_j;
                
                % compute V_t
                RV_vec_j = diag(RC_mat(:,:,j))';                             
                V_t_mat(j,:) = RV_vec_j./(H_ST_mat(j,:));
                
                if j<T                                        
                    
                    % update H_ST
                    H_ST_mat(j+1,:) = Omega_vec_ST + alpha_ST*RV_vec_j./Z_t_vec + beta_ST * H_ST_mat(j,:);
                    
                    % update Z_t
                    if j>=T_w
                        %Z_t_mat(:,:,j+1) = Omega_LT + A_mat_LT * mean(V_t_mat(:,:,j-T_w+1:j),3)*A_mat_LT' + para_HEAVY_LT(3)*Z_t_mat(:,:,j);
                        Z_t_mat(j+1,:) = omega_vec_LT + alpha_LT*mean(V_t_mat(j-T_w+1:j,:)) + beta_LT*Z_t_mat(j,:);                    
                    else
                        Z_t_mat(j+1,:) = Z_t_mat(j,:);
                       
                    end
                                                                             
               
                    temp = sqrt(diag(Qt_current));
                    Q_star_t = diag(temp(1:end));
                    Qt_new = Qbar_cDCC*(1 - alpha_DCC-beta_DCC) + alpha_DCC*(Q_star_t*u_vec_j')*u_vec_j*Q_star_t + beta_DCC*Qt_current;               
                    Qt_current = Qt_new;
                                                                                                                                                                             
                    
                    inv_Sigma_t = H_TOT_mat_j\eye(k);
                    logdet_Sigma_t = log(det(H_TOT_mat_j));
                    loglik_r2(j) = - 1/2 * logdet_Sigma_t  - 0.5* return_mat(j,:)*inv_Sigma_t*return_mat(j,:)';
                    
                end
                
                
            end
                                    
            if Flag==1
                LLF = 1e12;
            else
                
                loglik_r2 = loglik_r2  - 0.5*k*log(2*pi);
                LLF = -sum(loglik_r2(T_w+1:end));
                
                if isnan(LLF)
                    LLF = 1e12;
                end
            end
            
        end
        
        
        
        
        
        function [LLF,loglik_r2,H_TOT_mat,H_ST_mat, Z_t_mat,V_t_mat] = LogLik_MV_MF2_N(k,T,T_w,params,return_mat,ind_spec);
            
            
            
            para_HEAVY_ST = params(1:2);
            
            if ind_spec==1
                %alpha_LT = params(3);
                A_mat_LT = params(3);
                beta_LT  = params(end);
            elseif ind_spec==2
                alpha_LT_1  = params(3:3+k-1);
                beta_LT  = params(end);
                A_mat_LT = diag(alpha_LT_1);
            elseif ind_spec==3
                
                alpha_LT_1  = params(3:3+k-1);
                alpha_LT_2  = params(3+k);
                beta_LT  = params(end);
                
                [i_LT] = Admin.UnitLT(k);
                A_mat_LT = diag(alpha_LT_1) + alpha_LT_2 * i_LT;
                
            end
            
            
            I_k = eye(k);
            H_LT_bar = cov(return_mat);
            
            Omega_LT = H_LT_bar  - A_mat_LT * H_LT_bar  * A_mat_LT'  - beta_LT*H_LT_bar;
            %Omega_LT = H_LT_bar * (1-alpha_LT - beta_LT);
            chol_H_LT_bar = chol(H_LT_bar)';
            inv_chol_H_LT_bar = chol_H_LT_bar\I_k;
            
            chol_Z_new = chol_H_LT_bar;
            inv_chol_Z_new = inv_chol_H_LT_bar;
            
            Z_t_mat = zeros(k,k,T);
            Z_t_mat(:,:,1) = H_LT_bar;
            V_t_mat = zeros(k,k,T);
            
            H_ST_mat = zeros(k,k,T);
            H_TOT_mat = zeros(k,k,T);
            H_ST_mat(:,:,1) = eye(k);
            Omega_ST = (1-sum(para_HEAVY_ST))*I_k;
            loglik_r2= zeros(T,1);
            Flag = 0;
            
            for j =1 :T
                
                r_t_r_t_prime = return_mat(j,:)'*return_mat(j,:);
                
                % compute H_tot
                H_TOT_mat_j= chol_Z_new * H_ST_mat(:,:,j) * chol_Z_new';
                
                if rcond(H_TOT_mat_j)<1e-12
                    Flag = 1;
                    loglik_r2 = NaN;
                    break
                end
                
                H_TOT_mat(:,:,j) = H_TOT_mat_j;
                
                % compute V_t
                
                [chol_H_ST_mat_j,ind_pd_H_ST] = chol(H_ST_mat(:,:,j));
                if ind_pd_H_ST>0
                    Flag = 1;
                    loglik_r2 = NaN;
                    break
                end
                
                chol_H_ST_j = chol_H_ST_mat_j';
                inv_chol_H_ST_j = chol_H_ST_j\eye(k);
                V_t_mat(:,:,j) = (inv_chol_H_ST_j *r_t_r_t_prime)*inv_chol_H_ST_j';
                
                if j<T
                    
                    % update H_ST
                    H_ST_mat(:,:,j+1) = Omega_ST + para_HEAVY_ST(1)*(inv_chol_Z_new*r_t_r_t_prime)*inv_chol_Z_new' + para_HEAVY_ST(2) * H_ST_mat(:,:,j);
                    
                    % update Z_t
                    if j>=T_w
                        Z_t_mat(:,:,j+1) = Omega_LT + A_mat_LT * mean(V_t_mat(:,:,j-T_w+1:j),3)*A_mat_LT' + beta_LT*Z_t_mat(:,:,j);
                        %Z_t_mat(:,:,j+1) = Omega_LT + alpha_LT * mean(V_t_mat(:,:,j-T_w+1:j),3) + beta_LT*Z_t_mat(:,:,j);
                        %Z_t_mat(:,:,j+1) = Omega_LT + alpha_LT*mean(V_t_mat(:,:,j-T_w+1:j),3) + beta_LT*Z_t_mat(:,:,j);
                        Z_t_mat_new = Z_t_mat(:,:,j+1);
                    else
                        Z_t_mat(:,:,j+1) = H_LT_bar;
                        Z_t_mat_new = H_LT_bar;
                    end
                    
                    
                    [chol_Z_new_j,ind_pd_Z_j] = chol(Z_t_mat_new);
                    if ind_pd_Z_j>0
                        Flag = 1;
                        loglik_r2 = NaN;
                        break
                    end
                    
                    
                    chol_Z_new = chol_Z_new_j';
                    inv_chol_Z_new = chol_Z_new\eye(k);
                    
                    inv_Sigma_t = H_TOT_mat(:,:,j)\eye(k);
                    logdet_Sigma_t = log(det(H_TOT_mat(:,:,j)));
                    loglik_r2(j) = - 1/2 * logdet_Sigma_t  - 0.5* return_mat(j,:)*inv_Sigma_t*return_mat(j,:)';
                    
                end
                
                
            end
            
            if Flag==1
                LLF = 1e12;
            else
                
                loglik_r2 = loglik_r2  - 0.5*k*log(2*pi);
                LLF = -sum(loglik_r2(T_w+1:end));
                
                if isnan(LLF)
                    LLF = 1e12;
                end
            end
            
        end
        
        
        
        function [LLF,loglik_r2,V] = LogLik_BEKK_TRiesz_I_CT(k,T,params,y_mat,ind_BEKK);
            
            %ind BEKK = 1: scalar BEKK
            %ind BEKK = 2:  A eps eps A
            % CT_ means with covariance targeting
            
            
            
            
            if ind_BEKK==1
                para_BEKK = params(1:2);
                nu_vec = params(3:end);
                [V,Flag] = Filter_New.BEKK_CT(k,T,para_BEKK, y_mat);
            elseif ind_BEKK==2
                para_BEKK = params(1:2*k+2);
                nu_vec = params(2*k+3:end);
                [V,Flag] = Filter_New.BEKK_FULL_CT(k,T,para_BEKK, y_mat);
            elseif ind_BEKK==3
                para_BEKK = params(1:k+3);
                nu_vec = params(k+4:end);
                [V,Flag] = Filter_New.BEKK_RES_CT(k,T,para_BEKK, y_mat);
            elseif ind_BEKK==4
                para_BEKK = params(1:k+1);
                nu_vec = params(k+2:end);
                [V,Flag] = Filter_New.BEKK_DIAG_CT(k,T,para_BEKK, y_mat);
            elseif ind_BEKK==5
                para_BEKK = params(1:k+1);
                nu_vec = params(k+2:end);
                [V,Flag] = Filter_New.BEKK_SPILL_CT(k,T,para_BEKK, y_mat);
            end
            
            
            
            
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            
            
            if Flag ==1
                LLF = 1e14;
                loglik_r2= NaN;
            else
                loglik_r2 = zeros(T,1);
                
                for j = 1:T
                    
                    j_Sigma_0        = V(:,:,j);
                    Sigma_inv_j_0    = j_Sigma_0\eye(k);
                    
                    [L_0_trans,ind_chol] = chol(Sigma_inv_j_0);
                    if ind_chol>0
                        loglik_r2 = NaN;
                        break;
                    end
                    
                    L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                    inv_Sigma_j = L_trans' * L_trans;
                    
                    [~,ind_chol] = chol(inv_Sigma_j);
                    if ind_chol>0
                        loglik_r2 = NaN;
                        break;
                    end
                    
                    
                    if rcond(inv_Sigma_j)<=1e-12
                        loglik_r2 = NaN;
                        break;
                    end
                    
                    j_Sigma = inv(inv_Sigma_j);
                    
                    j_yy       = y_mat(j,:)'*y_mat(j,:);
                    j_yy_Sigma = j_yy + j_Sigma;
                    
                    ln_PWD_Sigma = LogLik.ln_PWD(inv_Sigma_j,0.5*nu_vec);
                    j_yy_Sigma_inv = j_yy_Sigma\eye(k);
                    
                    [~,ind_chol] = chol(j_yy_Sigma_inv);
                    if ind_chol>0
                        loglik_r2= NaN;
                        break
                    end
                    
                    %ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma,-0.5*(nu_vec+1));
                    ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma_inv,0.5*(nu_vec+1));
                    
                    loglik_r2(j) = -ln_PWD_Sigma  + ln_PWD_j_yy_Sigma;
                    
                    
                end
                
                loglik_r2 = loglik_r2 - (k/2)*log(pi) + LogLik.lnmultigamma_nu_vec(k, (nu_vec+1)/2 )   -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
                LLF       = -sum(loglik_r2);
                
                if isnan(LLF)
                    LLF = 1e12;
                end
                
            end
        end
        
        
        
        function [LLF,loglik_r2] = LogLik_BEKK_TRiesz_I_c_given_V_DAY(k,T,params,y_mat,V_DAY,ind_c);
            
            %ind BEKK = 1: scalar BEKK
            %ind BEKK = 2:  A eps eps A
            % CT_ means with covariance targeting
            
            if ind_c==1
                c_hat = params(1);
                nu_vec =params(2:end);
            elseif ind_c==2
                c_vec = params(1:k);
                nu_vec = params(k+1:end);
                c_hat = c_vec*c_vec';
                
            end
            
            
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            
            loglik_r2 = zeros(T,1);
            
            for j = 1:T
                
                j_Sigma_0        = c_hat.*V_DAY(:,:,j);
                Sigma_inv_j_0    = j_Sigma_0\eye(k);
                
                [L_0_trans,ind_chol] = chol(Sigma_inv_j_0);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                inv_Sigma_j = L_trans' * L_trans;
                
                [~,ind_chol] = chol(inv_Sigma_j);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                
                j_Sigma = inv(inv_Sigma_j);
                
                j_yy       = y_mat(j,:)'*y_mat(j,:);
                j_yy_Sigma = j_yy + j_Sigma;
                
                ln_PWD_Sigma = LogLik.ln_PWD(inv_Sigma_j,0.5*nu_vec);
                j_yy_Sigma_inv = j_yy_Sigma\eye(k);
                
                [~,ind_chol] = chol(j_yy_Sigma_inv);
                if ind_chol>0
                    loglik_r2= NaN;
                    break
                end
                
                %ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma,-0.5*(nu_vec+1));
                ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma_inv,0.5*(nu_vec+1));
                
                loglik_r2(j) = -ln_PWD_Sigma  + ln_PWD_j_yy_Sigma;
                
                
            end
            
            loglik_r2 = loglik_r2 - (k/2)*log(pi) + LogLik.lnmultigamma_nu_vec(k, (nu_vec+1)/2 )   -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
            LLF       = -sum(loglik_r2);
            
            if isnan(LLF)
                LLF = 1e12;
            end
            
        end
        
        
        
        function [LLF,loglik_r2,V] = LogLik_HEAVY_TRiesz_I_CT(k,T,params,y_mat,RC_mat,ind_HEAVY);
            
            %ind BEKK = 1: scalar BEKK
            %ind BEKK = 2:  A eps eps A
            % CT_ means with covariance targeting
            
            if ind_HEAVY==1
                para_HEAVY = params(1:2);
                nu_vec = params(3:end);
                [V,Flag] = Filter_New.HEAVY_CT(k,T,para_HEAVY, y_mat,RC_mat);
            elseif ind_HEAVY==2
                para_HEAVY = params(1:k+1);
                nu_vec = params(k+2:end);
                [V,Flag] = Filter_New.HEAVY_DIAG_CT(k,T,para_HEAVY, y_mat,RC_mat);
            end
            
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            if Flag ==1
                LLF = 1e14;
                loglik_r2= NaN;
            else
                loglik_r2 = zeros(T,1);
                
                for j = 1:T
                    
                    j_Sigma_0        = V(:,:,j);
                    Sigma_inv_j_0    = j_Sigma_0\eye(k);
                    
                    [L_0_trans,ind_chol] = chol(Sigma_inv_j_0);
                    if ind_chol>0
                        loglik_r2 = NaN;
                        break;
                    end
                    
                    L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                    inv_Sigma_j = L_trans' * L_trans;
                    
                    [~,ind_chol] = chol(inv_Sigma_j);
                    if ind_chol>0
                        loglik_r2 = NaN;
                        break;
                    end
                    
                    j_Sigma = inv(inv_Sigma_j);
                    
                    j_yy       = y_mat(j,:)'*y_mat(j,:);
                    j_yy_Sigma = j_yy + j_Sigma;
                    
                    ln_PWD_Sigma = LogLik.ln_PWD(inv_Sigma_j,0.5*nu_vec);
                    j_yy_Sigma_inv = j_yy_Sigma\eye(k);
                    
                    [~,ind_chol] = chol(j_yy_Sigma_inv);
                    if ind_chol>0
                        loglik_r2= NaN;
                        break
                    end
                    
                    ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma_inv,0.5*(nu_vec+1));
                    
                    loglik_r2(j) = -ln_PWD_Sigma  + ln_PWD_j_yy_Sigma;
                    
                    
                end
                
                loglik_r2 = loglik_r2 - (k/2)*log(pi) + LogLik.lnmultigamma_nu_vec(k, (nu_vec+1)/2 )   -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
                LLF       = -sum(loglik_r2);
                
                if isnan(LLF)
                    LLF = 1e12;
                end
                
            end
        end
        
        function [LLF,loglik_r2] = LogLik_BEKK_TRiesz_I_INDUSTRY_CT(k,T,params,y_mat,n_vec);
            
            % CT_ means with covariance targeting
            nr_groups = length(n_vec);
            para_BEKK = params(1:2);
            
            nu_vec = [];
            for g = 1:nr_groups
                nu_vec = [nu_vec; params(g+2)*ones(n_vec(g),1)];
            end
            
            %nu_vec = params(3:end);
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            [V,Flag] = Filter_New.BEKK_CT(k,T,para_BEKK, y_mat);
            
            if Flag ==1
                LLF = 1e14;
                loglik_r2= NaN;
            else
                loglik_r2 = zeros(T,1);
                
                for j = 1:T
                    
                    j_Sigma_0        = V(:,:,j);
                    Sigma_inv_j_0    = j_Sigma_0\eye(k);
                    
                    [L_0_trans,ind_chol] = chol(Sigma_inv_j_0);
                    if ind_chol>0
                        loglik_r2 = NaN;
                        break;
                    end
                    
                    L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                    inv_Sigma_j = L_trans' * L_trans;
                    
                    [~,ind_chol] = chol(inv_Sigma_j);
                    if ind_chol>0
                        loglik_r2 = NaN;
                        break;
                    end
                    
                    
                    j_Sigma = inv(inv_Sigma_j);
                    
                    j_yy       = y_mat(j,:)'*y_mat(j,:);
                    j_yy_Sigma = j_yy + j_Sigma;
                    
                    ln_PWD_Sigma = LogLik.ln_PWD(inv_Sigma_j,0.5*nu_vec);
                    j_yy_Sigma_inv = j_yy_Sigma\eye(k);
                    
                    [~,ind_chol] = chol(j_yy_Sigma_inv);
                    if ind_chol>0
                        loglik_r2= NaN;
                        break
                    end
                    
                    %ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma,-0.5*(nu_vec+1));
                    ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma_inv,0.5*(nu_vec+1));
                    
                    loglik_r2(j) = -ln_PWD_Sigma  + ln_PWD_j_yy_Sigma;
                    
                    
                end
                
                loglik_r2 = loglik_r2 - (k/2)*log(pi) + LogLik.lnmultigamma_nu_vec(k, (nu_vec+1)/2 )   -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
                LLF       = -sum(loglik_r2);
                
                if isnan(LLF)
                    LLF = 1e12;
                end
                
            end
        end
        
        
        function [LLF,loglik_r2] = LogLik_HEAVY_TRiesz_I_INDUSTRY_CT(k,T,params,y_mat,RC_mat,n_vec);
            
            % CT_ means with covariance targeting
            nr_groups = length(n_vec);
            para_HEAVY = params(1:2);
            
            nu_vec = [];
            for g = 1:nr_groups
                nu_vec = [nu_vec; params(g+2)*ones(n_vec(g),1)];
            end
            
            %nu_vec = params(3:end);
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            %[V,Flag] = Filter_New.BEKK_CT(k,T,para_HEAVY, y_mat);
            [V,Flag] = Filter_New.HEAVY_CT(k,T,para_HEAVY, y_mat,RC_mat);
            
            if Flag ==1
                LLF = 1e14;
                loglik_r2= NaN;
            else
                loglik_r2 = zeros(T,1);
                
                for j = 1:T
                    
                    j_Sigma_0        = V(:,:,j);
                    Sigma_inv_j_0    = j_Sigma_0\eye(k);
                    
                    [L_0_trans,ind_chol] = chol(Sigma_inv_j_0);
                    if ind_chol>0
                        loglik_r2 = NaN;
                        break;
                    end
                    
                    L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                    inv_Sigma_j = L_trans' * L_trans;
                    
                    [~,ind_chol] = chol(inv_Sigma_j);
                    if ind_chol>0
                        loglik_r2 = NaN;
                        break;
                    end
                    
                    
                    j_Sigma = inv(inv_Sigma_j);
                    
                    j_yy       = y_mat(j,:)'*y_mat(j,:);
                    j_yy_Sigma = j_yy + j_Sigma;
                    
                    ln_PWD_Sigma = LogLik.ln_PWD(inv_Sigma_j,0.5*nu_vec);
                    j_yy_Sigma_inv = j_yy_Sigma\eye(k);
                    
                    [~,ind_chol] = chol(j_yy_Sigma_inv);
                    if ind_chol>0
                        loglik_r2= NaN;
                        break
                    end
                    
                    %ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma,-0.5*(nu_vec+1));
                    ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma_inv,0.5*(nu_vec+1));
                    
                    loglik_r2(j) = -ln_PWD_Sigma  + ln_PWD_j_yy_Sigma;
                    
                    
                end
                
                loglik_r2 = loglik_r2 - (k/2)*log(pi) + LogLik.lnmultigamma_nu_vec(k, (nu_vec+1)/2 )   -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
                LLF       = -sum(loglik_r2);
                
                if isnan(LLF)
                    LLF = 1e12;
                end
                
            end
        end
        
        
        
        function [LLF,loglik_RC] = LogLik_Inv_Riesz_II(k,T,params,RC);
            
            %nr_l = k*(k+1)/2;
            %lower_vech_Sigma = params(1:nr_l);
            nu_vec = params(4:end);
            
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_II(nu_vec);
            
            U_00_trans = [params(1) params(2); 0 params(3)];
            Sigma_0 = U_00_trans * U_00_trans';
            [chol_Sigma_trans,ind_pd] = chol(Sigma_0);
            
            
            if ind_pd>0
                LLF = 1e14;
            else
                
                U_0_inv_Sigma = chol_Sigma_trans\eye(k);
                
                U_inv_Sigma = U_0_inv_Sigma *diag(sqrt(mu_vec));
                inv_Sigma = U_inv_Sigma * U_inv_Sigma';
                [~,ind_chol] = chol(inv_Sigma);
                
                loglik_RC= zeros(T,1);
                
                if ind_chol>0
                    LLF = 1e14;
                else
                    
                    Sigma = inv(inv_Sigma);
                    ln_PWD_Sigma_U = LogLik.ln_PWD_U(inv_Sigma,-0.5*nu_vec);
                    
                    if isnan(ln_PWD_Sigma_U)
                        LLF = NaN;
                        
                    else
                        
                        for j =1 :T
                            j_RC       = RC(:,:,j);
                            j_inv_RC   = inv(j_RC);
                            ln_PWD_j_RC_U  = LogLik.ln_PWD_U(j_inv_RC,0.5*(nu_vec+k+1));
                            
                            loglik_RC(j) = ln_PWD_Sigma_U  + ln_PWD_j_RC_U - 0.5 * trace(Sigma * j_inv_RC);
                        end
                        
                        loglik_RC = loglik_RC - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_U_nu_vec(k, nu_vec/2 );
                        LLF       = -sum(loglik_RC);
                        
                        if isnan(LLF)
                            LLF = 1e12;
                        end
                    end
                end
            end
        end
        
        
        function [LLF,loglik_RC] = LogLik_Inv_Riesz_II_CT(k,T,params,RC);
            
            % CT_ means with covariance targeting
            
            
            nu_vec = params;
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_II(nu_vec);
            
            Sigma_0 = mean(RC,3);
            [chol_Sigma_trans,ind_pd] = chol(Sigma_0);
            
            if ind_pd>0
                LLF = 1e14;
            else
                
                U_0_inv_Sigma = chol_Sigma_trans\eye(k);
                
                U_inv_Sigma = U_0_inv_Sigma *diag(sqrt(mu_vec));
                inv_Sigma = U_inv_Sigma * U_inv_Sigma';
                [~,ind_chol] = chol(inv_Sigma);
                
                loglik_RC= zeros(T,1);
                
                if ind_chol>0
                    LLF = 1e14;
                else
                    
                    Sigma = inv(inv_Sigma);
                    ln_PWD_Sigma_U = LogLik.ln_PWD_U(inv_Sigma,-0.5*nu_vec);
                    
                    if isnan(ln_PWD_Sigma_U)
                        LLF = NaN;
                        
                    else
                        
                        for j =1 :T
                            j_RC       = RC(:,:,j);
                            j_inv_RC   = inv(j_RC);
                            ln_PWD_j_RC_U  = LogLik.ln_PWD_U(j_inv_RC,0.5*(nu_vec+k+1));
                            
                            loglik_RC(j) = ln_PWD_Sigma_U  + ln_PWD_j_RC_U - 0.5 * trace(Sigma * j_inv_RC);
                        end
                        
                        loglik_RC = loglik_RC - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_U_nu_vec(k, nu_vec/2 );
                        LLF       = -sum(loglik_RC);
                        
                        if isnan(LLF)
                            LLF = 1e12;
                        end
                    end
                end
            end
        end
        
        
        function [LLF] = LogLik_FRiesz_RES(k,T, params,RC)
            
            %nr_l = k*(k+1)/2;
            %nu_vec = params(3:end);
            nu1   = params(4);
            nu2_vec = params(5:end);
            
            %corr_factor = params(4+k:end);
            
            mu_vec_0 = Admin.Compute_1st_moment_HIW(nu2_vec);
            mu_vec = nu1 *mu_vec_0;
            
            L_00_trans = [params(1) 0; params(2) params(3)];
            j_Sigma_0 = L_00_trans * L_00_trans';
            
            Sigma_inv_j_0    = j_Sigma_0\eye(k);
            
            
            [L_0_trans,ind_pd] = chol(Sigma_inv_j_0);
            if ind_pd>0
                LLF = 1e14;
                
            else
                
                L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                inv_Sigma_j = L_trans' * L_trans;
                j_Sigma = inv(inv_Sigma_j);
                
                
                [~,ind_chol] = chol(inv_Sigma_j);
                if ind_chol>0
                    LLF = 1e14;
                else
                    % get loglike
                    loglik_RC= zeros(T,1);
                    
                    for j=1:T
                        j_RC         = RC(:,:,j);
                        %j_inv_RC   = inv(j_RC);
                        
                        %ln_PWD_j_RC  = LogLik.ln_PWD(j_inv_RC,0.5*(nu_vec+k+1));
                        
                        j_X = j_RC + j_Sigma;
                        j_inv_X = j_X\eye(k);
                        
                        [~,ind_chol] = chol(j_inv_X);
                        if ind_chol>0
                            loglik_RC = nan;
                            break
                        end
                        
                        
                        
                        ln_PWD_X_j = LogLik.ln_PWD(j_inv_X,-0.5*(nu2_vec+nu1));
                        ln_PWD_Sigma_j = LogLik.ln_PWD(inv_Sigma_j,-0.5*nu2_vec);
                        j_ln_det_RC = log(det(j_RC));
                        
                        %j_RC        = RC(:,:,j);
                        %j_inv_RC    = inv(j_RC);
                        
                        loglik_RC(j) = ln_PWD_Sigma_j  + 0.5*(nu1-k-1)*j_ln_det_RC - ln_PWD_X_j;
                        
                    end
                    
                    loglik_RC = loglik_RC + LogLik.lnmultigamma_nu_vec(k, (nu2_vec+nu1)/2 ) ...
                        -  LogLik.lnmultigamma_nu_vec(k, nu2_vec/2 ) - LogLik.lnmultigamma(k, nu1/2 );
                    LLF  = -sum(loglik_RC);
                    
                    
                    
                    if  isnan(LLF) || isreal(LLF)==0
                        LLF=1e6;
                    end
                end
            end
        end
        
        
        function [LLF] = LogLik_FRiesz_I_v_gamma(k,T, params,RC)
            % This is the pdf from theorem 11 (Diaz Garcia version)
            
            nu1_vec   = params(4:3+k);
            nu2_vec = params(4+k:end);
            
            mu_vec = Admin.Compute_1st_moment_FRiesz_I(nu1_vec,nu2_vec);
            gamma_t = 0.5*(k+1) - [k:-1:1]';
            
            L_00_trans = [params(1) 0; params(2) params(3)];
            L_trans  = L_00_trans * diag(1./sqrt(mu_vec));
            j_Sigma  = L_trans * L_trans';
            
            
            %Sigma_inv_j_0    = j_Sigma_0\eye(k);
            
            [L_0_trans,ind_pd] = chol(j_Sigma);
            if ind_pd>0
                LLF = 1e14;
            else
                
                ln_PWD_Sigma_j = LogLik.ln_PWD(j_Sigma,-0.5*nu2_vec + gamma_t);
                loglik_RC= zeros(T,1);
                
                for j=1:T
                    j_RC         = RC(:,:,j);
                    ln_PWD_j_RC  = LogLik.ln_PWD(j_RC,0.5*(nu1_vec-k-1));
                    
                    j_X = j_RC + j_Sigma;
                    
                    [~,ind_chol] = chol(j_X);
                    if ind_chol>0
                        loglik_RC = nan;
                        break
                    end
                    
                    ln_PWD_X_j = LogLik.ln_PWD(j_X,-0.5*(nu2_vec+nu1_vec) + gamma_t);
                    loglik_RC(j) = -ln_PWD_Sigma_j + ln_PWD_j_RC + ln_PWD_X_j;
                    
                end
                
                loglik_RC = loglik_RC + LogLik.lnmultigamma_nu_vec(k, (nu2_vec+nu1_vec)/2 ) ...
                    -  LogLik.lnmultigamma_nu_vec(k, nu2_vec/2 ) - LogLik.lnmultigamma_nu_vec(k, nu1_vec/2 );
                LLF  = -sum(loglik_RC);
                
            end
            
            if  isnan(LLF) || isreal(LLF)==0
                LLF=1e14;
            end
            %end
        end
        
        function [LLF] = LogLik_FRiesz_II_v_gamma(k,T, params,RC)
            % This is the pdf from theorem 13 (Diaz Garcia version)
            
            %nr_l = k*(k+1)/2;
            %nu_vec = params(3:end);
            nu1_vec   = params(4:3+k);
            nu2_vec   = params(4+k:end);
            
            mu_vec = Admin.Compute_1st_moment_FRiesz_II(nu1_vec,nu2_vec);
            gamma_t = 0.5*(k+1) - [k:-1:1]';
            
            U_tilde = [params(1) params(2); 0 params(3)];
            U = U_tilde * diag(1./sqrt(mu_vec));
            Sigma = U * U';
            
            [~,ind_pd] = chol(Sigma);
            
            if ind_pd>0
                LLF = 1e14;
            else
                
                ln_PWD_Sigma = LogLik.ln_PWD_U(Sigma,-0.5*nu2_vec - gamma_t);
                loglik_RC= zeros(T,1);
                
                for j=1:T
                    j_RC         = RC(:,:,j);
                    ln_PWD_j_RC  = LogLik.ln_PWD_U(j_RC,0.5*(nu1_vec-k-1));
                    
                    j_X = j_RC + Sigma;
                    
                    [~,ind_chol] = chol(j_X);
                    if ind_chol>0
                        loglik_RC = nan;
                        break
                    end
                    
                    ln_PWD_X_j = LogLik.ln_PWD_U(j_X,-0.5*(nu2_vec+nu1_vec) - gamma_t);
                    loglik_RC(j) = -ln_PWD_Sigma + ln_PWD_j_RC + ln_PWD_X_j;
                    
                end
                
                loglik_RC = loglik_RC + LogLik.lnmultigamma_U_nu_vec(k, (nu2_vec+nu1_vec)/2 ) ...
                    -  LogLik.lnmultigamma_U_nu_vec(k, nu2_vec/2 ) - LogLik.lnmultigamma_U_nu_vec(k, nu1_vec/2 );
                LLF  = -sum(loglik_RC);
                
            end
            
            if  isnan(LLF) || isreal(LLF)==0
                LLF=1e14;
            end
            %end
        end
        
        
        function [LLF] = LogLik_FRiesz_I_v_gamma_CT(k,T, params,RC)
            % This is the pdf from theorem 11 (Diaz Garcia version)
            % with covariance targeting
            
            nu1_vec   = params(1:k);
            nu2_vec = params(1+k:end);
            
            mu_vec = Admin.Compute_1st_moment_FRiesz_I(nu1_vec,nu2_vec);
            gamma_t = 0.5*(k+1) - [k:-1:1]';
            
            Sigma_0 = mean(RC,3);
            [L_0_trans,ind_pd] = chol(Sigma_0);
            L_00_trans = L_0_trans';
            L_trans  = L_00_trans * diag(1./sqrt(mu_vec));
            j_Sigma  = L_trans * L_trans';
            
            
            %Sigma_inv_j_0    = j_Sigma_0\eye(k);
            
            [L_0_trans,ind_pd] = chol(j_Sigma);
            if ind_pd>0
                LLF = 1e14;
            else
                
                %L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                %inv_Sigma_j = L_trans' * L_trans;
                %j_Sigma = inv(inv_Sigma_j);
                
                %[~,ind_chol] = chol(inv_Sigma_j);
                %if ind_chol>0
                %    LLF = 1e14;
                %else
                
                ln_PWD_Sigma_j = LogLik.ln_PWD(j_Sigma,-0.5*nu2_vec + gamma_t);
                loglik_RC= zeros(T,1);
                
                for j=1:T
                    j_RC         = RC(:,:,j);
                    ln_PWD_j_RC  = LogLik.ln_PWD(j_RC,0.5*(nu1_vec-k-1));
                    
                    j_X = j_RC + j_Sigma;
                    
                    [~,ind_chol] = chol(j_X);
                    if ind_chol>0
                        loglik_RC = nan;
                        break
                    end
                    
                    ln_PWD_X_j = LogLik.ln_PWD(j_X,-0.5*(nu2_vec+nu1_vec) + gamma_t);
                    loglik_RC(j) = -ln_PWD_Sigma_j + ln_PWD_j_RC + ln_PWD_X_j;
                    
                end
                
                loglik_RC = loglik_RC + LogLik.lnmultigamma_nu_vec(k, (nu2_vec+nu1_vec)/2 ) ...
                    -  LogLik.lnmultigamma_nu_vec(k, nu2_vec/2 ) - LogLik.lnmultigamma_nu_vec(k, nu1_vec/2 );
                LLF  = -sum(loglik_RC);
                
            end
            
            if  isnan(LLF) || isreal(LLF)==0
                LLF=1e14;
            end
            %end
        end
        
        
        function [LLF] = LogLik_FRiesz_II_v_gamma_CT(k,T, params,RC)
            
            % This is the pdf from theorem 13 (Diaz Garcia version)
            % with covariance targeting
            
            nu1_vec   = params(1:k);
            nu2_vec   = params(1+k:end);
            
            mu_vec = Admin.Compute_1st_moment_FRiesz_II(nu1_vec,nu2_vec);
            gamma_t = 0.5*(k+1) - [k:-1:1]';
            
            Sigma_0 = mean(RC,3);
            U_tilde = chol(Sigma_0\eye(k))\eye(k);
            
            U = U_tilde * diag(1./sqrt(mu_vec));
            Sigma = U * U';
            
            [~,ind_pd] = chol(Sigma);
            
            if ind_pd>0
                LLF = 1e14;
            else
                
                ln_PWD_Sigma = LogLik.ln_PWD_U(Sigma,-0.5*nu2_vec - gamma_t);
                loglik_RC= zeros(T,1);
                
                for j=1:T
                    j_RC         = RC(:,:,j);
                    ln_PWD_j_RC  = LogLik.ln_PWD_U(j_RC,0.5*(nu1_vec-k-1));
                    
                    j_X = j_RC + Sigma;
                    
                    [~,ind_chol] = chol(j_X);
                    if ind_chol>0
                        loglik_RC = nan;
                        break
                    end
                    
                    ln_PWD_X_j = LogLik.ln_PWD_U(j_X,-0.5*(nu2_vec+nu1_vec) - gamma_t);
                    loglik_RC(j) = -ln_PWD_Sigma + ln_PWD_j_RC + ln_PWD_X_j;
                    
                end
                
                loglik_RC = loglik_RC + LogLik.lnmultigamma_U_nu_vec(k, (nu2_vec+nu1_vec)/2 ) ...
                    -  LogLik.lnmultigamma_U_nu_vec(k, nu2_vec/2 ) - LogLik.lnmultigamma_U_nu_vec(k, nu1_vec/2 );
                LLF  = -sum(loglik_RC);
                
            end
            
            if  isnan(LLF) || isreal(LLF)==0
                LLF=1e14;
            end
            %end
        end
        
        
        function [LLF] = LogLik_FRiesz_NO_SIGMA(k,T, params,RC)
            % old stuff
            
            
            %nr_l = k*(k+1)/2;
            %nu_vec = params(3:end);
            nu1_vec   = params(1:k);
            nu2_vec = params(1+k:end);
            
            loglik_RC= zeros(T,1);
            
            for j=1:T
                j_RC         = RC(:,:,j);
                %j_inv_RC   = inv(j_RC);
                
                ln_PWD_j_RC  = LogLik.ln_PWD(j_RC,0.5*(nu1_vec-k-1));
                
                j_X = j_RC + eye(k);
                ln_PWD_X_j = LogLik.ln_PWD(j_X,-0.5*(nu2_vec+nu1_vec));
                loglik_RC(j) = ln_PWD_j_RC + ln_PWD_X_j;
                
                
            end
            
            loglik_RC = loglik_RC + LogLik.lnmultigamma_nu_vec(k, (nu2_vec+nu1_vec)/2 ) ...
                -  LogLik.lnmultigamma_nu_vec(k, nu2_vec/2 ) - LogLik.lnmultigamma_nu_vec(k, nu1_vec/2 );
            LLF  = -sum(loglik_RC);
            
            
            
            if  isnan(LLF) || ~isreal(LLF)
                LLF=1e6;
            end
            %end
        end
        
        
        function [LLF,loglik_RC] = LogLik_CAIW(k,T,params,RC);
            
            %omega = params(1);
            %nr_l = k*(k+1)/2;
            %lower_vech_Omega = params(1:nr_l);
            %alfa_beta = [0.5246 0.4754]';
            
            nu = params(end);
            %Omega = Admin.Vech2Sym(k, lower_vech_Omega);
            
            [V,Flag] = Filter_New.CAW(k,1,T,params(1:2),RC);
            constant = (nu-k-1);
            
            if Flag==1
                LLF = 1e14;
            else
                
                loglik_RC= zeros(1,T);
                
                for j =1 :T
                    j_Sigma     = V(:,:,j);
                    j_detSigma = det(j_Sigma);
                    
                    j_RC        = RC(:,:,j);
                    j_detRC     = det(j_RC);
                    j_inv_RC    = inv(j_RC);
                    
                    loglik_RC(j) = (nu/2) * log(j_detSigma) -0.5*(nu+k+1)*log(j_detRC) -constant/2 * trace(j_inv_RC*j_Sigma);
                    
                end
                
                loglik_RC = loglik_RC + 0.5*nu*k*log(constant/2) - LogLik.lnmultigamma(k, 0.5*nu);
                LLF = -sum(loglik_RC);
                
                if isnan(LLF)
                    LLF = 1e12;
                    
                end
            end
        end
        
        
        
        function [LLF,V,loglik_mat, mv_beta] = LogLik_CAW_Inv_Riesz_I(k,p,T, params,RC)
            
            %nr_l = k*(k+1)/2;
            %nu_vec = params(3:end);
            nu_vec = params(3:3+k-1);
            %corr_factor = params(3+k:end);
            
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            %[V,Flag] = Filter_New.CAW_diag(k,p,T, params,RC);
            [V,Flag] = Filter_New.CAW(k,p,T,params(1:2),RC);
            
            if Flag==1
                LLF = 1e6;
                
            else
                % get loglike
                loglik_RC= zeros(T,1);
                loglik_mat = zeros(T,3);
                
                for j=1:T
                    j_RC         = RC(:,:,j);
                    j_inv_RC   = inv(j_RC);
                    ln_PWD_j_RC  = LogLik.ln_PWD(j_inv_RC,0.5*(nu_vec+k+1));
                    
                    j_Sigma_0        = V(:,:,j);
                    Sigma_inv_j_0    = j_Sigma_0\eye(k);
                    
                    L_0_trans = chol(Sigma_inv_j_0);
                    L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                    inv_Sigma_j = L_trans' * L_trans;
                    j_Sigma = inv(inv_Sigma_j);
                    
                    [~,ind_chol] = chol(inv_Sigma_j);
                    if ind_chol>0
                        Flag=1;
                        break;
                    end
                    
                    ln_PWD_Sigma_j = LogLik.ln_PWD(inv_Sigma_j,-0.5*nu_vec);
                    loglik_RC(j) = ln_PWD_Sigma_j  + ln_PWD_j_RC - 0.5 * trace(j_Sigma  * j_inv_RC);
                    
                    loglik_mat(j,:) = [ln_PWD_Sigma_j, ln_PWD_j_RC, -0.5*trace(j_Sigma  * j_inv_RC)];
                end
                
                loglik_RC = loglik_RC - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
                LLF       = -sum(loglik_RC);
                
                mv_beta= - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
                
                if  isnan(LLF) || ~isreal(LLF) || Flag==1
                    LLF=1e6;
                end
            end
        end
        
        
        function [LLF,loglik_RC] = LogLik_CAW_Inv_Riesz_II(k,p,T,params,RC);
            
            nu_vec = params(3:end);
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_II(nu_vec);
            
            [V,Flag] = Filter_New.CAW(k,p,T,params(1:2),RC);
            
            if Flag==1
                LLF = 1e6;
                
            else
                
                loglik_RC= zeros(T,1);
                
                for j =1 :T
                    j_RC       = RC(:,:,j);
                    j_inv_RC   = inv(j_RC);
                    ln_PWD_j_RC_U  = LogLik.ln_PWD_U(j_inv_RC,0.5*(nu_vec+k+1));
                    
                    j_Sigma_0        = V(:,:,j);
                    U_0_inv_Sigma_j  = chol(j_Sigma_0)\eye(k);
                    U_inv_Sigma_j    = U_0_inv_Sigma_j *diag(sqrt(mu_vec));
                    inv_Sigma_j        = U_inv_Sigma_j * U_inv_Sigma_j';
                    
                    [~,ind_chol] = chol(inv_Sigma_j);
                    
                    if ind_chol>0
                        Flag=1;
                        break;
                    end
                    
                    Sigma_j = inv(inv_Sigma_j);
                    ln_PWD_Sigma_U = LogLik.ln_PWD_U(inv_Sigma_j,-0.5*nu_vec);
                    
                    loglik_RC(j) = ln_PWD_Sigma_U  + ln_PWD_j_RC_U - 0.5 * trace(Sigma_j * j_inv_RC);
                end
                
                loglik_RC = loglik_RC - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_U_nu_vec(k, nu_vec/2 );
                LLF       = -sum(loglik_RC);
                
                if isnan(LLF) || ~isreal(LLF) || Flag ==1
                    LLF = 1e12;
                end
            end
        end
        
        
        function [LLF,V] = LogLik_CAW_FRiesz_RES(k,p,T, params,RC)
            
            %nr_l = k*(k+1)/2;
            %nu_vec = params(3:end);
            nu1   = params(3);
            nu2_vec = params(4:4+k-1);
            
            %corr_factor = params(4+k:end);
            
            
            mu_vec_0 = Admin.Compute_1st_moment_Inv_Riesz_I(nu2_vec);
            mu_vec = nu1 *mu_vec_0;
            
            %[V,Flag] = Filter_New.CAW_diag(k,p,T, params,RC);
            [V,Flag] = Filter_New.CAW(k,p,T,params(1:2),RC);
            
            if Flag==1
                LLF = 1e6;
                
            else
                % get loglike
                loglik_RC= zeros(T,1);
                
                for j=1:T
                    j_RC         = RC(:,:,j);
                    %j_inv_RC   = inv(j_RC);
                    
                    %ln_PWD_j_RC  = LogLik.ln_PWD(j_inv_RC,0.5*(nu_vec+k+1));
                    
                    j_Sigma_0        = V(:,:,j);
                    Sigma_inv_j_0    = j_Sigma_0\eye(k);
                    
                    L_0_trans = chol(Sigma_inv_j_0);
                    %L_trans = diag(sqrt(corr_factor)) * L_0_trans;
                    L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                    inv_Sigma_j = L_trans' * L_trans;
                    j_Sigma = inv(inv_Sigma_j);
                    
                    
                    [~,ind_chol] = chol(inv_Sigma_j);
                    if ind_chol>0
                        loglik_RC = nan;
                        break
                    end
                    
                    j_X = j_RC + j_Sigma;
                    j_inv_X = j_X\eye(k);
                    
                    [~,ind_chol] = chol(j_inv_X);
                    if ind_chol>0
                        loglik_RC = nan;
                        break
                    end
                    
                    
                    
                    ln_PWD_X_j = LogLik.ln_PWD(j_inv_X,-0.5*(nu2_vec+nu1));
                    ln_PWD_Sigma_j = LogLik.ln_PWD(inv_Sigma_j,-0.5*nu2_vec);
                    j_ln_det_RC = log(det(j_RC));
                    
                    %j_RC        = RC(:,:,j);
                    %j_inv_RC    = inv(j_RC);
                    
                    loglik_RC(j) = ln_PWD_Sigma_j  + 0.5*(nu1-k-1)*j_ln_det_RC - ln_PWD_X_j;
                    
                end
                
                loglik_RC = loglik_RC + LogLik.lnmultigamma_nu_vec(k, (nu2_vec+nu1)/2 ) ...
                    -  LogLik.lnmultigamma_nu_vec(k, nu2_vec/2 ) - LogLik.lnmultigamma(k, nu1/2 );
                LLF  = -sum(loglik_RC);
                
                
                
                if  isnan(LLF) || isreal(LLF)==0
                    LLF=1e6;
                end
            end
        end
        
        
        function [LLF,V,loglik_mat,mv_beta] = LogLik_CAW_FRiesz_I_v_gamma(k,p,T, params,RC)
            
            %nr_l = k*(k+1)/2;
            %nu_vec = params(3:end);
            
            nu1_vec   = params(3:3+k-1);
            nu2_vec = params(3+k:3+2*k-1);
            %3+2*k-1);
            
            gamma_t = 0.5*(k+1) - [k:-1:1]';
            
            %corr_factor = params(3+2*k:end);
            
            mu_vec = Admin.Compute_1st_moment_FRiesz_I(nu1_vec,nu2_vec);
            
            %mu_vec = params(3+2*k:end);
            
            %[V,Flag] = Filter_New.CAW_diag(k,p,T, params,RC);
            [V,Flag] = Filter_New.CAW(k,p,T,params(1:2),RC);
            
            if Flag==1
                LLF = 1e6;
                
            else
                % get loglike
                loglik_RC= zeros(T,1);
                loglik_mat= zeros(T,3);
                
                for j=1:T
                    j_RC         = RC(:,:,j);
                    
                    ln_PWD_j_RC  = LogLik.ln_PWD(j_RC,0.5*(nu1_vec-k-1));
                    
                    j_Sigma_0        = V(:,:,j);
                    L_trans_0        = chol(j_Sigma_0);
                    L_trans          = diag(1./sqrt(mu_vec))*L_trans_0;
                    j_Sigma          = L_trans' * L_trans;
                    
                    [~,ind_chol] = chol(j_Sigma);
                    if ind_chol>0
                        loglik_RC = nan;
                        break
                    end
                    
                    j_X = j_RC + j_Sigma;
                    
                    ln_PWD_X_j = LogLik.ln_PWD(j_X,-0.5*(nu2_vec+nu1_vec)+gamma_t);
                    ln_PWD_Sigma_j = LogLik.ln_PWD(j_Sigma,-0.5*nu2_vec+gamma_t);
                    
                    loglik_RC(j) = -ln_PWD_Sigma_j  + ln_PWD_j_RC + ln_PWD_X_j;
                    
                    loglik_mat(j,:) = [-ln_PWD_Sigma_j, ln_PWD_j_RC, ln_PWD_X_j];
                end
                
                loglik_RC = loglik_RC + LogLik.lnmultigamma_nu_vec(k, (nu2_vec+nu1_vec)/2 ) ...
                    -  LogLik.lnmultigamma_nu_vec(k, nu2_vec/2 ) - LogLik.lnmultigamma_nu_vec(k, nu1_vec/2 );
                LLF  = -sum(loglik_RC);
                
                mv_beta = LogLik.lnmultigamma_nu_vec(k, (nu2_vec+nu1_vec)/2 ) ...
                    -  LogLik.lnmultigamma_nu_vec(k, nu2_vec/2 ) - LogLik.lnmultigamma_nu_vec(k, nu1_vec/2 );
                
                if  isnan(LLF) || ~isreal(LLF)
                    LLF=1e6;
                end
            end
        end
        
        
        function [LLF] = LogLik_CAW_FRiesz_II_v_gamma(k,p,T, params,RC)
            
            nu1_vec   = params(3:k+2);
            nu2_vec   = params(3+k:end);
            
            mu_vec = Admin.Compute_1st_moment_FRiesz_II(nu1_vec,nu2_vec);
            gamma_t = 0.5*(k+1) - [k:-1:1]';
            
            
            [V,Flag] = Filter_New.CAW(k,p,T,params(1:2),RC);
            
            if Flag==1
                LLF = 1e14;
            else
                
                loglik_RC= zeros(T,1);
                
                for j=1:T
                    j_RC         = RC(:,:,j);
                    ln_PWD_j_RC  = LogLik.ln_PWD_U(j_RC,0.5*(nu1_vec-k-1));
                    
                    j_Sigma_0        = V(:,:,j);
                    U_tilde          = chol(j_Sigma_0\eye(k))\eye(k);
                    U                = U_tilde * diag(1./sqrt(mu_vec));
                    Sigma_j          = U * U';
                    
                    [~,ind_chol] = chol(Sigma_j);
                    if ind_chol>0
                        loglik_RC = nan;
                        break
                    end
                    
                    ln_PWD_Sigma_j = LogLik.ln_PWD_U(Sigma_j,-0.5*nu2_vec - gamma_t);
                    
                    j_X = j_RC + Sigma_j;
                    
                    %[~,ind_chol] = chol(j_X);
                    %if ind_chol>0
                    %    loglik_RC = nan;
                    %    break
                    %end
                    
                    ln_PWD_X_j = LogLik.ln_PWD_U(j_X,-0.5*(nu2_vec+nu1_vec) - gamma_t);
                    loglik_RC(j) = -ln_PWD_Sigma_j + ln_PWD_j_RC + ln_PWD_X_j;
                    
                end
                
                loglik_RC = loglik_RC + LogLik.lnmultigamma_U_nu_vec(k, (nu2_vec+nu1_vec)/2 ) ...
                    -  LogLik.lnmultigamma_U_nu_vec(k, nu2_vec/2 ) - LogLik.lnmultigamma_U_nu_vec(k, nu1_vec/2 );
                LLF  = -sum(loglik_RC);
                
            end
            
            if  isnan(LLF) || ~isreal(LLF)
                LLF=1e14;
            end
            %end
        end
        
        
        function [LLF,loglik_RC] = LogLik_CAW_Riesz_II(k,p,T,params,RC);
            
            nu_vec = params(3:end);
            [V,Flag] = Filter_New.CAW(k,p,T,params(1:2),RC);
            
            
            if Flag==1
                LLF = 1e14;
            else
                
                loglik_RC= zeros(T,1);
                
                
                for j=1:T
                    
                    j_RC       = RC(:,:,j);
                    ln_PWD_j_RC_U  = LogLik.ln_PWD_U(j_RC,0.5*(nu_vec-k-1));
                    
                    j_V          = V(:,:,j);
                    U_tilde      = chol(j_V\eye(k))\eye(k);
                    U = U_tilde * diag(1./sqrt(nu_vec));
                    Sigma_j = U * U';
                    
                    [~,ind_chol] = chol(Sigma_j);
                    if ind_chol>0
                        loglik_RC = nan;
                        break
                    end
                    
                    ln_PWD_Sigma_U = LogLik.ln_PWD_U(Sigma_j,0.5*nu_vec);
                    inv_Sigma_j = inv(Sigma_j);
                    loglik_RC(j) = -ln_PWD_Sigma_U  + ln_PWD_j_RC_U - 0.5 * trace(inv_Sigma_j * j_RC);
                    
                    
                end
                
                loglik_RC = loglik_RC - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_U_nu_vec(k, nu_vec/2 );
                LLF       = -sum(loglik_RC);
                
                if isnan(LLF)
                    LLF = 1e12;
                end
                
            end
        end
        
        
        function [LLF,V,loglik_mat,MV_beta] = LogLik_CAW_Riesz_I(k,p,T, params,RC)
            
            nu_vec = params(3:end);
            [V,Flag] = Filter_New.CAW(k,p,T,params(1:2),RC);
            
            if Flag==1
                LLF = 1e14;
                
            else
                
                % get loglike
                loglik_RC= zeros(T,1);
                loglik_mat = zeros(T,3);
                
                for j=1:T
                    j_RC         = RC(:,:,j);
                    ln_PWD_j_RC  = LogLik.ln_PWD(j_RC,0.5*(nu_vec-k-1));
                    
                    j_V          = V(:,:,j);
                    L_tilde      = chol(j_V)';
                    L            = L_tilde * diag(1./sqrt(nu_vec));
                    Sigma_j      = L * L';
                    
                    [~,ind_chol] = chol(Sigma_j);
                    if ind_chol>0
                        loglik_RC = nan;
                        break
                    end
                    
                    ln_PWD_j_Sigma = LogLik.ln_PWD(Sigma_j,0.5*nu_vec);
                    j_inv_Sigma       = Sigma_j\eye(k);
                    
                    loglik_RC(j,1) = -ln_PWD_j_Sigma  + ln_PWD_j_RC - 0.5 * trace(j_inv_Sigma * j_RC);
                    
                    loglik_mat(j,:) = [-ln_PWD_j_Sigma, ln_PWD_j_RC, -0.5*trace(j_inv_Sigma * j_RC)];
                    
                end
                
                loglik_RC = loglik_RC - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
                
                MV_beta = - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
                
                LLF= -sum(loglik_RC);
                
                if  isnan(LLF) || ~isreal(LLF)
                    LLF=1e6;
                end
            end
        end
        
        
        
        
        function [LLF,V] = LogLik_EWMA_Matrix_F_RC(k,T,params,RC);
            
            beta  = params(1);
            nu_f1 = params(2);
            nu_f2 = params(3);
            
            constant = nu_f1/ (nu_f2-k-1);
            [V] = Filter_New.EWMA(k, T, beta,RC);
            %[V,~,Flag]  = Filter_New.GAS_Gen_CT( k, T, [alpha_vec beta nu_t nu_f1 nu_f2], r2, RC );
            
            % get loglike
            loglik_RC= zeros(1,T);
            
            for j=1:T
                j_invV       = inv(V(:,:,j));
                j_RC         = RC(:,:,j);
                j_logdetRC   = log(det(RC(:,:,j)));
                
                %loglik_r2(j) = - 1/2 * j_logdetV  - 1/2 * trace( j_invV *r2(:,:,j) ); normal
                
                % matrix F
                loglik_RC(j) = 0.5*nu_f1 * log(det(constant*j_invV)) + 0.5*(nu_f1-k-1)*j_logdetRC - 0.5*(nu_f1+nu_f2)*log(det(eye(k)+ constant*j_invV*j_RC));
            end
            
            %loglik_r2= loglik_r2 - k/2 * log(2*pi); % normal
            
            loglik_RC = loglik_RC + LogLik.lnmultigamma(k, 0.5*(nu_f1+nu_f2)) - LogLik.lnmultigamma(k, 0.5*nu_f1) - LogLik.lnmultigamma(k, 0.5*nu_f2);
            
            LLF = -sum(loglik_RC);
            
            if  isnan(LLF)
                LLF=1e6;
            end
            
        end
        
        
        function [LLF,loglik_RC] = LogLik_Mult_T(k,T,params,r2,RC);
            
            nu_t = params(1);
            
            %vechDelta = params(2:end);
            %Delta = Admin.Vech2Sym(k,vechDelta);
            
            %[test1,test2] = chol(Delta);
            %if test2>0
            %    LLF = 1e12;
            %else
            
            
            loglik_RC= zeros(T,1);
            
            for j =1 :T
                inv_Delta = inv(RC(:,:,j));
                logdet_Delta = log(det(RC(:,:,j)));
                loglik_r2(j) = - 1/2 * logdet_Delta  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(inv_Delta*r2(:,:,j)));
                %loglik_r2(j) = - 1/2 * logdet_Delta  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*r(j,:)*inv_Delta*r(j,:)');
                
            end
            
            loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);
            LLF = -sum(loglik_r2);
            
            if isnan(LLF)
                LLF = 1e12;
                
            end
            % end
            
            
            
        end
        
        
        function [logL,Rt_mat,Qbar_cDCC] = LogLik_cDCC_PSEUDO_N(k,T,T_w,params, Stdresid)
            
            % not the true loglikelihood
            % only part that corresponds with the correlation part
            Flag = 0;
            
            Qbar_DCC = cov(Stdresid);
            Rt_mat = zeros(k,k,T);
            Qt_current = Qbar_DCC;
            alpha = params(1);
            beta =  params(2);
            
            logL_vec = zeros(T,1);
            
            q_ii_mat = zeros(T,k);
            q_ii_mat(1,:) = diag(Qbar_DCC);
            Qbar_new  = zeros(k,k);
            for m = 1:T
                if m<T
                    q_ii_mat(m+1,:) = (1-sum(params)) + alpha*q_ii_mat(m,:).*Stdresid(m,:).^2 + beta*q_ii_mat(m,:);
                end
                Q_square_root_star_t = diag(sqrt(q_ii_mat(m,:)));
                Qbar_new = Qbar_new +  Q_square_root_star_t * Stdresid(m,:)'* Stdresid(m,:) * Q_square_root_star_t;
            end
            
            Qbar_cDCC = Qbar_new/T;
            
            for j = 1:T
                
                Rt_DCC= Qt_current./(sqrt(diag(Qt_current))*sqrt(diag(Qt_current))');
                
                if j<T
                    temp = sqrt(diag(Qt_current));
                    Q_star_t = diag(temp(1:end));
                    Qt_new = Qbar_cDCC*(1 - alpha-beta) + alpha*(Q_star_t*Stdresid(j,:)'*Stdresid(j,:)*Q_star_t) + beta*Qt_current;
                end
                
                if rcond(Rt_DCC)<1e-12                   
                    %LogL_vec(j,1) = NaN;
                    Flag = 1;
                    break
                end
                
                Rt_mat(:,:,j) = Rt_DCC;
                logL_vec(j,1) = -log(det(Rt_DCC)) - Stdresid(j,:)*inv(Rt_DCC)*Stdresid(j,:)' + Stdresid(j,:)*Stdresid(j,:)';
                
                Qt_current = Qt_new;
            end
            
            %Q_T_plus_1 = Qt_current;
            %R_T_plus_1 = Rt_current;
            if Flag==1
                logL = 1e14;
            else
                logL = -(1/2)*sum(logL_vec(63+1:end));
            end
        end
        
        
        function [logL,Rt_mat,Q_bar_t] = LogLik_cDCC_PSEUDO_N_RC(k,T,T_w,params, Stdresid,RCORR_mat)
            
            % not the true loglikelihood
            % only part that corresponds with the correlation part
            
            Qbar_DCC = cov(Stdresid);
            Qbar_DCC_t = mean(RCORR_mat,3);
            Rt_mat = zeros(k,k,T);
            Q_bar_t = zeros(k,k,T);
            Qt_current = Qbar_DCC;
            alpha = params(1);
            beta =  params(2);
            
            logL_vec = zeros(T,1);
            
%             q_ii_mat = zeros(T,k);
%             q_ii_mat(1,:) = diag(Qbar_DCC);
%             Qbar_new  = zeros(k,k);
%             for m = 1:T
%                 if m<T
%                     q_ii_mat(m+1,:) = (1-sum(params)) + alpha*q_ii_mat(m,:).*Stdresid(m,:).^2 + beta*q_ii_mat(m,:);
%                 end
%                 Q_square_root_star_t = diag(sqrt(q_ii_mat(m,:)));
%                 Qbar_new = Qbar_new +  Q_square_root_star_t * Stdresid(m,:)'* Stdresid(m,:) * Q_square_root_star_t;
%             end
%             
%             Qbar_cDCC = Qbar_new/T;
            
            for j = 1:T
                
                Rt_DCC= Qt_current./(sqrt(diag(Qt_current))*sqrt(diag(Qt_current))');
                
                if j<T
                    
                    if j>=T_w
                         Qbar_DCC_t = mean(RCORR_mat(:,:,j-T_w+1:j),3);         
                    end
                    
                    temp = sqrt(diag(Qt_current));
                    Q_star_t = diag(temp(1:end));
                    Qt_new = Qbar_DCC_t*(1 - alpha-beta) + alpha*(Q_star_t*Stdresid(j,:)'*Stdresid(j,:)*Q_star_t) + beta*Qt_current;
                end
                
                Q_bar_t(:,:,j) = Qbar_DCC_t;
                
                Rt_mat(:,:,j) = Rt_DCC;
                logL_vec(j,1) = -log(det(Rt_DCC)) - Stdresid(j,:)*inv(Rt_DCC)*Stdresid(j,:)' + Stdresid(j,:)*Stdresid(j,:)';
               
                Qt_current = Qt_new;
            end
            
            %Q_T_plus_1 = Qt_current;
            %R_T_plus_1 = Rt_current;
            logL = -(1/2)*sum(logL_vec(63+1:end));
        end
        
        
         
        function [logL,Rt_mat,Qbar_DCC] = LogLik_REALIZED_DCC_PSEUDO_N(k,T,T_w,params, Stdresid,RCORR_mat)
            
            % not the true loglikelihood
            % only part that corresponds with the correlation part
            
           % Qbar_DCC = cov(Stdresid);
            Qbar_DCC = mean(RCORR_mat,3);
            Rt_mat = zeros(k,k,T);
            %Q_bar_t = zeros(k,k,T);
            Qt_current = Qbar_DCC;
            alpha = params(1);
            beta =  params(2);
            
            logL_vec = zeros(T,1);
            
%             q_ii_mat = zeros(T,k);
%             q_ii_mat(1,:) = diag(Qbar_DCC);
%             Qbar_new  = zeros(k,k);
%             for m = 1:T
%                 if m<T
%                     q_ii_mat(m+1,:) = (1-sum(params)) + alpha*q_ii_mat(m,:).*Stdresid(m,:).^2 + beta*q_ii_mat(m,:);
%                 end
%                 Q_square_root_star_t = diag(sqrt(q_ii_mat(m,:)));
%                 Qbar_new = Qbar_new +  Q_square_root_star_t * Stdresid(m,:)'* Stdresid(m,:) * Q_square_root_star_t;
%             end
%             
%             Qbar_cDCC = Qbar_new/T;
            
            for j = 1:T
                
                %Rt_DCC= Qt_current./(sqrt(diag(Qt_current))*sqrt(diag(Qt_current))');
                Rt_DCC = Qt_current;
                if j<T                                        
                    
                    %temp = sqrt(diag(Qt_current));
                    %Q_star_t = diag(temp(1:end));
                    Qt_new = Qbar_DCC*(1 - alpha-beta) + alpha*RCORR_mat(:,:,j) + beta*Qt_current;
                end
                
                %Q_bar_t(:,:,j) = Qbar_DCC_t;
                
                Rt_mat(:,:,j) = Rt_DCC;
                %Rt_mat(:,:,j) = Qt_current;
                logL_vec(j,1) = -log(det(Rt_DCC)) - Stdresid(j,:)*inv(Rt_DCC)*Stdresid(j,:)' + Stdresid(j,:)*Stdresid(j,:)';
               
                Qt_current = Qt_new;
            end
            
            %Q_T_plus_1 = Qt_current;
            %R_T_plus_1 = Rt_current;
            logL = -(1/2)*sum(logL_vec(63+1:end));
        end
        
        
         function [logL,Rt_mat,Qbar_t] = LogLik_REALIZED_DCC_PSEUDO_N2(k,T,T_w,params, Stdresid,RCORR_mat)
            
            % not the true loglikelihood
            % only part that corresponds with the correlation part
            
           % Qbar_DCC = cov(Stdresid);
            Qbar_DCC_t = mean(RCORR_mat,3);
            Rt_mat = zeros(k,k,T);
            Qbar_t = zeros(k,k,T);            
            Qt_current = Qbar_DCC_t;
            alpha = params(1);
            beta =  params(2);
            
            logL_vec = zeros(T,1);
            
%             q_ii_mat = zeros(T,k);
%             q_ii_mat(1,:) = diag(Qbar_DCC);
%             Qbar_new  = zeros(k,k);
%             for m = 1:T
%                 if m<T
%                     q_ii_mat(m+1,:) = (1-sum(params)) + alpha*q_ii_mat(m,:).*Stdresid(m,:).^2 + beta*q_ii_mat(m,:);
%                 end
%                 Q_square_root_star_t = diag(sqrt(q_ii_mat(m,:)));
%                 Qbar_new = Qbar_new +  Q_square_root_star_t * Stdresid(m,:)'* Stdresid(m,:) * Q_square_root_star_t;
%             end
%             
%             Qbar_cDCC = Qbar_new/T;
            
            for j = 1:T
                
                %Rt_DCC= Qt_current./(sqrt(diag(Qt_current))*sqrt(diag(Qt_current))');
                Rt_DCC = Qt_current;
                if j<T                                        
                    
                    if j>=T_w
                        Qbar_DCC_t = mean(RCORR_mat(:,:,j-T_w+1:j),3);                        
                    end
                    
                    %temp = sqrt(diag(Qt_current));
                    %Q_star_t = diag(temp(1:end));
                    Qt_new = Qbar_DCC_t*(1 - alpha-beta) + alpha*RCORR_mat(:,:,j) + beta*Qt_current;
                end
                
                Qbar_t(:,:,j) = Qbar_DCC_t;
                
                Rt_mat(:,:,j) = Rt_DCC;
                %Rt_mat(:,:,j) = Qt_current;
                logL_vec(j,1) = -log(det(Rt_DCC)) - Stdresid(j,:)*inv(Rt_DCC)*Stdresid(j,:)' + Stdresid(j,:)*Stdresid(j,:)';
               
                Qt_current = Qt_new;
            end
            
            %Q_T_plus_1 = Qt_current;
            %R_T_plus_1 = Rt_current;
            logL = -(1/2)*sum(logL_vec(63+1:end));
        end
        
        
        
        
        
        function [LLF,Rt_mat,Qbar_cDCC] = LogLik_cDCC_N(k,T,params, Stdresid)
            
            % not the true loglikelihood
            % only part that corresponds with the correlation part
            
            Qbar_DCC = cov(Stdresid);
            D_mat = diag(sqrt(diag(Qbar_DCC)));
            
            %Rt_mat = [];
            Rt_mat = zeros(k,k,T);
            Qt_current = Qbar_DCC;
            alpha = params(1);
            beta =  params(2);
            
            loglik_vec = zeros(T,1);
            
            q_ii_mat = zeros(T,k);
            q_ii_mat(1,:) = diag(Qbar_DCC);
            Qbar_new  = zeros(k,k);
            for m = 1:T
                if m<T
                    q_ii_mat(m+1,:) = (1-sum(params)) + alpha*q_ii_mat(m,:).*Stdresid(m,:).^2 + beta*q_ii_mat(m,:);
                end
                Q_square_root_star_t = diag(sqrt(q_ii_mat(m,:)));
                Qbar_new = Qbar_new +  Q_square_root_star_t * Stdresid(m,:)'* Stdresid(m,:) * Q_square_root_star_t;
            end
            
            Qbar_cDCC = Qbar_new/T;
            
            for j = 1:T
                
                Rt_DCC= Qt_current./(sqrt(diag(Qt_current))*sqrt(diag(Qt_current))');
                
                if j<T
                    temp = sqrt(diag(Qt_current));
                    Q_star_t = diag(temp(1:end));
                    Qt_new = Qbar_cDCC*(1 - alpha-beta) + alpha*(Q_star_t*Stdresid(j,:)'*Stdresid(j,:)*Q_star_t) + beta*Qt_current;
                end
                
                Rt_mat(:,:,j) = Rt_DCC;
                
                Ht_DCC = D_mat*Rt_DCC*D_mat;
                [~,ind_chol] = chol(Ht_DCC);
                if ind_chol>0
                    loglik_vec = NaN;
                    break;
                end
                
                loglik_vec(j,1) = -0.5*log(det(Ht_DCC)) - 0.5 * Stdresid(j,:)*inv(Ht_DCC)*Stdresid(j,:)';
                
                Qt_current = Qt_new;
            end
            
            loglik_vec = loglik_vec  - 0.5*k*log(2*pi);
            LLF = -sum(loglik_vec);
            
            if isnan(LLF)
                LLF = 1e12;
            end
            
        end
        
        
        function [LLF,Rt_mat] = LogLik_cDCC_X_N(k,T,params, Stdresid,X_mat)
            
            % not the true loglikelihood
            % only part that corresponds with the correlation part
            
            Qbar_DCC = cov(Stdresid);
            X_bar = mean(X_mat,3);
            %Qbar_Anne = [1 params(4); params(4) 1];
            D_mat = diag(sqrt(diag(Qbar_DCC)));
            
            %Rt_mat = [];
            Rt_mat = zeros(k,k,T);
            Qt_current = Qbar_DCC;
            alpha = params(1);
            beta =  params(2);
            delta = params(3);
            
            loglik_vec = zeros(T,1);
            
            %             q_ii_mat = zeros(T,k);
            %             q_ii_mat(1,:) = diag(Qbar_DCC);
            %             Qbar_new  = zeros(k,k);
            %             for m = 1:T
            %                 if m<T
            %                     q_ii_mat(m+1,:) = (1-sum(params)) + alpha*q_ii_mat(m,:).*Stdresid(m,:).^2 + beta*q_ii_mat(m,:);
            %                 end
            %                 Q_square_root_star_t = diag(sqrt(q_ii_mat(m,:)));
            %                 Qbar_new = Qbar_new +  Q_square_root_star_t * Stdresid(m,:)'* Stdresid(m,:) * Q_square_root_star_t;
            %             end
            %
            %             Qbar_cDCC = Qbar_new/T;
            
            for j = 1:T
                
                Rt_DCC= Qt_current./(sqrt(diag(Qt_current))*sqrt(diag(Qt_current))');
                
                if j<T
                    %temp = sqrt(diag(Qt_current));
                    %Q_star_t = diag(temp(1:end));
                    
                    %Qt_new = Qbar_DCC*(1 - alpha-beta) + alpha*(Q_star_t*Stdresid(j,:)'*Stdresid(j,:)*Q_star_t) + delta*RCORR+ beta*Qt_current;
                    Qt_new = Qbar_DCC*(1 - alpha-beta) - delta*X_bar + alpha*Stdresid(j,:)'*Stdresid(j,:)+ delta*X_mat(:,:,j)+ beta*Qt_current;
                    %Qt_new = Qbar_Anne + alpha*Stdresid(j,:)'*Stdresid(j,:)+ delta*X_mat(:,:,j)+ beta*Qt_current;
                end
                
                Rt_mat(:,:,j) = Rt_DCC;
                
                Ht_DCC = D_mat*Rt_DCC*D_mat;
                [~,ind_chol] = chol(Ht_DCC);
                if ind_chol>0
                    loglik_vec = NaN;
                    break;
                end
                
                loglik_vec(j,1) = -0.5*log(det(Ht_DCC)) - 0.5 * Stdresid(j,:)*inv(Ht_DCC)*Stdresid(j,:)';
                
                Qt_current = Qt_new;
            end
            
            loglik_vec = loglik_vec  - 0.5*k*log(2*pi);
            LLF = -sum(loglik_vec);
            
            if isnan(LLF)
                LLF = 1e12;
            end
            
        end
        
        
        function [LLF,Rt_mat,Qbar_cDCC] = LogLik_cDCC_t(k,T,params, Stdresid)
            
            
            Qbar_DCC = cov(Stdresid);
            D_mat = diag(sqrt(diag(Qbar_DCC)));
            
            Rt_mat = [];
            %Rt_mat = zeros(k,k,T);
            Qt_current = Qbar_DCC;
            alpha = params(1);
            beta =  params(2);
            nu_t = params(3);
            
            loglik_r2 = zeros(T,1);
            
            q_ii_mat = zeros(T,k);
            q_ii_mat(1,:) = diag(Qbar_DCC);
            Qbar_new  = zeros(k,k);
            for m = 1:T
                if m<T
                    q_ii_mat(m+1,:) = (1 - alpha-beta) + alpha*q_ii_mat(m,:).*Stdresid(m,:).^2 + beta*q_ii_mat(m,:);
                end
                Q_square_root_star_t = diag(sqrt(q_ii_mat(m,:)));
                Qbar_new = Qbar_new +  Q_square_root_star_t * Stdresid(m,:)'* Stdresid(m,:) * Q_square_root_star_t;
            end
            
            Qbar_cDCC = Qbar_new/T;
            
            for j = 1:T
                
                Rt_DCC= Qt_current./(sqrt(diag(Qt_current))*sqrt(diag(Qt_current))');
                
                if j<T
                    temp = sqrt(diag(Qt_current));
                    Q_star_t = diag(temp(1:end));
                    Qt_new = Qbar_cDCC*(1 - alpha-beta) + alpha*(Q_star_t*Stdresid(j,:)'*Stdresid(j,:)*Q_star_t) + beta*Qt_current;
                end
                
                %Rt_mat(:,:,j) = Rt_DCC;
                
                Ht_DCC = D_mat*Rt_DCC*D_mat;
                
                [~,ind_chol] = chol(Ht_DCC);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                
                logdet_Ht = log(det(Ht_DCC));
                inv_Ht = Ht_DCC\eye(k);
                loglik_r2(j) = - 1/2 * logdet_Ht  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*Stdresid(j,:)*inv_Ht*Stdresid(j,:)');
                
                
                Qt_current = Qt_new;
            end
            
            loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);
            LLF = -sum(loglik_r2);
            
            if isnan(LLF)
                LLF = 1e12;
            end
            
        end
        
        
        function [LLF,Rt_mat,Qbar_cDCC] = LogLik_cDCC_TRIESZ(k,T,params, Stdresid)
            
            % not the true loglikelihood
            % only part that corresponds with the correlation part
            
            Qbar_DCC = cov(Stdresid);
            D_mat = diag(sqrt(diag(Qbar_DCC)));
            
            Rt_mat = [];
            %Rt_mat = zeros(k,k,T);
            Qt_current = Qbar_DCC;
            alpha = params(1);
            beta =  params(2);
            
            nu_vec = params(3:end);
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            
            loglik_r2 = zeros(T,1);
            
            q_ii_mat = zeros(T,k);
            q_ii_mat(1,:) = diag(Qbar_DCC);
            Qbar_new  = zeros(k,k);
            for m = 1:T
                if m<T
                    q_ii_mat(m+1,:) = (1 - alpha-beta) + alpha*q_ii_mat(m,:).*Stdresid(m,:).^2 + beta*q_ii_mat(m,:);
                end
                Q_square_root_star_t = diag(sqrt(q_ii_mat(m,:)));
                Qbar_new = Qbar_new +  Q_square_root_star_t * Stdresid(m,:)'* Stdresid(m,:) * Q_square_root_star_t;
            end
            
            Qbar_cDCC = Qbar_new/T;
            
            for j = 1:T
                
                Rt_DCC= Qt_current./(sqrt(diag(Qt_current))*sqrt(diag(Qt_current))');
                
                if j<T
                    temp = sqrt(diag(Qt_current));
                    Q_star_t = diag(temp(1:end));
                    Qt_new = Qbar_cDCC*(1 - alpha-beta) + alpha*(Q_star_t*Stdresid(j,:)'*Stdresid(j,:)*Q_star_t) + beta*Qt_current;
                end
                
                %Rt_mat(:,:,j) = Rt_DCC;
                
                Ht_DCC = D_mat*Rt_DCC*D_mat;
                [~,ind_chol] = chol(Ht_DCC);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                
                H_inv_j_0    = Ht_DCC\eye(k);
                
                [L_0_trans,ind_chol] = chol(H_inv_j_0);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                inv_H_j = L_trans' * L_trans;
                
                [~,ind_chol] = chol(inv_H_j);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                
                j_H = inv(inv_H_j);
                j_yy       = Stdresid(j,:)'*Stdresid(j,:);
                j_yy_H = j_yy + j_H;
                
                ln_PWD_H = LogLik.ln_PWD(inv_H_j,0.5*nu_vec);
                j_yy_H_inv = j_yy_H\eye(k);
                
                [~,ind_chol] = chol(j_yy_H_inv);
                if ind_chol>0
                    loglik_r2= NaN;
                    break
                end
                
                ln_PWD_j_yy_H  = LogLik.ln_PWD(j_yy_H_inv,0.5*(nu_vec+1));
                loglik_r2(j) = -ln_PWD_H  + ln_PWD_j_yy_H;
                
                Qt_current = Qt_new;
                
            end
            
            loglik_r2 = loglik_r2 - (k/2)*log(pi) + LogLik.lnmultigamma_nu_vec(k, (nu_vec+1)/2 )   -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
            LLF = -sum(loglik_r2);
            
            if isnan(LLF)
                LLF = 1e12;
            end
            
        end
        
        function [LLF,Rt_mat,Qbar_cDCC] = LogLik_cDCC_TRIESZ_INDUSTRY(k,T,params, Stdresid,n_vec)
            
            % not the true loglikelihood
            % only part that corresponds with the correlation part
            
            nr_groups = length(n_vec);
            
            nu_vec = [];
            for g = 1:nr_groups
                nu_vec = [nu_vec; params(g+2)*ones(n_vec(g),1)];
            end
            
            
            
            Qbar_DCC = cov(Stdresid);
            D_mat = diag(sqrt(diag(Qbar_DCC)));
            
            Rt_mat = [];
            %Rt_mat = zeros(k,k,T);
            Qt_current = Qbar_DCC;
            alpha = params(1);
            beta =  params(2);
            
            %nu_vec = params(3:end);
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            
            loglik_r2 = zeros(T,1);
            
            q_ii_mat = zeros(T,k);
            q_ii_mat(1,:) = diag(Qbar_DCC);
            Qbar_new  = zeros(k,k);
            for m = 1:T
                if m<T
                    q_ii_mat(m+1,:) = (1 - alpha-beta) + alpha*q_ii_mat(m,:).*Stdresid(m,:).^2 + beta*q_ii_mat(m,:);
                end
                Q_square_root_star_t = diag(sqrt(q_ii_mat(m,:)));
                Qbar_new = Qbar_new +  Q_square_root_star_t * Stdresid(m,:)'* Stdresid(m,:) * Q_square_root_star_t;
            end
            
            Qbar_cDCC = Qbar_new/T;
            
            for j = 1:T
                
                Rt_DCC= Qt_current./(sqrt(diag(Qt_current))*sqrt(diag(Qt_current))');
                
                if j<T
                    temp = sqrt(diag(Qt_current));
                    Q_star_t = diag(temp(1:end));
                    Qt_new = Qbar_cDCC*(1 - alpha-beta) + alpha*(Q_star_t*Stdresid(j,:)'*Stdresid(j,:)*Q_star_t) + beta*Qt_current;
                end
                
                %Rt_mat(:,:,j) = Rt_DCC;
                
                Ht_DCC = D_mat*Rt_DCC*D_mat;
                [~,ind_chol] = chol(Ht_DCC);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                
                H_inv_j_0    = Ht_DCC\eye(k);
                
                [L_0_trans,ind_chol] = chol(H_inv_j_0);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                inv_H_j = L_trans' * L_trans;
                
                [~,ind_chol] = chol(inv_H_j);
                if ind_chol>0
                    loglik_r2 = NaN;
                    break;
                end
                
                
                j_H = inv(inv_H_j);
                j_yy       = Stdresid(j,:)'*Stdresid(j,:);
                j_yy_H = j_yy + j_H;
                
                ln_PWD_H = LogLik.ln_PWD(inv_H_j,0.5*nu_vec);
                j_yy_H_inv = j_yy_H\eye(k);
                
                [~,ind_chol] = chol(j_yy_H_inv);
                if ind_chol>0
                    loglik_r2= NaN;
                    break
                end
                
                ln_PWD_j_yy_H  = LogLik.ln_PWD(j_yy_H_inv,0.5*(nu_vec+1));
                loglik_r2(j) = -ln_PWD_H  + ln_PWD_j_yy_H;
                
                Qt_current = Qt_new;
                
            end
            
            loglik_r2 = loglik_r2 - (k/2)*log(pi) + LogLik.lnmultigamma_nu_vec(k, (nu_vec+1)/2 )   -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
            LLF = -sum(loglik_r2);
            
            if isnan(LLF)
                LLF = 1e12;
            end
            
        end
        
        function [LLF,loglike_vec,Rt_mat,Qbar_cDCC,sigma_vec] = LogLik_cDCC_Copula(k,T,params, u_mat,ind_t_dist,ind_Rt)
            
            if ind_t_dist ==1
                nu = params(end);
                Stdresid_0 = tinv(u_mat,nu);
                Stdresid = sqrt((nu-2)/nu) * Stdresid_0;
            else
                Stdresid = norminv(u_mat);
            end
            
            if ind_Rt==1
                Rt_mat = zeros(k,k,T);
            else
                Rt_mat = [];
            end
            
            Qbar_DCC = cov(Stdresid);
            
            
            Qt_current = Qbar_DCC;
            alpha = params(1);
            beta =  params(2);
            
            loglike_vec = zeros(T,1);
            sigma_vec = zeros(T,1);
            
            q_ii_mat = zeros(T,k);
            q_ii_mat(1,:) = diag(Qbar_DCC);
            Qbar_new  = zeros(k,k);
            for m = 1:T
                if m<T
                    q_ii_mat(m+1,:) = (1-alpha-beta) + alpha*q_ii_mat(m,:).*Stdresid(m,:).^2 + beta*q_ii_mat(m,:);
                end
                Q_square_root_star_t = diag(sqrt(q_ii_mat(m,:)));
                Qbar_new = Qbar_new +  Q_square_root_star_t * Stdresid(m,:)'* Stdresid(m,:) * Q_square_root_star_t;
            end
            
            Qbar_cDCC = Qbar_new/T;
            
            for j = 1:T
                
                Rt_DCC= Qt_current./(sqrt(diag(Qt_current))*sqrt(diag(Qt_current))');
                
                if isnan(rcond(Rt_DCC))
                    loglike_vec = 1e14;
                    break;
                end
                
                if ind_Rt ==1
                    Rt_mat(:,:,j) =  Rt_DCC; %[Rt_DCC(1,3) Rt_DCC(2,9) Rt_DCC(5,15) Rt_DCC(6,14) Rt_DCC(12,18) Rt_DCC(17,19) Rt_DCC(8,25)];
                end
                if j<T
                    temp = sqrt(diag(Qt_current));
                    Q_star_t = diag(temp(1:end));
                    Qt_new = Qbar_cDCC*(1 - alpha-beta) + alpha*(Q_star_t*Stdresid(j,:)'*Stdresid(j,:)*Q_star_t) + beta*Qt_current;
                end
                
                sigma_vec(j,1) = sqrt((1/k)*ones(1,k)*Rt_DCC * ones(k,1)*(1/k));
                
                
                %Rt_mat(:,:,j) = Rt_DCC;
                x_R_inv_x = Stdresid(j,:)*inv(Rt_DCC)*Stdresid(j,:)';
                log_det_Rt   = log(det(Rt_DCC));
                if ind_t_dist ==1
                    % tdist
                    loglike_vec(j,1) = -0.5 *(nu+k)*log(1 + x_R_inv_x/(nu-2)) - 0.5* log_det_Rt;
                else
                    loglike_vec(j,1) = -0.5 * x_R_inv_x - 0.5*log_det_Rt;
                end
                
                Qt_current = Qt_new;
            end
            
            if isreal(loglike_vec) && sum(loglike_vec) ~=1e14
                
                if ind_t_dist==1
                    loglike_vec = loglike_vec + gammaln(0.5*(nu+k)) + (k-1)*gammaln(nu/2) - k * gammaln(0.5*(nu+1)) ...
                        + 0.5*(nu+1)*sum(log(ones(T,k)+(1/(nu-2))*Stdresid.^2),2);
                else
                    loglike_vec = loglike_vec + 0.5*sum(Stdresid.^2,2);
                end
                LLF = -sum(loglike_vec);
            else
                LLF = 1e14;
            end
            
        end
        
        
        function [LLF,loglik_RV] = LogLik_Uni_F_RV_T_r2(T,params,r2,RV);
            
            omega = params(1);
            alpha = params(2);
            beta  = params(3);
            nu_t  = params(4);
            nu_f1 = params(5);
            nu_f2 = params(6);
            
            constant = nu_f1/ (nu_f2-2);
            [h]  = Filter_New.GAS_Uni(T, [omega alpha beta nu_t nu_f1 nu_f2], r2, RV );
            
            
            loglik_RV = (0.5*nu_f1 - 1) * log(RV) -0.5*nu_f1 * log(h) -0.5*(nu_f1+nu_f2)*log(1 + (1./h)*constant.*RV) ...
                + 0.5*nu_f1 * log(constant);
            loglik_RV = loglik_RV + gammaln(0.5*(nu_f1+nu_f2)) - gammaln(0.5*nu_f1) - gammaln(0.5*nu_f2);
            
            loglik_r2 = -0.5 * log(h) - 0.5*(nu_t+1)*log(1 + r2./((nu_t-2)*h));
            loglik_r2 = loglik_r2 + gammaln(0.5*(nu_t+1)) - gammaln(0.5*nu_t) - 0.5 * log(pi*(nu_t-2));
            
            
            LLF= -sum(loglik_r2 + loglik_RV);
            
            if  isnan(LLF)
                LLF=1e6;
            end
            
            
        end
        
        function [LLF,loglik_RV] = LogLik_Uni_VT_F_RV_T_r2(T,params,r2,RV);
            
            alpha = params(1);
            beta  = params(2);
            nu_t  = params(3);
            nu_f1 = params(4);
            nu_f2 = params(5);
            
            constant = nu_f1/ (nu_f2-2);
            [h]  = Filter_New.GAS_VT_Uni(T, [alpha beta nu_t nu_f1 nu_f2], r2, RV );
            
            
            loglik_RV = (0.5*nu_f1 - 1) * log(RV) -0.5*nu_f1 * log(h) -0.5*(nu_f1+nu_f2)*log(1 + (1./h)*constant.*RV) ...
                + 0.5*nu_f1 * log(constant);
            loglik_RV = loglik_RV + gammaln(0.5*(nu_f1+nu_f2)) - gammaln(0.5*nu_f1) - gammaln(0.5*nu_f2);
            
            loglik_r2 = -0.5 * log(h) - 0.5*(nu_t+1)*log(1 + r2./((nu_t-2)*h));
            loglik_r2 = loglik_r2 + gammaln(0.5*(nu_t+1)) - gammaln(0.5*nu_t) - 0.5 * log(pi*(nu_t-2));
            
            
            LLF= -sum(loglik_r2 + loglik_RV);
            
            if  isnan(LLF)
                LLF=1e6;
            end
            
            
        end
        
        
        function [LLF,h] = LogLik_Heavy_N(T,params,r2,RV);
            
            omega = params(1);
            alpha = params(2);
            beta  = params(3);
            
            h = zeros(T,1);
            h(1) = mean(r2);
            for j = 2:T
                h(j,1) = omega + alpha * RV(j-1) + beta * h(j-1,1);
            end
            
            loglik_r2 = -0.5*log(2*pi)-0.5*(log(h) + r2./h);
            
            LLF= -sum(loglik_r2);
            
            if  isnan(LLF)
                LLF=1e6;
            end
            
            
        end
        
        
        
        function [LLF,loglik_F] = LogLik_Uni_F(T,params,RV);
            
            nu1 = params(1);
            nu2 = params(2);
            sigma2 = params(3);
            
            
            loglik_F= zeros(T,1);
            
            for j =1 :T
                %loglik_F(j) = (0.5*nu1 - 1) * log(X(j)) -0.5*(nu1+nu2)*log(1 + (nu1/nu2)*X(j));
                loglik_F(j) = (0.5*nu1 - 1) * log(RV(j)) -0.5*(nu1+nu2)*log(1 + (1/sigma2)*(nu1/(nu2-2))*RV(j));
                
            end
            
            loglik_F = loglik_F + gammaln(0.5*(nu1+nu2)) - gammaln(0.5*nu1) - gammaln(0.5*nu2);
            %loglik_F = loglik_F + 0.5*nu1 * log(nu1/nu2);
            loglik_F = loglik_F + 0.5*nu1 * log(nu1/(nu2-2)) - (nu1/2)*log(sigma2);
            
            LLF = -sum(loglik_F);
            
            if isnan(LLF)
                LLF = 1e12;
                
            end
            
            
            
            
        end
        
        
        function [min_LL_sum,h_mat,eps_mat] = Loglik_GARCH_t_aggregated(T,params,return_mat);
            
            k = size(return_mat,2);
            
            c = params(1);
            omega = params(2);
            alpha = params(3);
            beta  = params(4);
            nu_vec = params(5:end);
            
            min_LL_vec = zeros(k,1);
            h_mat = zeros(T,k);
            eps_mat = zeros(T,k);
            
            for i = 1:k
                [min_LL_i,h_vec_i,eps_vec_i] =  LogLik.Loglik_GARCH_t(T,[c omega alpha beta nu_vec(i)]',return_mat(:,i),0);
                min_LL_vec(i) = min_LL_i;
                h_mat(:,i) = h_vec_i;
                eps_mat(:,i) = eps_vec_i;
            end
            
            min_LL_sum = sum(min_LL_vec);
            
        end
        
        function [min_LL,h_vec,eps_vec] =  Loglik_GARCH_t(T,params,ret_vec,nr_lags)
            
            % garch (verschillende opties mogelijk) model met scheve t verdeling
            c     = params(1);
            omega = params(end-3);
            alpha = params(end-2);
            beta  = params(end-1);
            
            nu_t = params(end);
            
            if nr_lags>0
                lag_matrix = lagmatrix(ret_vec,1:nr_lags);
                lag_matrix(1:nr_lags,:) = [];
                phi_vec = params(2:2+nr_lags-1);
                eps_vec = ret_vec(nr_lags+1:end) - c - lag_matrix * phi_vec;
                T = T-nr_lags;
            else
                eps_vec = ret_vec - c;
            end
            
            eps2_vec =eps_vec.^2;
            h_vec = zeros(T,1);
            h_vec(1) = var(ret_vec);
            for i = 1:T+1
                if i<T
                    h_vec(i+1) = omega  + alpha * eps2_vec(i) + beta*h_vec(i);
                end
            end
            % h_last = h_vec(end);
            
            loglik_r2 = -0.5 * log(h_vec) - 0.5*(nu_t+1)*log(1 + eps2_vec./((nu_t-2)*h_vec));
            loglik_r2 = loglik_r2 + gammaln(0.5*(nu_t+1)) - gammaln(0.5*nu_t) - 0.5 * log(pi*(nu_t-2));
            
            min_LL= -sum(loglik_r2);
            
            if  isnan(min_LL) || isreal(loglik_r2)==0
                min_LL=1e6;
            end
            
        end
        
        function [LLF,loglik_C] = LogLik_Uni_Chi(T,params,RV);
            
            nu1 = params(1);
            sigma2 = params(2);
            
            
            loglik_C= zeros(T,1);
            
            for j =1 :T
                loglik_C(j) = (0.5*nu1 - 1) * log(RV(j)) -0.5*nu1*RV(j)/sigma2;
                
            end
            
            loglik_C = loglik_C - gammaln(0.5*nu1) - 0.5*nu1 * log(2/nu1) - (nu1/2)*log(sigma2);
            %loglik_F = loglik_F + 0.5*nu1 * log(nu1/nu2);
            
            LLF = -sum(loglik_C);
            
            if isnan(LLF)
                LLF = 1e12;
                
            end
            
            
            
            
        end
        
        
        
        function log_pdf = log_pdf_Matrix_F(RC,V,nu_vec_F);
            
            import LogLik.*
            
            nu_f1    = nu_vec_F(1);
            nu_f2    = nu_vec_F(2);
            k = size(RC,1);
            T = size(RC,3);
            log_pdf = zeros(T,1);
            
            constant = nu_f1/ (nu_f2-k-1);
            
            for i = 1:T
                RC_i       = RC(:,:,i);
                if size(V,3)==1
                    V_i        = V;
                else
                    V_i        = V(:,:,i);
                end
                rcond_i    = rcond(V_i);
                if rcond_i < 1e-15
                    apie = 1;
                end
                invV_i     = inv(V_i);
                logdetRC_i = log(det(RC_i));
                
                log_pdf(i,1) = 0.5*nu_f1 * log(det(constant*invV_i)) + 0.5*(nu_f1-k-1)*logdetRC_i - 0.5*(nu_f1+nu_f2)*log(det(eye(k)+ constant*invV_i*RC_i)) ...
                    + LogLik.lnmultigamma(k, 0.5*(nu_f1+nu_f2)) - LogLik.lnmultigamma(k, 0.5*nu_f1) - LogLik.lnmultigamma(k, 0.5*nu_f2);
                
            end
        end
        
        function log_pdf = log_pdf_Wishart(RC,V,nu_W);
            
            import LogLik.*
            k = size(RC,1);
            T = size(RC,3);
            log_pdf = zeros(T,1);
            
            for i = 1:T
                
                if size(V,3)==1
                    V_i = V;
                else
                    V_i = V(:,:,i);
                end
                RC_i = RC(:,:,i);
                invV_i     = inv(V_i);
                logdetV_i  = log(det(V_i));
                logdetRC_i = log(det(RC_i));
                
                log_pdf(i,1) =  - nu_W/2 * logdetV_i   +  (nu_W-k-1)/2 * logdetRC_i  - (nu_W/2) * trace(invV_i * RC_i) + ...
                    (k*nu_W)/2*log(nu_W/2) -  LogLik.lnmultigamma(k, nu_W/2 );
                
                
            end
        end
        
        function [log_pdf,C_q_mat] = log_pdf_Mult_T(return_mat,V_vech,nu_t,ind_csl,q_vec,M);
            
            C_q_mat = [];
            
            
            [T,k] = size(return_mat);
            I_k = eye(k);
            log_pdf = zeros(T,1);
            return_mat_dem = return_mat - ones(T,1)*mean(return_mat);
            
            if ind_csl==1
                C_q_mat = nan*ones(T,length(q_vec));
            end
            
            
            for i = 1:T
                V_i       = Admin.Vech2Sym(k,V_vech(i,:));
                invV_i    = V_i\I_k;
                logdetV_i = log(det(V_i));
                
                log_pdf(i,1) = gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi) ...
                    - 1/2 * logdetV_i  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*(return_mat_dem(i,:)*invV_i*return_mat_dem(i,:)'));
                
                
                if ind_csl==1
                    for m = 1:length(q_vec)
                        q_i = q_vec(m)*ones(k,1);
                        if sum(return_mat(i,:)<=q_i)==k
                            C_q_mat(i,m) = qsimvt(M,nu_t,V_i,-inf*ones(k,1),q_i);
                        end
                        
                    end
                end
            end
        end
        
        
        function log_pdf = log_pdf_Mult_TRIESZ(return_mat,vech_V,nu_vec);
            
            [T,k] = size(return_mat);
            log_pdf = zeros(T,1);
            y_mat = return_mat - ones(T,1)*mean(return_mat);
            
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            I_k = eye(k);
            
            for j = 1:T
                
                j_Sigma_0        = Admin.Vech2Sym(k,vech_V(j,:));
                Sigma_inv_j_0    = j_Sigma_0\I_k;
                
                [L_0_trans,ind_chol] = chol(Sigma_inv_j_0);
                %                     if ind_chol>0
                %                         loglik_r2 = NaN;
                %                         break;
                %                     end
                
                L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                inv_Sigma_j = L_trans' * L_trans;
                
                %                     [~,ind_chol] = chol(inv_Sigma_j);
                %                     if ind_chol>0
                %                         loglik_r2 = NaN;
                %                         break;
                %                     end
                %
                %
                j_Sigma = inv(inv_Sigma_j);
                
                j_yy       = y_mat(j,:)'*y_mat(j,:);
                j_yy_Sigma = j_yy + j_Sigma;
                
                ln_PWD_Sigma = LogLik.ln_PWD(inv_Sigma_j,0.5*nu_vec);
                j_yy_Sigma_inv = j_yy_Sigma\I_k;
                
                %                     [~,ind_chol] = chol(j_yy_Sigma_inv);
                %                     if ind_chol>0
                %                         loglik_r2= NaN;
                %                         break
                %                     end
                
                %ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma,-0.5*(nu_vec+1));
                ln_PWD_j_yy_Sigma  = LogLik.ln_PWD(j_yy_Sigma_inv,0.5*(nu_vec+1));
                
                log_pdf(j) = -ln_PWD_Sigma  + ln_PWD_j_yy_Sigma;
                
                
            end
            
            log_pdf = log_pdf - (k/2)*log(pi) + LogLik.lnmultigamma_nu_vec(k, (nu_vec+1)/2 )   -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
            %log_pdf = -log_pdf;
            
        end
        
        function [log_pdf,C_q_mat] = log_pdf_Mult_N(return_mat,V_vech,ind_csl,q_vec,M);
            
            C_q_mat = [];
            
            
            [T,k] = size(return_mat);
            log_pdf = zeros(T,1);
            return_mat_dem = return_mat - ones(T,1)*mean(return_mat);
            I_k = eye(k);
            
            if ind_csl==1
                C_q_mat = nan*ones(T,length(q_vec));
            end
            
            for i = 1:T
                V_i       = Admin.Vech2Sym(k,V_vech(i,:));
                invV_i    = V_i\I_k;
                logdetV_i = log(det(V_i));
                
                log_pdf(i,1) = - k/2 * log(2*pi)- 1/2 * logdetV_i - 1/2 * (return_mat_dem(i,:)*invV_i*return_mat_dem(i,:)');
                
                if ind_csl==1
                    for m = 1:length(q_vec)
                        q_i = q_vec(m)*ones(k,1);
                        if sum(return_mat(i,:)<=q_i)==k
                            C_q_mat(i,m) = qsimvn(M,V_i,-inf*ones(k,1),q_i);
                        end
                        
                    end
                end
                
                
            end
        end
        
        function log_pdf = log_pdf_IW(RC,V,nu_IW);
            k = size(RC,1);
            T = size(RC,3);
            
            log_pdf = zeros(T,1);
            constant = (nu_IW-k-1);
            for i = 1:T
                V_i       = V(:,:,i);
                logdetV_i = log(det(V_i));
                
                RC_i       = RC(:,:,i);
                logdetRC_i = log(det(RC_i));
                inv_RC_i   = inv(RC_i);
                
                log_pdf(i) = (nu_IW/2) * logdetV_i -0.5*(nu_IW+k+1)*logdetRC_i - constant/2 * trace(inv_RC_i*V_i)...
                    + 0.5*nu_IW*k*log(constant/2) - LogLik.lnmultigamma(k, 0.5*nu_IW);
            end
        end
        
        
    end
end