classdef LogLik_UV % in LogLik.m
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
        
        
        function [logL,Rt_mat,rho_vec,sigma_vec]=LogLik_DECO(k,T,params, Stdresid_0,asset_group_vec)
            
            %Stdresid_0 = norminv(u_mat);
            
            if isempty(asset_group_vec)
                g = 1;
                Rt_mat = [];
                rho_vec = zeros(T,1);
                Stdresid = Stdresid_0;
            else
                g = max(asset_group_vec);
                %nr_block_pairs = g * (g-1)/2;
                Rt_mat = zeros(k,k,T);
                rho_vec = [];
                
                % sort assets
                asset_group_mat = [asset_group_vec [1:k]'];
                asset_group_mat_sort = sortrows(asset_group_mat,1);
                Stdresid = Stdresid_0(:,asset_group_mat_sort(:,2));
                n_vec = sum((asset_group_vec*ones(1,g))== (ones(k,1)*(1:g)));
            end
            
            alpha = params(1);
            beta = params(2);
            iota = ones(k,1);
            I_n = eye(k);
            J_n = ones(k,k);
            logL = 0;
            
            Qbar_DCC = cov(Stdresid);
            Qt_current = Qbar_DCC;
            q_ii_mat = zeros(T,k);
            q_ii_mat(1,:) = diag(Qbar_DCC);
            Qbar_new  = zeros(k,k);
            
            sigma_vec = []; %zeros(T,1);
            
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
                    %Rt_DCC_new= Qt_new./(sqrt(diag(Qt_new))*sqrt(diag(Qt_new))');
                end
                
                if g==1
                    rho_t = (iota'*Rt_DCC*iota - k)/(k*(k-1));
                    R_inv_t = 1/(1-rho_t) * I_n - rho_t/((1-rho_t)*(1+(k-1)*rho_t)) * J_n;
                    det_R_t = ((1-rho_t)^(k-1))*(1+(k-1)*rho_t);
                    
                    rho_vec(j,1) = rho_t;
                    logL = logL + log(det_R_t) + Stdresid(j,:)*R_inv_t*Stdresid(j,:)' - Stdresid(j,:)*Stdresid(j,:)';
                    
                    %rho_mat = rho_t*ones(n,n)-eye(n);
                    %sigma_vec(j,1) = sqrt((1/n)*ones(1,n)*rho_mat * ones(n,1)*(1/n));
                    
                else
                    % build r_ij
                    teller_i = 1;
                    Rho_block_mat = zeros(g,g);
                    ingred_mat_lik = zeros(g,3);
                    Rt_0 = zeros(k,k);
                    Rt_1 = zeros(k,k);
                    for i = 1:g
                        n_i = n_vec(i);
                        Rt_DCC_ii = Rt_DCC(teller_i:teller_i+n_i-1,teller_i:teller_i+n_i-1);
                        if n_i>1
                            rho_ii = (ones(n_i,1)'*Rt_DCC_ii*ones(n_i,1) - n_i)/(n_i*(n_i-1));
                        else
                            rho_ii = 0;
                        end
                        Rho_block_mat(i,i) = rho_ii;
                        
                        Rt_0(teller_i:teller_i+n_i-1,teller_i:teller_i+n_i-1) = (1- rho_ii)*eye(n_i);
                        Rt_1(teller_i:teller_i+n_i-1,teller_i:teller_i+n_i-1) = rho_ii;
                        
                        teller_j = sum(n_vec(1:i))+1;
                        for b = i+1:g
                            Rt_DCC_ij = Rt_DCC(teller_i:teller_i+n_i-1,teller_j:teller_j+n_vec(b)-1);
                            rho_ij = mean(Rt_DCC_ij(:));
                            Rho_block_mat(i,b) = rho_ij;
                            
                            Rt_1(teller_i:teller_i+n_i-1,teller_j:teller_j+n_vec(b)-1) = rho_ij;
                            Rt_1(teller_j:teller_j+n_vec(b)-1,teller_i:teller_i+n_i-1) = rho_ij;
                            
                            teller_j = teller_j+n_vec(b);
                        end
                        
                        r_i_vec = Stdresid(j,teller_i:teller_i+n_i-1);
                        ingred_mat_lik(i,:) = [sum(r_i_vec.^2) sum(r_i_vec)^2 sum(r_i_vec)];
                        
                        teller_i = teller_i+n_vec(i);
                    end
                    
                    Rt_DECO_gr = Rt_0 + Rt_1;
                    
                    %Rt_mat(:,:,j) = Rt_DECO_gr;
                    
                    
                    % now compute the likelihood;
                    % ga alle pairs af
                    % let o
                    
                    b_vec = 1./(1-diag(Rho_block_mat));
                    for i = 1:g
                        b_i      = b_vec(i);
                        ingred_i = ingred_mat_lik(i,:);
                        rho_ii   = Rho_block_mat(i,i);
                        n_i      = n_vec(i);
                        for b = i+1:g
                            rho_bb = Rho_block_mat(b,b);
                            rho_ib = Rho_block_mat(i,b);
                            b_b    = b_vec(b);
                            n_b    = n_vec(b);
                            ingred_b = ingred_mat_lik(b,:);
                            
                            temp_c = (rho_ii*(n_i-1)+1)*(rho_bb*(n_b-1)+1) - n_i*n_b*rho_ib^2;
                            
                            c_1    = (rho_ii*(rho_bb*(n_b-1)+1) - rho_ib^2*n_b)/((rho_ii-1)*temp_c);
                            c_2    = (rho_bb*(rho_ii*(n_i-1)+1) - rho_ib^2*n_i)/((rho_bb-1)*temp_c);
                            c_3    = rho_ib/(n_i*n_b*rho_ib^2 - (rho_ii*(n_i-1)+1)*(rho_bb*(n_b-1)+1));
                            
                            logL   = logL + (n_i-1)*log(1-rho_ii) + (n_b-1)*log(1-rho_bb) + log(temp_c) + b_i*ingred_i(1) + b_b*ingred_b(1) + ...
                                c_1 * ingred_i(2) + c_2 * ingred_b(2) + 2*c_3*ingred_i(3)*ingred_b(3) - Stdresid(j,:)*Stdresid(j,:)';
                        end
                    end
                    
                    
                end
                
                Qt_current = Qt_new;
            end
            
            if isnan(logL) || isinf(abs(logL)) || isreal(logL)~=1
                logL = 1e14;
            else
                logL = (1/2)*logL;
            end
            
        end
        
        function [logL,Rt_mat,Qbar_cDCC] = LogLik_cDCC_N(k,T,params, Stdresid)
            
            Qbar_DCC = cov(Stdresid);
            Rt_mat = zeros(k,k,T);
            Qt_current = Qbar_DCC;
            alpha = params(1);
            beta =  params(2);
            
            logL = 0;
            
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
                logL = logL + log(det(Rt_DCC)) + Stdresid(j,:)*inv(Rt_DCC)*Stdresid(j,:)' - Stdresid(j,:)*Stdresid(j,:)';
                
                Qt_current = Qt_new;
            end
            
            %Q_T_plus_1 = Qt_current;
            %R_T_plus_1 = Rt_current;
            logL = (1/2)*logL;
        end
        
        
        function CL = LogLik_cDCC_CL(k,T,params,Stdresid)
            
            
            nr_comb = k-1;
            CL_vec = zeros(nr_comb,1);
            %tic;
            for i = 1:nr_comb
                Stdresid_pair_i = Stdresid(:,i:i+1);
                CL_vec(i,1) =  LogLik.LogLik_cDCC_N(2,T,params,Stdresid_pair_i);
            end
            %toc;
            
            CL = sum(CL_vec);
            
        end
        
        
        function LogL = LogLik_N_Pseudo_DECO_BLOCK(Stdresid_0,Rt_mat,asset_group_vec)
            
            %Stdresid_0 = norminv(u_mat);
            [T,k] = size(Stdresid_0);
            
            asset_group_mat = [asset_group_vec [1:k]'];
            asset_group_mat_sort = sortrows(asset_group_mat,1);
            Stdresid = Stdresid_0(:,asset_group_mat_sort(:,2));
            
            loglik_vec = zeros(T,1);
            for j = 1:T
                loglik_vec(j,1) = log(det(Rt_mat(:,:,j))) + Stdresid(j,:)*inv(Rt_mat(:,:,j))*Stdresid(j,:)'  - Stdresid(j,:)*Stdresid(j,:)';
            end
            
            LogL = -0.5*sum(loglik_vec);
            
            
        end
        
        
        
        function [LLF,LLF_r2,LLF_RC,V] = LogLik_GAS_tF_CT(k,T,params,r2,RC);
            
            
            alpha = params(1);
            beta  = params(2);
            nu_t  = params(3);
            nu_f1 = params(4);
            nu_f2 = params(5);
            
            constant = nu_f1/ (nu_f2-k-1);
            
            [V,~,Flag]  = Filter_New.GAS_CT( k, T, [alpha beta nu_t nu_f1 nu_f2], r2, RC );
            %[V,~,Flag]  = Filter_New.GAS_Gen_CT( k, T, [alpha_vec beta nu_t nu_f1 nu_f2], r2, RC );
            
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
                loglik_RC = loglik_RC + LogLik.lnmultigamma(k, 0.5*(nu_f1+nu_f2)) - LogLik.lnmultigamma(k, 0.5*nu_f1) - LogLik.lnmultigamma(k, 0.5*nu_f2);
                
                LLF_r2 = sum(loglik_r2);
                LLF_RC = sum(loglik_RC);
                
                LLF= -sum(loglik_r2 + loglik_RC);
                
                if  isnan(LLF)
                    LLF=1e6;
                end
                
            end
            
        end
        
        
        function [LLF,LLF_r2,V] = LogLik_Mult_t(k,T,params,r2_mat);
            
            
            nu_t  = params(1);
            
            V = mean(r2_mat,3);
                    j_invV       = inv(V);
                    j_logdetV    = log(det(V));                        
               loglik_r2 = zeros(T,1);
               
                for j=1:T
                    
                    %loglik_r2(j) = - 1/2 * j_logdetV  - 1/2 * trace( j_invV *r2(:,:,j) ); normal
                    
                    % Student T
                    loglik_r2(j) = - 1/2 * j_logdetV  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(j_invV*r2_mat(:,:,j)));
                    
                    
                end
                
                %loglik_r2= loglik_r2 - k/2 * log(2*pi); % normal
                
                loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);                               
                LLF_r2 = sum(loglik_r2);                                
                LLF= -sum(loglik_r2);
                
                if  isnan(LLF)
                    LLF=1e6;
                end                            
            
        end
        
        
        function [LLF,LLF_RC,V] = LogLik_GAS_F_MV(k,T,params,RC);
            
            alpha = params(1);
            beta  = params(2);
            %nu_t  = params(3);
            nu_f1 = params(3);
            nu_f2 = params(4);
            
            V0  = mean(RC, 3);
            fbar= V0;
            [s, f, V] = Filter_New.Store(k, T, V0);
            constant = nu_f1/ (nu_f2-k-1);
            %loglik_r2= zeros(1,T);
            loglik_RC= zeros(1,T);
            
            for j=1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    LLF = 1e12;
                    %LLF_r2 = 1e12;
                    LLF_RC = 1e12;
                    break
                end
                
                [~,nr_pd] = chol(V(:,:,j));
                if nr_pd>0
                    LLF = 1e12;
                    LLF_RC = 1e12;
                    break
                end
                
                % update the score and cov
                V_j   = V(:,:,j);
                j_RC  = RC(:,:,j);
                j_invV  = inv( V_j );
                j_logdetRC   = log(det(j_RC));
                %w_j   = (nu_t + k)/ (nu_t - 2 + trace(invV*r2(:,:,j)));
                %d1    =  1*w_j * r2(:,:,j) - 1*V_j;
                d2    = 1*nu_f1*((nu_f1 +nu_f2)/(nu_f2-k-1)*j_RC*inv(eye(k)+ constant*j_invV*j_RC) - V_j) ;
                
                %nabla = (d1 + d2)/(1+nu_f1);
                nabla = (d2)/(1+nu_f1);
                s(:,:,j) = nabla;
                
                if j<T
                    f(:,:,j+1) = (1-beta)*fbar + beta * f(:,:,j) + alpha * s(:,:,j);
                    V(:,:,j+1) = f(:,:,j+1);
                end
                
                loglik_RC(j) = 0.5*nu_f1 * log(det(constant*j_invV)) + 0.5*(nu_f1-k-1)*j_logdetRC - 0.5*(nu_f1+nu_f2)*log(det(eye(k)+ constant*j_invV*j_RC));
                
            end
            
            loglik_RC = loglik_RC + LogLik.lnmultigamma(k, 0.5*(nu_f1+nu_f2)) - LogLik.lnmultigamma(k, 0.5*nu_f1) - LogLik.lnmultigamma(k, 0.5*nu_f2);
            
            %LLF_r2 = sum(loglik_r2);
            LLF_RC = sum(loglik_RC);
            LLF= -LLF_RC;
            
            if  isnan(LLF)
                LLF=1e14;
            end
            
            
            
        end
        
        
        
        function [LLF,LLF_r2,V,V_N] = LogLik_GAS_t_MV_given_RK(k,T,params,V_DAY,r2_mat,r2_mat_N);
            
            %c  = params(1);
            alpha = params(end-2);
            beta  = params(end-1);
            nu_t  = params(end);
            
            
            V0 = mean(r2_mat,3);
            fbar= (1-beta)* mean(r2_mat_N,3);
            V_N = zeros(k,k,T);
            V = zeros(k,k,T);
            V_N(:,:,1) =mean(r2_mat_N, 3);
            V(:,:,1) = V0;
            
            loglik_r2= zeros(T,1);
            Flag = 0;
            for j=1:T
                
                check_ind = rcond(V(:,:,j));
                [~,nr_pd] = chol(V(:,:,j));
                if check_ind < 1e-15 || nr_pd>0
                    LLF = 1e12;
                    LLF_r2 = 1e12;
                    Flag = 1;
                    %LLF_RC = 1e12;
                    break
                else
                    
                    % update the score and cov
                    V_j   = V(:,:,j);
                    j_invV  = inv( V_j );
                    j_logdetV  = log(det(V_j));
                    
                    w_j   = (nu_t + k)/ (nu_t - 2 + trace(j_invV*r2_mat(:,:,j)));
                    d1    =  1*w_j * r2_mat(:,:,j) - 1*V_j;
                    s_j = d1;
                    
                    if j<T
                        V_N(:,:,j+1) = fbar + beta * V_N(:,:,j) + alpha * s_j;
                        V(:,:,j+1) = V_N(:,:,j+1) + V_DAY(:,:,j+1);
                    end
                    loglik_r2(j) = - 1/2 * j_logdetV  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(j_invV*r2_mat(:,:,j)));
                    
                    
                end
                
            end
            
            if Flag==0
                loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);
                
                LLF_r2 = sum(loglik_r2);
                %LLF_RC = sum(loglik_RC);
                LLF= -LLF_r2;
                
                if  isnan(LLF)
                    LLF=1e14;
                end
                
            end
            
        end
        
        
        function [LLF,LLF_r2,V] = LogLik_GAS_t_MV(k,T,params,r2_mat);
            
            %c  = params(1);
            alpha = params(1);
            beta  = params(2);
            nu_t  = params(3);
            
            
            V0 = mean(r2_mat,3);
            %fbar= c * eye(k);
            fbar = V0*(1-beta);
            V = zeros(k,k,T);
            V(:,:,1) = V0;
            FLAG=0;
            loglik_r2= zeros(T,1);
            
            for j=1:T
                
                check_ind = rcond(V(:,:,j));
                [~,nr_pd] = chol(V(:,:,j));
                if check_ind < 1e-15 || nr_pd>0
                    LLF = 1e12;
                    LLF_r2 = 1e12;
                    FLAG=1;
                    %LLF_RC = 1e12;
                    break
                else
                    
                    
                    % update the score and cov
                    V_j   = V(:,:,j);
                    j_invV  = inv( V_j );
                    j_logdetV  = log(det(V_j));
                    
                    w_j   = (nu_t + k)/ (nu_t - 2 + trace(j_invV*r2_mat(:,:,j)));
                    d1    =  1*w_j * r2_mat(:,:,j) - 1*V_j;
                    s_j = d1;
                    
                    if j<T
                        V(:,:,j+1) = fbar + beta * V(:,:,j) + alpha * s_j;
                    end
                    loglik_r2(j) = - 1/2 * j_logdetV  - 0.5*(nu_t + k) * log(1 + (1/(nu_t-2))*trace(j_invV*r2_mat(:,:,j)));
                    
                end
                
            end
               
            if FLAG~=1
                loglik_r2 = loglik_r2  + gammaln(0.5*(nu_t + k)) - gammaln(nu_t/2) - 0.5*k*log((nu_t-2)*pi);
                LLF_r2 = sum(loglik_r2);
                %LLF_RC = sum(loglik_RC);
                LLF= -LLF_r2;
                
                if  isnan(LLF)
                    LLF=1e14;
                end
                
            end
            
        end
        
        
        
        
        
        function [LLF,h_RV] = LogLik_GAS_F(T,params,RV,ON_return);
            % GAS on RV only (HENCE OPEN TO CLOSE VOL)
            
            %omega = (1-beta)*fbar
            omega = params(1);
            alpha = params(2);
            beta  = params(3);
            nu_f1 = params(4);
            nu_f2 = params(5);
            
            if size(ON_return,1)>0
                gamma = params(6);
            else
                gamma = 0;
                ON_return = zeros(T,1);
            end
            
            
            %RV_r2_factor = mean(eps2_vec)/mean(RV);
            %constant =  RV_r2_factor * nu_f1/ (nu_f2-2);
            constant =   nu_f1/ (nu_f2-2);
            % score objects: score and GAS
            
            h0 = mean(RV);
            h_RV = zeros(T,1);
            h_RV(1) = h0;
            s_RV = zeros(T,1);
            %fbar = h0;
            
            
            for j=1:T
                
                % update the score and cov
                
                d_RV = nu_f1 *(((nu_f1 +nu_f2)/(nu_f2-2))*RV(j,1)/(1 + constant * RV(j,1)/h_RV(j,1)) - h_RV(j,1));
                %d_RV = nu_f1 *(((nu_f1 +nu_f2)/(nu_f2-2))*RV(j,1)/(1 + constant * RV(j,1)/h_RV(j,1)));
                nabla = (d_RV)/(nu_f1+1);
                s_RV(j,1) = nabla;
                
                if j<T
                    h_RV(j+1,1) = omega + beta * h_RV(j,1) + alpha * s_RV(j,1) + gamma * ON_return(j+1,1)^2;
                end
                
            end
            
            loglik_RV = (0.5*nu_f1 - 1) * log(RV) -0.5*nu_f1 * log(h_RV) -0.5*(nu_f1+nu_f2)*log(1 + (1./h_RV)*constant.*RV) ...
                + 0.5*nu_f1 * log(constant);
            loglik_RV = loglik_RV +  gammaln(0.5*(nu_f1+nu_f2)) - gammaln(0.5*nu_f1) - gammaln(0.5*nu_f2);
            
            LLF= -sum(loglik_RV);
            
            if  isnan(LLF) || isreal(LLF)==0|| isreal(loglik_RV)==0 || any(h_RV<=0)
                LLF=1e14;
            end
            
            
        end
        
        
        function [LLF,h,c_vec,c_vec_tot,s_c,eps_vec] = LogLik_GAS_t_CTC_given_h_RV(T,params,ret_vec,h_OC,X_vec);
            % nr_lags
            
            c      = params(1);
            eps_vec = ret_vec - c;
            eps2_vec = eps_vec.^2;
            
            %omega = (1-beta)*fbar
            omega = params(2);
            alpha = params(3);
            beta  = params(4);
            nu_t  = params(5);
            
            if ~isempty(X_vec)
                gamma = params(6);
            else
                gamma = 0;
                X_vec = zeros(T,1);
            end
            
            RV_r2_factor = mean(eps2_vec)/mean(h_OC);
            % score objects: score and GAS
            
            c_vec =  RV_r2_factor *ones(T,1);
            c_vec_tot = c_vec;
            h = zeros(T,1);
            h(1) = mean(eps2_vec);
            s_c = zeros(T,1);
            
            for j=1:T
                
                % update the score and cov
                w_j = (nu_t + 1)/ (nu_t - 2 +  eps2_vec(j,1)/h(j,1));
                d_c = w_j *  eps2_vec(j,1)/h_OC(j,1) - c_vec(j,1);
                s_c(j,1) = d_c;
                
                if j<T
                    c_vec(j+1,1) = omega + beta * c_vec(j,1) + alpha * s_c(j,1);                    
                    c_vec_tot(j+1,1) = c_vec(j+1,1) + gamma* X_vec(j+1);
                    h(j+1,1) = h_OC(j+1,1)* c_vec_tot(j+1);
                end
                
            end
            
            % h = h_RV.* c_vec;
            
            loglik_r2 = -0.5 * log(h) - 0.5*(nu_t+1)*log(1 + eps2_vec./((nu_t-2)*h));
            loglik_r2 = loglik_r2 + gammaln(0.5*(nu_t+1)) - gammaln(0.5*nu_t) - 0.5 * log(pi*(nu_t-2));
            LLF= -sum(loglik_r2);
            
            if  isnan(LLF) || isreal(LLF)==0||any(h<=0)
                LLF=1e14;
            end
            
            
        end
        
        
        function [LLF,LLF_r2,LLF_RV,h,h_RV,c_vec] = LogLik_GAS_t_CTC_given_h_RV_1_step(T,params,ret_vec,RV);
            % nr_lags
            % ret_vec = close to close return
            
            c      = params(1);
            eps_vec = ret_vec - c;
            eps2_vec = eps_vec.^2;
            
            %omega = (1-beta)*fbar
            omega_c = params(2);
            alpha_c = params(3);
            beta_c  = params(4);
            nu_t  = params(5);
            omega_D = params(6);
            alpha_D = params(7);
            beta_D  = params(8);            
            nu_f1 = params(9);
            nu_f2 = params(10);
                        
            constant =   nu_f1/ (nu_f2-2);
            h0 = mean(RV);
            h_RV = zeros(T,1);
            h_RV(1) = h0;
            s_RV = zeros(T,1);
            
            RV_r2_factor = mean(eps2_vec)/h0;
            % score objects: score and GAS
            
            c_vec =  RV_r2_factor *ones(T,1);            
            h = zeros(T,1);
            h(1) = mean(eps2_vec);
            s_c = zeros(T,1);
            
            for j=1:T
                % update the score and cov
                
                d_RV = nu_f1 *(((nu_f1 +nu_f2)/(nu_f2-2))*RV(j,1)/(1 + constant * RV(j,1)/h_RV(j,1)) - h_RV(j,1));
                s_RV(j,1) = (d_RV)/(nu_f1+1);
                                                            
                w_j = (nu_t + 1)/ (nu_t - 2 +  eps2_vec(j,1)/h(j,1));
                d_c = w_j *  eps2_vec(j,1)/h_RV(j,1) - c_vec(j,1);
                s_c(j,1) = d_c;
                
                if j<T
                    h_RV(j+1,1) = omega_D + beta_D * h_RV(j,1) + alpha_D * s_RV(j,1);
                    c_vec(j+1,1) = omega_c + beta_c * c_vec(j,1) + alpha_c * s_c(j,1);                                                            
                    h(j+1,1) = h_RV(j+1,1)* c_vec(j+1);
                end
                
            end
                                                            
            loglik_r2 = -0.5 * log(h) - 0.5*(nu_t+1)*log(1 + eps2_vec./((nu_t-2)*h));
            loglik_r2 = loglik_r2 + gammaln(0.5*(nu_t+1)) - gammaln(0.5*nu_t) - 0.5 * log(pi*(nu_t-2));
            
            loglik_RV = (0.5*nu_f1 - 1) * log(RV) -0.5*nu_f1 * log(h_RV) -0.5*(nu_f1+nu_f2)*log(1 + (1./h_RV)*constant.*RV) ...
                + 0.5*nu_f1 * log(constant);
            loglik_RV = loglik_RV +  gammaln(0.5*(nu_f1+nu_f2)) - gammaln(0.5*nu_f1) - gammaln(0.5*nu_f2);
            
            LLF_r2 = sum(loglik_r2);
            LLF_RV = sum(loglik_RV);
            
            LLF= -sum(loglik_r2 + loglik_RV);
            
            if  isnan(LLF) || ~isreal(LLF) ||any(h<=0)
                LLF=1e14;
            end
            
            
        end
        
        
        
        function [LLF,h,c_vec,eps_vec] = LogLik_VIX_t_CTC_given_h_RV(T,params,ret_vec,h_OC,VIX_vec);
            % nr_lags
            
            c      = params(1);
            eps_vec = ret_vec - c;
            eps2_vec = eps_vec.^2;
            
            %omega = (1-beta)*fbar
            omega = params(2);          
            nu_t  = params(3);
            gamma = params(4);
                                                                       
            c_vec =  omega + gamma*VIX_vec;
            h     = h_OC.*c_vec;                                                            
            
            loglik_r2 = -0.5 * log(h) - 0.5*(nu_t+1)*log(1 + eps2_vec./((nu_t-2)*h));
            loglik_r2 = loglik_r2 + gammaln(0.5*(nu_t+1)) - gammaln(0.5*nu_t) - 0.5 * log(pi*(nu_t-2));
            LLF= -sum(loglik_r2);
            
            if  isnan(LLF) || isreal(LLF)==0||any(h<=0)
                LLF=1e14;
            end
            
            
        end
        
        
        
        function [LLF,h,h_N,s_N,eps_vec] = LogLik_GAS_t_CTC_given_h_RV_ADD(T,params,ret_vec,h_D,ret_N);
            % nr_lags
            
            c      = params(1);
            eps_vec = ret_vec - c;
            eps2_vec = eps_vec.^2;
            
            %omega = (1-beta)*fbar
            omega = params(2);
            alpha = params(3);
            beta  = params(4);
            nu_t  = params(5);
            
            
            % score objects: score and GAS
            
            h_N =  zeros(T,1);
            h_N(1) = var(ret_N);
            h = zeros(T,1);
            h(1) = mean(eps2_vec);
            s_N = zeros(T,1);
            
            for j=1:T
                
                % update the score and cov
                w_j = (nu_t + 1)/ (nu_t - 2 +  eps2_vec(j,1)/h(j,1));
                d_N = w_j *  eps2_vec(j,1) - h(j,1);
                s_N(j,1) = d_N;
                
                if j<T
                    h_N(j+1,1) = omega + beta * h_N(j,1) + alpha * s_N(j,1);                                                            
                    h(j+1,1) = h_D(j+1,1) + h_N(j+1,1);
                end
                
            end
            
            if any(h_N<0)
                LLF = 1e14;
            
            else
           
            
            % h = h_RV.* c_vec;
            
            loglik_r2 = -0.5 * log(h) - 0.5*(nu_t+1)*log(1 + eps2_vec./((nu_t-2)*h));
            loglik_r2 = loglik_r2 + gammaln(0.5*(nu_t+1)) - gammaln(0.5*nu_t) - 0.5 * log(pi*(nu_t-2));
            LLF= -sum(loglik_r2);
            
            if  isnan(LLF) || isreal(LLF)==0||any(h<=0)
                LLF=1e14;
            end
            
            end
        end
        
        
        
        function [LLF,LLF_r2,LLF_RV,h,h_N] = LogLik_GAS_t_CTC_ADD_1_step(T,params,ret_vec,ret_N,RV);
            % nr_lags
            % ret_vec = close to close return
            
            c      = params(1);
            eps_vec = ret_vec - c;
            eps2_vec = eps_vec.^2;
            
            %omega = (1-beta)*fbar
            omega_D = params(2);
            alpha_D = params(3);
            beta_D  = params(4);
            omega_N = params(5);
            alpha_N = params(6);
            beta_N  = params(7);
            nu_t  = params(8);
            nu_f1 = params(9);
            nu_f2 = params(10);
            
            
            constant =   nu_f1/ (nu_f2-2);
            h0 = mean(RV);
            h_RV = zeros(T,1);
            h_RV(1) = h0;
            s_RV = zeros(T,1);
            
            h_N =  zeros(T,1);
            h_N(1) = var(ret_N);
            h = zeros(T,1);
            h(1) = mean(eps2_vec);
            s_N = zeros(T,1);
            
            for j=1:T
                
                % update the score and cov
                w_j = (nu_t + 1)/ (nu_t - 2 +  eps2_vec(j,1)/h(j,1));
                d_N = w_j *  eps2_vec(j,1) - h(j,1);
                s_N(j,1) = d_N;
                
                d_RV = nu_f1 *(((nu_f1 +nu_f2)/(nu_f2-2))*RV(j,1)/(1 + constant * RV(j,1)/h_RV(j,1)) - h_RV(j,1));
                s_RV(j,1) = (d_RV)/(nu_f1+1);
                
                if j<T
                    h_N(j+1,1) = omega_N + beta_N * h_N(j,1) + alpha_N * s_N(j,1);
                    h_RV(j+1,1) = omega_D + beta_D * h_RV(j,1) + alpha_D * s_RV(j,1);
                    h(j+1,1) = h_RV(j+1,1) + h_N(j+1,1);
                end
                
                
            end
            
            % h = h_RV.* c_vec;
            
            loglik_r2 = -0.5 * log(h) - 0.5*(nu_t+1)*log(1 + eps2_vec./((nu_t-2)*h));
            loglik_r2 = loglik_r2 + gammaln(0.5*(nu_t+1)) - gammaln(0.5*nu_t) - 0.5 * log(pi*(nu_t-2));
            
            loglik_RV = (0.5*nu_f1 - 1) * log(RV) -0.5*nu_f1 * log(h_RV) -0.5*(nu_f1+nu_f2)*log(1 + (1./h_RV)*constant.*RV) ...
                + 0.5*nu_f1 * log(constant);
            loglik_RV = loglik_RV +  gammaln(0.5*(nu_f1+nu_f2)) - gammaln(0.5*nu_f1) - gammaln(0.5*nu_f2);
            
            LLF_r2 = sum(loglik_r2);
            LLF_RV = sum(loglik_RV);
            
            LLF= -sum(loglik_r2 + loglik_RV);
            
            if  isnan(LLF) || isreal(LLF)==0||any(h<=0)
                LLF=1e14;
            end
            
            
        end
        
        
        
        
        function [LLF,h,eps_vec] = LogLik_t_CTC_given_h_RV(T,params,ret_vec,h_OC,c_fix_given);
            % nr_lags
            
            c      = params(1);
            eps_vec = ret_vec - c;
            eps2_vec = eps_vec.^2;
            
            %omega = (1-beta)*fbar
            if length(params) == 3
                c_fixed = params(2);
                nu_t  = params(3);
            else
                c_fixed = c_fix_given;
                nu_t  = params(2);
            end
            
            %RV_r2_factor = mean(eps2_vec)/mean(h_OC);
            
            %c_vec =  RV_r2_factor *ones(T,1);
            h = h_OC * c_fixed;
            
            % h = h_RV.* c_vec;
            
            loglik_r2 = -0.5 * log(h) - 0.5*(nu_t+1)*log(1 + eps2_vec./((nu_t-2)*h));
            loglik_r2 = loglik_r2 + gammaln(0.5*(nu_t+1)) - gammaln(0.5*nu_t) - 0.5 * log(pi*(nu_t-2));
            LLF= -sum(loglik_r2);
            
            if  isnan(LLF) || isreal(LLF)==0||any(h<=0)
                LLF=1e14;
            end
            
            
        end
        
        
        
        
        function [LLF,LLF_t,LLF_RV,h,eps_vec] = LogLik_GAS_tF_VT_Uni(T,params,ret_vec,RV,nr_lags);
            % nr_lags
            
            c      = params(1);
            if nr_lags>0
                lag_matrix = lagmatrix(ret_vec,1:nr_lags);
                lag_matrix(1:nr_lags,:) = [];
                phi_vec = params(2:2+nr_lags-1);
                eps_vec = ret_vec(nr_lags+1:end) - c - lag_matrix * phi_vec;
                T = T-nr_lags;
                RV = RV(nr_lags+1:end);
            else
                eps_vec = ret_vec - c;
            end
            
            eps2_vec = eps_vec.^2;
            
            %omega = (1-beta)*fbar
            omega = params(end-5);
            alpha = params(end-4);
            beta  = params(end-3);
            nu_t  = params(end-2);
            nu_f1 = params(end-1);
            nu_f2 = params(end);
            
            %RV_r2_factor = mean(eps2_vec)/mean(RV);
            %constant =  RV_r2_factor * nu_f1/ (nu_f2-2);
            constant =   nu_f1/ (nu_f2-2);
            % score objects: score and GAS
            
            h0 = mean(RV);
            h = zeros(T,1);
            h(1) = h0;
            s = zeros(T,1);
            %fbar = h0;
            
            for j=1:T
                
                % update the score and cov
                w_j = (nu_t + 1)/ (nu_t - 2 +  eps2_vec(j,1)/h(j,1));
                d1 = w_j *  eps2_vec(j,1) - h(j,1);
                d2 = nu_f1 *(((nu_f1 +nu_f2)/(nu_f2-2))*RV(j,1)/(1 + constant * RV(j,1)/h(j,1)) - h(j,1));
                
                nabla = (d1 + d2)/(nu_f1+1);
                s(j,1) = nabla;
                
                if j<T
                    h(j+1,1) = omega + beta * h(j,1) + alpha * s(j,1);
                end
                
            end
            
            
            loglik_RV = (0.5*nu_f1 - 1) * log(RV) -0.5*nu_f1 * log(h) -0.5*(nu_f1+nu_f2)*log(1 + (1./h)*constant.*RV) ...
                + 0.5*nu_f1 * log(constant);
            loglik_RV = loglik_RV +  gammaln(0.5*(nu_f1+nu_f2)) - gammaln(0.5*nu_f1) - gammaln(0.5*nu_f2);
            
            loglik_r2 = -0.5 * log(h) - 0.5*(nu_t+1)*log(1 + eps2_vec./((nu_t-2)*h));
            loglik_r2 = loglik_r2 + gammaln(0.5*(nu_t+1)) - gammaln(0.5*nu_t) - 0.5 * log(pi*(nu_t-2));
            
            
            LLF= -sum(loglik_r2 + loglik_RV);
            LLF_t= -sum(loglik_r2);
            LLF_RV = -sum(loglik_RV);
            
            if  isnan(LLF) || isreal(LLF)==0|| isreal(loglik_RV)==0 || any(h<=0)
                LLF=1e14;
            end
            
            
        end
        
        
        function [LLF,h,eps_vec] = LogLik_GAS_t_VT_Uni(T,params,ret_vec,nr_lags,X_vec);
            % nr_lags
            
            c      = params(1);
            if nr_lags>0
                lag_matrix = lagmatrix(ret_vec,1:nr_lags);
                lag_matrix(1:nr_lags,:) = [];
                phi_vec = params(2:2+nr_lags-1);
                eps_vec = ret_vec(nr_lags+1:end) - c - lag_matrix * phi_vec;
                T = T-nr_lags;
            else
                eps_vec = ret_vec - c;
            end
            
            eps2_vec = eps_vec.^2;
            
            %omega = (1-beta)*fbar
            omega = params(2);
            alpha = params(3);
            beta  = params(4);
            nu_t  = params(5);
            
            
            if ~isempty(X_vec)
                gamma = params(6); 
            else
                gamma = 0;
                X_vec = zeros(T,1);
            end
                
            
            % score objects: score and GAS
            
            h0 = var(eps_vec);
            h = zeros(T,1);
            h(1) = h0;
            s = zeros(T,1);
            %fbar = h0;
            
            for j=1:T
                
                % update the score and cov
                w_j = (nu_t + 1)/ (nu_t - 2 +  eps2_vec(j,1)/h(j,1));
                d1 = w_j *  eps2_vec(j,1) - h(j,1);
                nabla = d1;
                s(j,1) = nabla;
                
                if j<T
                    h(j+1,1) = omega + beta * h(j,1) + alpha * s(j,1) + gamma * X_vec(j,1);
                end
                
            end
            
            loglik_r2 = -0.5 * log(h) - 0.5*(nu_t+1)*log(1 + eps2_vec./((nu_t-2)*h));
            loglik_r2 = loglik_r2 + gammaln(0.5*(nu_t+1)) - gammaln(0.5*nu_t) - 0.5 * log(pi*(nu_t-2));
            
            LLF= -sum(loglik_r2);
            
            if  isnan(LLF) || isreal(loglik_r2)==0
                LLF=1e6;
            end
            
            
        end
        
        function [LLF,loglik_RV,h,eps_vec] = LogLik_GAS_HAR_tF_VT_Uni(T,params,ret_vec,RV,nr_lags);
            % nr_lags
            l_2 = 12;
            l_3 = 60;
            
            c      = params(1);
            if nr_lags>0
                lag_matrix = lagmatrix(ret_vec,1:nr_lags);
                lag_matrix(1:nr_lags,:) = [];
                phi_vec = params(2:2+nr_lags-1);
                eps_vec = ret_vec(nr_lags+1:end) - c - lag_matrix * phi_vec;
                T = T-nr_lags;
                RV = RV(nr_lags+1:end);
            else
                eps_vec = ret_vec - c;
            end
            
            eps2_vec = eps_vec.^2;
            
            alpha = params(end-6);
            beta_1  = params(end-5);
            beta_2  = params(end-4);
            beta_3  = params(end-3);
            nu_t  = params(end-2);
            nu_f1 = params(end-1);
            nu_f2 = params(end);
            
            constant = nu_f1/ (nu_f2-2);
            % score objects: score and GAS
            
            h0 = mean(RV);
            h = zeros(T,1);
            h(1) = h0;
            s = zeros(T,1);
            fbar = h0;
            omega_zero = (1-beta_1-beta_2-beta_3)*fbar;
            
            for j=1:T
                
                % update the score and cov
                w_j = (nu_t + 1)/ (nu_t - 2 +  eps2_vec(j,1)/h(j,1));
                d1 = w_j *  eps2_vec(j,1) - h(j,1);
                d2 = nu_f1 *(((nu_f1 +nu_f2)/(nu_f2-2))*RV(j,1)/(1 + constant * RV(j,1)/h(j,1)) - h(j,1));
                
                nabla = (d1 + d2)/(nu_f1+1);
                s(j,1) = nabla;
                
                if j<T
                    if j < l_2
                        h(j+1,1) =  omega_zero + beta_1* h(j,1) + beta_2* mean(h(1:j)) + beta_3*  mean(h(1:j)) + alpha * s(j,1);
                    elseif j<l_3
                        h(j+1,1) =  omega_zero + beta_1* h(j,1) + beta_2* mean(h(j-l_2+1:j)) + beta_3* mean(h(1:j)) + alpha * s(j,1);
                    else
                        h(j+1,1) =  omega_zero + beta_1* h(j,1) + beta_2*  mean(h(j-l_2+1:j)) + beta_3*  mean(h(j-l_3+1:j)) + alpha * s(j,1);
                    end
                    
                end
                
            end
            
            loglik_RV = (0.5*nu_f1 - 1) * log(RV) -0.5*nu_f1 * log(h) -0.5*(nu_f1+nu_f2)*log(1 + (1./h)*constant.*RV) ...
                + 0.5*nu_f1 * log(constant);
            loglik_RV = loglik_RV + gammaln(0.5*(nu_f1+nu_f2)) - gammaln(0.5*nu_f1) - gammaln(0.5*nu_f2);
            
            loglik_r2 = -0.5 * log(h) - 0.5*(nu_t+1)*log(1 + eps2_vec./((nu_t-2)*h));
            loglik_r2 = loglik_r2 + gammaln(0.5*(nu_t+1)) - gammaln(0.5*nu_t) - 0.5 * log(pi*(nu_t-2));
            
            
            LLF= -sum(loglik_r2 + loglik_RV);
            
            if  isnan(LLF) || isreal(loglik_RV)==0
                LLF=1e6;
            end
            
            
        end
        
        function [LLF,h,eps_vec] = LogLik_Heavy_N(T,params,ret_vec,RV,nr_lags);
            
            c      = params(1);
            if nr_lags>0
                lag_matrix = lagmatrix(ret_vec,1:nr_lags);
                lag_matrix(1:nr_lags,:) = [];
                phi_vec = params(2:2+nr_lags-1);
                eps_vec = ret_vec(nr_lags+1:end) - c - lag_matrix * phi_vec;
                T = T-nr_lags;
                RV = RV(nr_lags+1:end);
            end
            
            eps2_vec = eps_vec.^2;
            %omega = (1-beta)*fbar
            omega = params(end-3);
            alpha = params(end-2);
            beta  = params(end-1);
            
            
            h = zeros(T,1);
            h(1) = mean(RV);
            
            for j = 1:T-1
                h(j+1,1) = omega + beta * h(j,1) + alpha * RV(j,1);
            end
            
            loglik_r2 = -0.5*log(2*pi)-0.5*(log(h) + eps2_vec./h);
            
            LLF= -sum(loglik_r2);
            
            if isnan(LLF) || isreal(LLF) == 0 || isinf(abs(LLF))==1
                LLF = 1e14;
            end
            
            
            
        end
        
        function [LLF,h,eps_vec] = Loglike_HEAVY_t(T,params,ret_vec,RV,nr_lags);
            
            c      = params(1);
            if nr_lags>0
                lag_matrix = lagmatrix(ret_vec,1:nr_lags);
                lag_matrix(1:nr_lags,:) = [];
                phi_vec = params(2:2+nr_lags-1);
                eps_vec = ret_vec(nr_lags+1:end) - c - lag_matrix * phi_vec;
                T = T-nr_lags;
                RV = RV(nr_lags+1:end);
            else
                eps_vec = ret_vec-c;
            end
            
            eps2_vec = eps_vec.^2;
            %omega = (1-beta)*fbar
            omega = params(end-3);
            alpha = params(end-2);
            beta  = params(end-1);
            nu    = params(end);
            
            h = zeros(T,1);
            h(1) = mean(RV);
            
            for j = 1:T
                
                if j<T
                    h(j+1,1) = omega + beta * h(j,1) + alpha * RV(j,1);
                end
                
            end
            
            LLF = T*gammaln(0.5*(nu+1)) - T*gammaln(nu/2) - T/2*log(pi*(nu-2));
            LLF = LLF - 0.5*sum(log(h)) - ((nu+1)/2)*sum(log(1 + eps2_vec./(h*(nu-2)) ));
            LLF = -LLF;
            
            if isnan(LLF) || isreal(LLF) == 0 || isinf(abs(LLF))==1
                LLF = 1e14;
            end
            
        end
        
        
        function [LLF,h,eps_vec] = Loglike_GARCH_t(T,params,ret_vec);
            
            c      = params(1);
            
            eps_vec = ret_vec-c;
            eps2_vec = eps_vec.^2;
            %omega = (1-beta)*fbar
            omega = params(end-3);
            alpha = params(end-2);
            beta  = params(end-1);
            nu    = params(end);
            
            h = zeros(T,1);
            h(1) = mean(eps2_vec);
            
            for j = 1:T
                
                if j<T
                    h(j+1,1) = omega + beta * h(j,1) + alpha * eps2_vec(j,1);
                end
                
            end
            
            LLF = T*gammaln(0.5*(nu+1)) - T*gammaln(nu/2) - T/2*log(pi*(nu-2));
            LLF = LLF - 0.5*sum(log(h)) - ((nu+1)/2)*sum(log(1 + eps2_vec./(h*(nu-2)) ));
            LLF = -LLF;
            
            if isnan(LLF) || isreal(LLF) == 0 || isinf(abs(LLF))==1
                LLF = 1e14;
            end
            
        end
        
        
        function [LLF,h_vec_TOT,eps_vec] = Loglike_GARCH_N(T,params,ret_vec,T_w);
            
            mu = params(1);
            alpha = params(end-2);            
            %gamma = params(end-2);
            beta  = params(end-1);            
            tau = params(end);
                      
            eps_vec = ret_vec-mu;
            %I_neg = eps_vec<=0;
            eps2_vec = eps_vec.^2;            
            
            h_vec_ST = zeros(T,1);
            h_vec_ST(1) = 1;
            
            for j = 1:T
                
                if j<T
                    h_vec_ST(j+1,1) =  (1-alpha-beta) + beta * h_vec_ST(j,1) + alpha* eps2_vec(j,1)/tau;
                end
                
            end
            
            h_vec_TOT = h_vec_ST * tau;
            loglik_r2 = -0.5*log(2*pi) -0.5 * log(h_vec_TOT) - 0.5*eps2_vec./(h_vec_TOT);            
            LLF= -sum(loglik_r2(T_w+1:end));
            
            if isnan(LLF) || ~isreal(LLF) || isinf(abs(LLF))
                LLF = 1e14;
            end
            
        end
        
        function [LLF,h_vec_TOT,eps_vec] = Loglike_GJR_GARCH_N(T,params,ret_vec,T_w);
            
            mu = params(1);
            alpha = params(end-3);            
            gamma = params(end-2);
            beta  = params(end-1);            
            tau = params(end);
                      
            eps_vec = ret_vec-mu;
            I_neg = eps_vec<=0;
            eps2_vec = eps_vec.^2;            
            
            h_vec_ST = zeros(T,1);
            h_vec_ST(1) = 1;
            
            for j = 1:T
                
                if j<T
                    h_vec_ST(j+1,1) =  (1-alpha-beta-gamma/2) + beta * h_vec_ST(j,1) + (alpha+gamma*I_neg(j,1)) * eps2_vec(j,1)/tau;
                end
                
            end
            h_vec_TOT = h_vec_ST * tau;
            loglik_r2 = -0.5*log(2*pi) -0.5 * log(h_vec_TOT) - 0.5*eps2_vec./(h_vec_TOT);                             
            LLF= -sum(loglik_r2(T_w+1:end));
            
            if isnan(LLF) || ~isreal(LLF) || isinf(abs(LLF))
                LLF = 1e14;
            end
            
        end
        
        function [LLF,h_tot,eps_vec] = Loglike_HEAVY_N(T,params,ret_vec,RV_vec,T_w);
            
                        
            mu     = params(1);          
            alpha = params(end-2);
            beta  = params(end-1);            
            tau = params(end);
            
            
            eps_vec = ret_vec-mu;
            eps2_vec = eps_vec.^2;                        
            
            
            h_vec = zeros(T,1);
            h_vec(1) = 1;
            
            for j = 1:T
                
                if j<T
                    h_vec(j+1,1) =  (1-alpha-beta) + beta * h_vec(j,1) + alpha * RV_vec(j,1)/tau;
                end
                
            end

            h_tot =h_vec * tau;
            loglik_r2 = -0.5*log(2*pi) -0.5 * log(h_tot) - 0.5*eps2_vec./(h_tot);            
            LLF= -sum(loglik_r2(T_w+1:end));
            
            if isnan(LLF) || ~isreal(LLF) || isinf(abs(LLF))
                LLF = 1e14;
            end
            
        end
        
        
         function [LLF,h_tot,eps_vec] = Loglike_HEAVY_GJR_N(T,params,ret_vec,RV_vec,T_w);
            
                        
            mu     = params(1);            
            alpha = params(end-3);
            gamma = params(end-2);
            beta  = params(end-1);            
            tau = params(end);
            
            
            eps_vec = ret_vec-mu;
            eps2_vec = eps_vec.^2;                        
            I_vec = eps_vec<=0;
            
            h_vec = zeros(T,1);
            h_vec(1) = 1;
            
            for j = 1:T
                
                if j<T
                    h_vec(j+1,1) =  (1-alpha-beta) + beta * h_vec(j,1) + (alpha + gamma*I_vec(j)) * RV_vec(j,1)/tau;
                end
                
            end

            h_tot =h_vec * tau;
            loglik_r2 = -0.5*log(2*pi) -0.5 * log(h_tot) - 0.5*eps2_vec./(h_tot);            
            LLF= -sum(loglik_r2(T_w+1:end));
            
            if isnan(LLF) || ~isreal(LLF) || isinf(abs(LLF))
                LLF = 1e14;
            end
            
        end
        
        
        function [LLF,h_vec,tau_vec,eps_vec,v_vec] = Loglike_GJR_GARCH_MF_N(T,params,ret_vec,T_w);
            
                        
            
            mu    = params(end-6);
            alpha = params(end-5);
            gamma = params(end-4);
            beta  = params(end-3);            
            lambda_0= params(end-2);           
            lambda_1= params(end-1);
            lambda_2 = params(end);
            
            
            eps_vec = ret_vec-mu;
            I_neg = eps_vec<0;
            
            eps2_vec = eps_vec.^2;
            
            omega = lambda_0;
            
            
            h_vec = zeros(T,1);
            h_vec(1) = 1;
            tau_vec = zeros(T,1);
            v_vec = zeros(T,1);
            tau_vec(1:T_w) = var(eps_vec);
            
            for j = 1:T
                
                if j<T                                                      
                        v_vec(j,1) = eps2_vec(j,1)/h_vec(j,1);                    
                        %h_vec(j+1,1) =  (1-alpha-beta) + beta * h_vec(j,1) + alpha * eps2_vec(j,1)/tau_vec(j,1);
                        h_vec(j+1,1) =  (1-alpha-beta-0.5*gamma) + beta * h_vec(j,1) + (alpha+gamma*I_neg(j))* eps2_vec(j,1)/tau_vec(j,1);
                        
                        if j>=T_w
                            tau_vec(j+1,1) = omega + lambda_1*mean(v_vec(j-T_w+1:j)) + lambda_2 * tau_vec(j,1);
                        end
                        
                end
                
            end
            
            loglik_r2 = -0.5*log(2*pi) -0.5 * log(h_vec) -0.5*log(tau_vec) - 0.5*eps2_vec./(h_vec.*tau_vec);              
            LLF= -sum(loglik_r2(T_w+1:end));
            %LLF= -sum(loglik_r2((2*252+1):end));
            
            
            
            
            
            if isnan(LLF) || ~isreal(LLF) || isinf(abs(LLF))
                LLF = 1e14;
            end
            
        end
        
        
        function [LLF,h_vec,tau_vec,eps_vec,v_vec] = Loglike_GJR_GARCH_MF_N_2(T,params,ret_vec,T_w);
            
                        
            
            mu    = params(end-6);
            alpha = params(end-5);
            gamma = params(end-4);
            beta  = params(end-3);            
            lambda_0= params(end-2);           
            lambda_1= params(end-1);
            lambda_2 = params(end);
            
            
            eps_vec = ret_vec-mu;
            I_neg = eps_vec<0;
            
            eps2_vec = eps_vec.^2;
            
            omega = lambda_0*(1-lambda_1-lambda_2);
            
            
            h_vec = zeros(T,1);
            h_vec(1) = 1;
            tau_vec = zeros(T,1);
            v_vec = zeros(T,1);
            tau_vec(1:T_w) = var(eps_vec);
            
            for j = 1:T
                
                if j<T                                                      
                        v_vec(j,1) = eps2_vec(j,1)/h_vec(j,1);                    
                        %h_vec(j+1,1) =  (1-alpha-beta) + beta * h_vec(j,1) + alpha * eps2_vec(j,1)/tau_vec(j,1);
                        h_vec(j+1,1) =  (1-alpha-beta-0.5*gamma) + beta * h_vec(j,1) + (alpha+gamma*I_neg(j))* eps2_vec(j,1)/tau_vec(j,1);
                        
                        if j>=T_w
                            tau_vec(j+1,1) = omega + lambda_1*mean(v_vec(j-T_w+1:j)) + lambda_2 * tau_vec(j,1);
                        end
                        
                end
                
            end
            
            loglik_r2 = -0.5*log(2*pi) -0.5 * log(h_vec) -0.5*log(tau_vec) - 0.5*eps2_vec./(h_vec.*tau_vec);              
            LLF= -sum(loglik_r2(T_w+1:end));
            %LLF= -sum(loglik_r2((2*252+1):end));
            
            
            
            
            
            if isnan(LLF) || ~isreal(LLF) || isinf(abs(LLF))
                LLF = 1e14;
            end
            
        end
        
        
        
        function [LLF,h_vec,tau_vec,eps_vec] = Loglike_GJR_GARCH_MF_N_VT(T,params,ret_vec,T_w);
            
                        
            
            mu    = params(end-5);
            alpha = params(end-4);
            gamma = params(end-3);
            beta  = params(end-2);                        
            lambda_1= params(end-1);
            lambda_2 = params(end);
            
            
            eps_vec = ret_vec-mu;
            I_neg = eps_vec<0;
            
            eps2_vec = eps_vec.^2;
            
            omega = var(eps_vec) * (1-lambda_1-lambda_2);
            
            
            h_vec = zeros(T,1);
            h_vec(1) = 1;
            tau_vec = zeros(T,1);
            v_vec = zeros(T,1);
            tau_vec(1:T_w) = var(eps_vec);
            
            for j = 1:T
                
                if j<T                                                      
                        v_vec(j,1) = eps2_vec(j,1)/h_vec(j,1);                    
                        %h_vec(j+1,1) =  (1-alpha-beta) + beta * h_vec(j,1) + alpha * eps2_vec(j,1)/tau_vec(j,1);
                        h_vec(j+1,1) =  (1-alpha-beta-0.5*gamma) + beta * h_vec(j,1) + (alpha+gamma*I_neg(j))* eps2_vec(j,1)/tau_vec(j,1);
                        
                        if j>=T_w
                            tau_vec(j+1,1) = omega + lambda_1*mean(v_vec(j-T_w+1:j)) + lambda_2 * tau_vec(j,1);
                        end
                        
                end
                
            end
            
            loglik_r2 = -0.5*log(2*pi) -0.5 * log(h_vec) -0.5*log(tau_vec) - 0.5*eps2_vec./(h_vec.*tau_vec);              
            LLF= -sum(loglik_r2(T_w+1:end));
            
            if isnan(LLF) || ~isreal(LLF) || isinf(abs(LLF))
                LLF = 1e14;
            end
            
        end
        
        
        function [LLF,h_vec,tau_vec,eps_vec] = Loglike_GARCH_RV_MF_N(T,params,ret_vec,RV_vec,T_w);
            
            %Aleen met otc data           
            
            mu    = params(1);                        
            alpha = params(end-4);           
            beta  = params(end-3);            
            lambda_0= params(end-2);
            lambda_1= params(end-1);
            lambda_2 = params(end);
            
            
            eps_vec = ret_vec-mu;
            eps2_vec = eps_vec.^2;            
            
            I_vec = eps_vec<=0;
            
            h_vec = zeros(T,1);
            h_vec(1) = 1;
            tau_vec = zeros(T,1);
            v_vec = zeros(T,1);
            tau_vec(1:T_w) = var(eps_vec);
            
            for j = 1:T
                
                if j<T                                                      
                        v_vec(j,1) = RV_vec(j,1)/h_vec(j,1);                    
                        h_vec(j+1,1) =  (1-alpha-beta) + beta * h_vec(j,1) + alpha* RV_vec(j,1)/tau_vec(j,1);
                        
                        if j>=T_w
                            tau_vec(j+1,1) = lambda_0 + lambda_1*mean(v_vec(j-T_w+1:j)) + lambda_2 * tau_vec(j,1);
                        end
                        
                end
                
            end
            
            h_tot= h_vec.*tau_vec;
            loglik_r2 = -0.5*log(2*pi) -0.5 * log(h_tot) - 0.5*eps2_vec./h_tot;              
            LLF= -sum(loglik_r2(T_w+1:end));
            
            if isnan(LLF) || ~isreal(LLF) || isinf(abs(LLF))
                LLF = 1e14;
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
        
        
        function [LLF,mu_vec,s_vec] = LogLik_GAS_X_N_TV_mu(T,params,y,X);
            % nr_lags
            
            omega   = params(1);
            alpha   = params(2);
            beta    = params(3);
            sigma2  = params(4);
            
            if size(X,1)>0
                gamma = params(end);
            else
                X = zeros(T,1);
                gamma = 0;
            end
                                                         
            
            mu_vec = zeros(T,1);
            eps2_vec = zeros(T,1);
            mu_vec(1) = mean(y);
            s_vec = zeros(T,1);
            
            for j=1:T
                
                % update the score and cov
                eps_t = y(j,1)-mu_vec(j,1);              
                eps2_vec(j,1) = eps_t^2;
                s_mu = eps_t;
                s_vec(j,1) = s_mu;
                                
                if j<T
                    mu_vec(j+1,1) = omega + beta * mu_vec(j,1) + alpha * s_mu + gamma * X(j,1);
                end
                
            end
                        
            loglik_r2 = -0.5*log(2*pi) -0.5 * log(sigma2) - 0.5*eps2_vec./sigma2;            
            LLF= -sum(loglik_r2);
            
            if  isnan(LLF) || isreal(LLF)==0
                LLF=1e14;
            end
                                  
            
       end
       
       
        
        
        
        
        
    end
end