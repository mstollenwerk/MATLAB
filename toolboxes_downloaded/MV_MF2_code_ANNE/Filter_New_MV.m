classdef Filter_New_MV % in Filter.m
    methods (Static = true)
        
        
        %%
        function [s, f, V] = Store(k, T, V0);
            
            %s  = zeros(k*(k+1)/2, 1);
            %f  = zeros(k*(k+1)/2, 1); f(:,1)= Admin.LowerTr2Vech(k, transpose(chol(V0)));
            %V  = zeros(k, k, T); V(:,:,1) = V0;
            
            s = zeros(k,k,T);
            f = zeros(k,k,T); f(:,:,1) = V0;
            V  = zeros(k, k, T); V(:,:,1) = V0;
            
        end
        
        
        %%
        function [ V, s, f,Flag] = GAS( k, T, params, r2, RC )
            %Simulate the cont. process with n intraday data for T days with V covariance
            
            % assign params
            omega = params(1);
            alpha = params(2);
            beta  = params(3);
            nu_t  = params(4);
            nu_f1 = params(5);
            nu_f2 = params(6);
            
            % score objects: score and GAS
            V0 = mean(RC, 3);
            [s, f, V] = Filter_New.Store(k, T, V0);
            constant = nu_f1/ (nu_f2-k-1);
            
            for j=1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                % update the score and cov
                V_j      = V(:,:,j);
                RC_j     = RC(:,:,j);
                invV     = inv( V_j );
                w_j      = (nu_t + k)/ (nu_t - 2 + trace(invV*r2(:,:,j)));
                d1       =  w_j * r2(:,:,j) - V_j;
                d2       = nu_f1 * 1*( (nu_f1 +nu_f2)/(nu_f2-k-1)*RC_j*inv(eye(k)+ constant*invV*RC_j) - V_j);   % geschaalde score uit F verdeling
                
                nabla = (d1 + d2)/(nu_f1 + 1);
                s(:,:,j) = nabla;
                
                if j<T
                    f(:,:,j+1) = omega + beta * f(:,:,j) + alpha * s(:,:,j);
                    V(:,:,j+1) = f(:,:,j+1);
                end
                
                
            end
        end
        
        
        %% GAS filter with CT = Covariance Targeting
        function [ V, s,Flag] = GAS_tF_CT( k, T, params, return_mat, RC )
            %Simulate the cont. process with n intraday data for T days with V covariance
            
            
            % assign params
            alpha   = params(1);
            beta    = params(2);
            nu_t    = params(3);
            nu_f1   = params(4);
            nu_f2   = params(5);
            Flag = 0;
            
            % score objects: score and GAS
            V0  = mean(RC, 3);
            fbar= V0;
            [s, f, V] = Filter_New.Store(k, T, V0);
            constant = nu_f1/ (nu_f2-k-1);
            for j=1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                [~,nr_pd] = chol(V(:,:,j));
                if nr_pd>0
                    Flag = 1;
                    break
                end
                
                % update the score and cov
                V_j   = V(:,:,j);
                RC_j  = RC(:,:,j);
                invV  = inv( V_j );
                r2_j = return_mat(j,:)'*return_mat(j,:);
                
                w_j   = (nu_t + k)/ (nu_t - 2 + trace(invV*r2_j));
                d1    =  1*w_j * r2_j - 1*V_j;
                d2    = 1*nu_f1*((nu_f1 +nu_f2)/(nu_f2-k-1)*RC_j*inv(eye(k)+ constant*invV*RC_j) - V_j) ;
                
                %w_vec(j,1) = w_j;
                %d1_mat(:,:,j) = (alpha * d1) / (nu_f1 + 1);
                %d2_mat(:,:,j) = (alpha * d2)/ (nu_f1+1);
                nabla = (d1 + d2)/(1+nu_f1);
                %nabla = (d2)/(1+nu_f1);
                s(:,:,j) = nabla;
                
                if j<T
                    f(:,:,j+1) = (1-beta)*fbar + beta * f(:,:,j) + alpha * s(:,:,j);
                    V(:,:,j+1) = f(:,:,j+1);
                end
                
            end
            
            
        end
        
        
        %% GAS filter with CT = Covariance Targeting
        function [ V, s,Flag] = GAS_t_CT( k, T, params, return_mat)
            %Simulate the cont. process with n intraday data for T days with V covariance
            
            
            % assign params
            alpha   = params(1);
            beta    = params(2);
            nu_t    = params(3);
            Flag = 0;
            
            % score objects: score and GAS
            V0  = cov(return_mat);
            fbar= V0;
            [s, f, V] = Filter_New.Store(k, T, V0);
            
            for j=1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                [~,nr_pd] = chol(V(:,:,j));
                if nr_pd>0
                    Flag = 1;
                    break
                end
                
                % update the score and cov
                V_j   = V(:,:,j);
                invV  = inv( V_j );
                r2_j = return_mat(j,:)'*return_mat(j,:);
                
                w_j   = (nu_t + k)/ (nu_t - 2 + trace(invV*r2_j));
                nabla    =  w_j * r2_j - V_j;
                s(:,:,j) = nabla;
                
                if j<T
                    f(:,:,j+1) = (1-beta)*fbar + beta * f(:,:,j) + alpha * s(:,:,j);
                    V(:,:,j+1) = f(:,:,j+1);
                end
                
            end
            
            
        end
        
        
        function [ V, s,Flag,nu_vec,f_nu_vec_score] = GAS_t_CT_tv_nu( k, T, params, return_mat)
            %Simulate the cont. process with n intraday data for T days with V covariance
            
            
            % assign params
            alpha   = params(1);
            beta    = params(2);
            omega_nu    = params(3);
            
            alpha_nu    = params(4);
            beta_nu    = params(5);
            
            %omega_nu = nu_bar*(1-beta_nu);
            
            Flag = 0;
            
            f_nu_vec = zeros(T,1);
            f_nu_vec_score = zeros(T,1);
            f_nu_vec(1) = params(6);
            
            % score objects: score and GAS
            V0  = cov(return_mat);
            fbar= V0;
            [s, f, V] = Filter_New.Store(k, T, V0);
            
            for j=1:T
                
                nu_t = 2 + exp(f_nu_vec(j));
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15 || isnan(nu_t) || isinf(nu_t)|| ~isreal(nu_t)
                    Flag = 1;
                    break
                end
                
                [~,nr_pd] = chol(V(:,:,j));
                if nr_pd>0
                    Flag = 1;
                    break
                end
                
                % update the score and cov
                V_j   = V(:,:,j);
                invV  = inv( V_j );
                r2_j = return_mat(j,:)'*return_mat(j,:);
                
                w_j   = (nu_t + k)/ (nu_t - 2 + trace(invV*r2_j));
                nabla    =  w_j * r2_j - V_j;
                s(:,:,j) = nabla;
                
                                              
                eps_V_inv_eps_t = (return_mat(j,:)*invV) * return_mat(j,:)';
                nu_t_min_2 = nu_t -2;
                                        
                    
                nu_score_t = 0.5*psi(0.5*(nu_t+k)) - 0.5*psi(0.5 *nu_t) -0.5 * log(1 +  eps_V_inv_eps_t/nu_t_min_2)...
                        +0.5 * ((nu_t+1)/nu_t_min_2) * 1/(nu_t_min_2 + eps_V_inv_eps_t) - 0.5 * k/(pi*nu_t_min_2);
                
                    
                nu_score_f_t = nu_score_t * exp(f_nu_vec(j));                    
                f_nu_vec_score(j,1) = nu_score_f_t;
                
                
                
                if j<T
                    f(:,:,j+1) = (1-beta)*fbar + beta * f(:,:,j) + alpha * s(:,:,j);
                    V(:,:,j+1) = f(:,:,j+1);
                    
                   f_nu_vec(j+1,1) = omega_nu + alpha_nu * nu_score_f_t + beta_nu * f_nu_vec(j,1);                                                            
                    
                end
                
            end
            
            nu_vec = 2 + exp(f_nu_vec);
            
        end
        
        
         
        
        
        function [V_mat_vech,Flag] = GAS_t_CT_fcst( k, T, params, return_mat,T_IS)
            %Simulate the cont. process with n intraday data for T days with V covariance
            
            
            % assign params
            alpha   = params(1);
            beta    = params(2);
            nu_t    = params(3);
            
            
            % score objects: score and GAS
            Vbar  = cov(return_mat(1:T_IS,:));
            
            Flag = 0;
            nr_vech = k*(k+1)/2;
            V_old =  Vbar;
            V_mat_vech     = zeros(T-T_IS,nr_vech);
            
            
            for j=1:T
                
                check_ind = rcond(V_old);
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                [~,nr_pd] = chol(V_old);
                if nr_pd>0
                    Flag = 1;
                    break
                end
                
                
                if j>T_IS
                    V_mat_vech(j-T_IS,:) = Admin.Sym2Vech(k,V_old);
                end
                
                
                % update the score and cov
                
                invV  = inv(V_old);
                r2_j = return_mat(j,:)'*return_mat(j,:);
                w_j   = (nu_t + k)/ (nu_t - 2 + trace(invV*r2_j));
                nabla    =  w_j * r2_j - V_old;
                
                if j<T
                    Vnew = (1-beta)*Vbar + beta * V_old + alpha * nabla;
                    V_old = Vnew;
                end
                
            end
            
            
        end
        
        
        function [ V, s,Flag] = GAS_CT_c( k, T, params, r2, RC )
            %Simulate the cont. process with n intraday data for T days with V covariance
            
            
            % assign params
            alpha   = params(1);
            beta    = params(2);
            nu_t    = params(3);
            nu_f1   = params(4);
            nu_f2   = params(5);
            c       = params(6);
            
            Flag = 0;
            
            % score objects: score and GAS
            V0  = mean(RC, 3);
            fbar= V0;
            [s, f, V] = Filter_New.Store(k, T, V0);
            constant = nu_f1/ (nu_f2-k-1);
            
            for j=1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                [~,nr_pd] = chol(V(:,:,j));
                if nr_pd>0
                    Flag = 1;
                    break
                end
                
                % update the score and cov
                V_j   = V(:,:,j);
                RC_j  = RC(:,:,j);
                invV  = inv( V_j );
                w_j   = (nu_t + k)/ (nu_t - 2 + trace(invV*r2(:,:,j)));
                d1    =  1*w_j * r2(:,:,j) - 1 * V_j;
                d2    = nu_f1*1*((nu_f1 +nu_f2)/(nu_f2-k-1)*RC_j*inv(eye(k)+ constant*invV*RC_j) - V_j) ;
                
                %w_vec(j,1) = w_j;
                %d1_mat(:,:,j) = (alpha * d1) / (nu_f1 + 1);
                %d2_mat(:,:,j) = (alpha * d2)/ (nu_f1+1);
                nabla = (d1 + d2)/(1+nu_f1);
                %nabla = (d2)/(1+nu_f1);
                s(:,:,j) = nabla;
                
                if j<T
                    f(:,:,j+1) = c*(1-beta)*fbar + beta * f(:,:,j) + alpha * s(:,:,j);
                    V(:,:,j+1) = f(:,:,j+1);
                end
                
            end
            
            
        end
        
        function [V,s,Flag]  = GAS_HAR_CT(k, T,params,r2,RC);
            alfa   = params(1);
            beta_1 = params(2);
            beta_2 = params(3);
            beta_3 = params(4);
            nu_t    = params(5);
            nu_f1   = params(6);
            nu_f2   = params(7);
            
            %l_2 = 5;
            %l_3 = 22;
            
            l_2 = 12;
            l_3 = 60;
            
            
            Flag = 0;
            V0  = mean(RC, 3);
            B0   = (1-beta_1-beta_2-beta_3)*V0;
            [s,~,V] = Filter_New.Store(k, T, V0);
            constant = nu_f1/ (nu_f2-k-1);
            for j = 1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                [~,nr_pd] = chol(V(:,:,j));
                if nr_pd>0
                    Flag = 1;
                    break
                end
                
                % update the score and cov
                V_j   = V(:,:,j);
                RC_j  = RC(:,:,j);
                invV  = inv( V_j );
                w_j   = (nu_t + k)/ (nu_t - 2 + trace(invV*r2(:,:,j)));
                d1    =  1*w_j * r2(:,:,j) - 1 * V_j;
                d2    = nu_f1*1*((nu_f1 +nu_f2)/(nu_f2-k-1)*RC_j*inv(eye(k)+ constant*invV*RC_j) - V_j) ;
                nabla = (d1 + d2)/(1+nu_f1);
                s(:,:,j) = nabla;
                
                if j<T
                    if j<l_2
                        V(:,:,j+1) = B0 + alfa*nabla + beta_1 * V(:,:,j) + beta_2 * mean(V(:,:,1:j),3) ...
                            + beta_3 * mean(V(:,:,1:j),3);
                    elseif j<l_3
                        
                        V(:,:,j+1) = B0 + alfa*nabla + beta_1 * V(:,:,j) + beta_2 * mean(V(:,:,j-l_2+1:j),3) ...
                            + beta_3 * mean(V(:,:,1:j),3);
                    else
                        V(:,:,j+1) = B0 + alfa*nabla + beta_1 * V(:,:,j) + beta_2 * mean(V(:,:,j-l_2+1:j),3)...
                            + beta_3 * mean(V(:,:,j-l_3+1:j),3);
                    end
                end
            end
        end
        
        
        
        
        
        function [ V, s,Flag] = GAS_Gen_CT( k, T, params, r2, RC )
            %Simulate the cont. process with n intraday data for T days with V covariance
            
            
            % assign params
            alpha_vec   = params(1:k);
            A_mat = diag(alpha_vec);
            %A_mat = alpha_vec'*alpha_vec;
            
            beta    = params(k+1);
            nu_t    = params(end-2);
            nu_f1   = params(end-1);
            nu_f2   = params(end);
            Flag = 0;
            
            % score objects: score and GAS
            V0  = mean(RC, 3);
            fbar= V0;
            [s, f, V] = Filter_New.Store(k, T, V0);
            constant = nu_f1/ (nu_f2-k-1);
            
            initial_A_cond = 2*beta*ones(k,1) - diag(A_mat^2);
            
            if any(initial_A_cond<0)
                Flag = 1;
                
            else
                
                for j=1:T
                    
                    check_ind = rcond(V(:,:,j));
                    if check_ind < 1e-15
                        Flag = 1;
                        break
                    end
                    
                    [~,temp_chol_ind] = chol(V(:,:,j));
                    if temp_chol_ind>0
                        Flag = 1;
                        break
                    end
                    
                    % update the score and cov
                    V_j   = V(:,:,j);
                    RC_j  = RC(:,:,j);
                    invV  = inv( V_j );
                    w_j   = (nu_t + k)/ (nu_t - 2 + trace(invV*r2(:,:,j)));
                    %d1    = 1 * w_j * r2(:,:,j) - 1 * V_j;
                    d1    = 1*w_j * r2(:,:,j) - 1*V_j;
                    d2    = 1*nu_f1*((nu_f1 +nu_f2)/(nu_f2-k-1)*RC_j*inv(eye(k)+ constant*invV*RC_j) - V_j) ;
                    
                    %w_vec(j,1) = w_j;
                    %d1_mat(:,:,j) = (alpha * d1) / (nu_f1 + 1);
                    %d2_mat(:,:,j) = (alpha * d2)/ (nu_f1+1);
                    nabla = (d1 + d2)/(1+nu_f1);
                    %nabla = (d2)/(1+nu_f1);
                    s(:,:,j) = nabla;
                    
                    if j<T
                        f(:,:,j+1) = (1-beta)*fbar + beta * f(:,:,j) + A_mat * s(:,:,j) * A_mat;
                        V(:,:,j+1) = f(:,:,j+1);
                    end
                    
                end
                
                
            end
        end
        
        
        function [ V, s,Flag] = GAS_CT_W_dyn( k, T, params, r2, RC )
            %Simulate the cont. process with n intraday data for T days with V covariance
            
            
            % if size(V, 1) ~= k
            %     display('check dims');
            % end
            
            
            % assign params
            alpha   = params(1);
            beta    = params(2);
            nu_f1   = params(3);
            Flag = 0;
            
            % score objects: score and GAS
            V0  = mean(RC, 3);
            fbar= V0;
            [s, f, V] = Filter_New.Store(k, T, V0);
            for j=1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                % update the score and cov
                V_j   = V(:,:,j);
                RC_j  = RC(:,:,j);
                d1    = 0.5 * (r2(:,:,j) - V_j);
                d2    = 0.5*nu_f1*(RC_j - V_j);
                
                nabla = (d1 + d2)/(1+nu_f1);
                s(:,:,j) = nabla;
                
                if j<T
                    f(:,:,j+1) = (1-beta)*fbar + beta * f(:,:,j) + alpha * s(:,:,j);
                    V(:,:,j+1) = f(:,:,j+1);
                end
                
            end
            
            
        end
        
        
        
        
        function [ V_mat_vech,Flag,V_mat_vech_IS,V_mat_vech_TOT] = GAS_tF_CT_fcst( k, T, params, return_mat, RC,T_IS )
            
            % general idea forecasting:
            % run de filter vanaf obs 1, maar wel met omega based on T_IS
            % gooi daarna de eerste T_IS observaties van V weg.
            
            % assign params
            alpha   = params(1);
            beta    = params(2);
            nu_t    = params(3);
            nu_f1   = params(4);
            nu_f2   = params(5);
            Flag = 0;
            
            % score objects: score and GAS
            Vbar  = mean(RC(:,:,1:T_IS), 3);
            
            nr_vech = k*(k+1)/2;
            V_old =  Vbar;
            V_mat_vech     = zeros(T,nr_vech);
            constant = nu_f1/ (nu_f2-k-1);
            
            for j=1:T
                
                check_ind = rcond(V_old);
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                
               V_mat_vech(j,:) = Admin.Sym2Vech(k,V_old);
              
                
                % update the score and cov
                
                RC_j  = RC(:,:,j);
                invV  = inv(V_old);
                r2_j  = return_mat(j,:)'*return_mat(j,:);
                
                
                w_j   = (nu_t + k)/ (nu_t - 2 + trace(invV*r2_j));
                d1 = 1 * w_j * r2_j - 1 * V_old;
                d2    = nu_f1*1*((nu_f1 +nu_f2)/(nu_f2-k-1)*RC_j*inv(eye(k)+ constant*invV*RC_j) - V_old) ;
                
                nabla_t = (d1 + d2)/(1+nu_f1);
                
                if j<T
                    V_new = (1-beta)*Vbar + beta * V_old + alpha * nabla_t;
                    V_old = V_new;
                end
                
            end
            
            V_mat_vech_IS  = V_mat_vech(1:T_IS,:);
            V_mat_vech_TOT =   V_mat_vech;
            V_mat_vech(1:T_IS,:) = [];
            
               
            
            
            
        end
        
        
        function [V,s,Flag]  = GAS_tF_HAR_CT_fcst(k, T,params,r2,RC,T_IS);
            alfa   = params(1);
            beta_1 = params(2);
            beta_2 = params(3);
            beta_3 = params(4);
            nu_t    = params(5);
            nu_f1   = params(6);
            nu_f2   = params(7);
            
            l_2 = 12;
            l_3 = 60;
            
            Flag = 0;
            V0  = mean(RC(:,:,1:T_IS), 3);
            B0   = (1-beta_1-beta_2-beta_3)*V0;
            [s,~,V] = Filter_New.Store(k, T, V0);
            constant = nu_f1/ (nu_f2-k-1);
            for j = 1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                [~,nr_pd] = chol(V(:,:,j));
                if nr_pd>0
                    Flag = 1;
                    break
                end
                
                % update the score and cov
                V_j   = V(:,:,j);
                RC_j  = RC(:,:,j);
                invV  = inv( V_j );
                w_j   = (nu_t + k)/ (nu_t - 2 + trace(invV*r2(:,:,j)));
                d1    =  1*w_j * r2(:,:,j) - 1 * V_j;
                d2    = nu_f1*1*((nu_f1 +nu_f2)/(nu_f2-k-1)*RC_j*inv(eye(k)+ constant*invV*RC_j) - V_j) ;
                nabla = (d1 + d2)/(1+nu_f1);
                s(:,:,j) = nabla;
                
                if j<T
                    if j<l_2
                        V(:,:,j+1) = B0 + alfa*nabla + beta_1 * V(:,:,j) + beta_2 * mean(V(:,:,1:j),3) ...
                            + beta_3 * mean(V(:,:,1:j),3);
                        
                        [~,nr_pd] = chol(V(:,:,j+1));
                        if nr_pd>0 && alfa>1
                            alfa_new = 0.99;
                            V(:,:,j+1) = B0 + alfa_new*nabla + beta_1 * V(:,:,j) + beta_2 * mean(V(:,:,1:j),3) ...
                                + beta_3 * mean(V(:,:,1:j),3);
                        end
                        
                        
                    elseif j<l_3
                        
                        V(:,:,j+1) = B0 + alfa*nabla + beta_1 * V(:,:,j) + beta_2 * mean(V(:,:,j-l_2+1:j),3) ...
                            + beta_3 * mean(V(:,:,1:j),3);
                        
                        [~,nr_pd] = chol(V(:,:,j+1));
                        if nr_pd>0 && alfa>1
                            alfa_new = 0.99;
                            V(:,:,j+1) = B0 + alfa_new*nabla + beta_1 * V(:,:,j) + beta_2 * mean(V(:,:,j-l_2+1:j),3) ...
                                + beta_3 * mean(V(:,:,1:j),3);
                            
                        end
                    else
                        V(:,:,j+1) = B0 + alfa*nabla + beta_1 * V(:,:,j) + beta_2 * mean(V(:,:,j-l_2+1:j),3)...
                            + beta_3 * mean(V(:,:,j-l_3+1:j),3);
                        
                        [~,nr_pd] = chol(V(:,:,j+1));
                        if nr_pd>0 && alfa>1
                            alfa_new = 0.99;
                            V(:,:,j+1) = B0 + alfa_new*nabla + beta_1 * V(:,:,j) + beta_2 * mean(V(:,:,j-l_2+1:j),3)...
                                + beta_3 * mean(V(:,:,j-l_3+1:j),3);
                            
                        end
                        
                        
                    end
                end
            end
            
            V(:,:,1:T_IS) = [];
        end
        
        
        
        function [ V, s,Flag] = GAS_CT_c_fcst( k, T, params, r2, RC,T_IS )
            
            % assign params
            alpha   = params(1);
            beta    = params(2);
            nu_t    = params(3);
            nu_f1   = params(4);
            nu_f2   = params(5);
            c       = params(6);
            
            Flag = 0;
            
            % score objects: score and GAS
            V0   = mean(RC(:,:,1:T_IS), 3);
            fbar = V0;
            [s, f, V] = Filter_New.Store(k, T, V0);
            constant = nu_f1/ (nu_f2-k-1);
            for j=1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                % update the score and cov
                V_j   = V(:,:,j);
                RC_j  = RC(:,:,j);
                invV  = inv( V_j );
                w_j   = (nu_t + k)/ (nu_t - 2 + trace(invV*r2(:,:,j)));
                d1 = 1 * w_j * r2(:,:,j) - 1 * V_j;
                d2    = nu_f1*1*((nu_f1 +nu_f2)/(nu_f2-k-1)*RC_j*inv(eye(k)+ constant*invV*RC_j) - V_j) ;
                
                nabla = (d1 + d2)/(1+nu_f1);
                s(:,:,j) = nabla;
                
                if j<T
                    f(:,:,j+1) = c*(1-beta)*fbar + beta * f(:,:,j) + alpha * s(:,:,j);
                    V(:,:,j+1) = f(:,:,j+1);
                end
                
            end
            
            V(:,:,1:T_IS) = [];
            
        end
        
        
        
        
        function [ V, s] = GAS_CT_2( k, T, params, r2, RC )
            %Simulate the cont. process with n intraday data for T days with V covariance
            
            
            % if size(V, 1) ~= k
            %     display('check dims');
            % end
            
            
            % assign params
            alpha1   = params(1);
            alpha2   = params(2);
            beta    = params(3);
            nu_t    = params(4);
            nu_f1   = params(5);
            nu_f2   = params(6);
            % score objects: score and GAS
            V0  = mean(RC, 3);
            fbar= V0;
            [s, f, V] = Filter_New.Store(k, T, V0);
            constant = nu_f1/ (nu_f2-k-1);
            
            for j=1:T
                
                % update the score and cov
                V_j   = V(:,:,j);
                RC_j  = RC(:,:,j);
                invV  = inv( V_j );
                w_j   = (nu_t + k)/ (nu_t - 2 + trace(invV*r2(:,:,j)));
                %V_inv_half = Admin.sqroot_invsymmat(V_j);
                d1 = 1 * w_j * r2(:,:,j) - 1 * V_j;
                %d1 = 0.5 * w_j * V_inv_half * r2(:,:,j)*V_inv_half - 0.5 * eye(k);
                
                %d2    = 0.5*(1 + nu_f1/nu_f2)*RC_j*inv(eye(k)+(nu_f1/nu_f2)*invV*RC_j) - 0.5*V_j ;
                d2    = nu_f1*1*((nu_f1 +nu_f2)/(nu_f2-k-1)*RC_j*inv(eye(k)+ constant*invV*RC_j) - V_j) ;
                
                %w_vec(j,1) = w_j;
                %d1_mat(:,:,j) = d1/(1+nu_f1);
                %d2_mat(:,:,j) = d2/(1+nu_f1);
                nabla1 = (d1)/(1+nu_f1);
                nabla2 = (d2)/(1+nu_f1);
                s(:,:,j) = nabla1 + nabla2;
                
                if j<T
                    f(:,:,j+1) = (1-beta)*fbar + beta * f(:,:,j) + alpha1*nabla1 + alpha2*nabla2;
                    V(:,:,j+1) = f(:,:,j+1);
                end
                
            end
            
            
        end
        
        
        
        function [V] = EWMA(k, T, params,RC)
            
            V0  = mean(RC, 3);
            b_EWMA = params;
            V = zeros(k,k,T);
            V(:,:,1) = V0;
            for j = 1:T
                if j<T
                    V(:,:,j+1) = b_EWMA * V(:,:,j) + (1-b_EWMA)*RC(:,:,j);
                end
            end
            
            
            
        end
        
        function [V,Flag] = CAW(k,p,T, params,RC)
            
            V0  = mean(RC, 3);
            V = zeros(k,k,T);
            V(:,:,1) = V0;
            Flag = 0;
            if p ==1
                A = params(1);
                B = params(2);
            elseif p ==2
                A_1 = params(1);
                A_2 = params(2);
                B_1 = params(3);
                B_2 = params(4);
            end
            
            if p ==1
                C_tar = (1-A-B)*V0;
                
                [~,a] = chol(C_tar);
                
                if a>0 || (A + B)>1
                    Flag = 1;
                else
                    
                    for j = 1:T
                        if j<T
                            V(:,:,j+1) =  C_tar + B * V(:,:,j) + A*RC(:,:,j);
                        end
                    end
                end
                
            elseif p ==2
                V(:,:,2) = V0;
                
                C_tar = (1-A_1 - A_2 - B_1 - B_2)*V0;
                [~,a] = chol(C_tar);
                
                if a>0  || (A_1+ A_2+ B_1+ B_2)>1
                    Flag = 1;
                else
                    for j = 2:T
                        if j<T
                            V(:,:,j+1) =  C_tar + B_1 * V(:,:,j) + B_2 * V(:,:,j-1)  + A_1*RC(:,:,j)  + A_2*RC(:,:,j-1);
                        end
                    end
                end
            end
            
        end
        
        function [V,Flag] = CAW_IW(k,p,T,Omega,params,RC)
            
            V0  = Omega;
            V = zeros(k,k,T);
            V(:,:,1) = V0;
            Flag = 0;
            if p ==1
                A = params(1);
                B = params(2);
            elseif p ==2
                A_1 = params(1);
                A_2 = params(2);
                B_1 = params(3);
                B_2 = params(4);
            end
            
            if p ==1
                [~,ind_chol] = chol(Omega);
                if (A + B)>1 || ind_chol>0
                    Flag = 1;
                else
                    for j = 1:T
                        if j<T
                            V(:,:,j+1) =  Omega + B * V(:,:,j) + A*RC(:,:,j);
                        end
                    end
                end
            elseif p ==2
                V(:,:,2) = V0;
                if (A_1+ A_2+ B_1+ B_2)>1
                    Flag = 1;
                else
                    for j = 2:T
                        if j<T
                            V(:,:,j+1) =  Omega + B_1 * V(:,:,j) + B_2 * V(:,:,j-1)  + A_1*RC(:,:,j)  + A_2*RC(:,:,j-1);
                        end
                    end
                end
            end
            
        end
        
        
        function [V,Flag] = CAW_diag(k,p,T, params,RC)
            
            V0  = mean(RC, 3);
            V = zeros(k,k,T);
            V(:,:,1) = V0;
            Flag = 0;
            if p ==1
                A = diag(params(1:k));
                B = diag(params(k+1:end-1));
            elseif p ==2
                A_1 = diag(params(1:k));
                A_2 = diag(params(k+1:2*k));
                B_1 = diag(params(2*k+1:3*k));
                B_2 = diag(params(3*k+1:4*k));
            end
            
            if p ==1
                C_tar = V0 - A*V0*A - B*V0*B;
                
                [~,a] = chol(C_tar);
                
                if a>0 || any(diag(A*A + B*B)>1)
                    Flag = 1;
                else
                    
                    for j = 1:T
                        if j<T
                            V(:,:,j+1) =  C_tar + B * V(:,:,j)*B + A*RC(:,:,j)*A;
                        end
                    end
                end
                
            elseif p ==2
                V(:,:,2) = V0;
                
                C_tar = V0 - A_1*V0*A_1 - A_2*V0*A_2 - B_1*V0*B_1 - B_2*V0*B_2;
                [~,a] = chol(C_tar);
                
                if a>0  || any(diag(A_1*A_1 + A_2*A_2 + B_1*B_1 + B_2*B_2)>1)
                    Flag = 1;
                else
                    for j = 2:T
                        if j<T
                            V(:,:,j+1) =  C_tar + B_1 * V(:,:,j)*B_1 + B_2 * V(:,:,j-1)*B_2  + A_1*RC(:,:,j)*A_1  + A_2*RC(:,:,j-1)*A_2;
                        end
                    end
                end
            end
            
        end
        
        
        function [loglik_RC,V,Flag] = HEAVY_F_fcst(k,T, params,RC,T_IS);
            
            % assign params
            alpha   = params(1);
            beta    = params(2);
            nu_f1    = params(3);
            nu_f2    = params(4);
            Flag = 0;
            
            constant = nu_f1/ (nu_f2-k-1);
            
            % score objects: score and GAS
            V0  = mean(RC(:,:,1:T_IS), 3);
            V = zeros(k,k,T);
            loglik_RC = zeros(T,1);
            V(:,:,1) = V0;
            for j=1:T
                j_V          = V(:,:,j);
                j_RC         = RC(:,:,j);
                
                check_ind = rcond(j_V);
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                % update cov
                if j<T
                    V(:,:,j+1) = (1-beta-alpha)*V0 + beta * j_V + alpha * j_RC;
                end
                
                if j>T_IS
                    j_invV       = inv(j_V);
                    j_logdetRC   = log(det(j_RC));
                    
                    loglik_RC(j) = 0.5*nu_f1 * log(det(constant*j_invV)) + 0.5*(nu_f1-k-1)*j_logdetRC - 0.5*(nu_f1+nu_f2)*log(det(eye(k)+ constant*j_invV*j_RC));
                end
                
            end
            
            loglik_RC = loglik_RC + LogLik_MV.lnmultigamma(k, 0.5*(nu_f1+nu_f2)) - LogLik_MV.lnmultigamma(k, 0.5*nu_f1) - LogLik_MV.lnmultigamma(k, 0.5*nu_f2);
            
            V(:,:,1:T_IS) = [];
            loglik_RC(1:T_IS) = [];
            
        end
        
        
        
        function [loglik_vec,V,Flag] = CAW_fcst(k,T, params,RC,T_IS);
            
            % assign params
            alpha   = params(1);
            beta    = params(2);
            nu      = params(3);
            Flag = 0;
            
            % score objects: score and GAS
            V0  = mean(RC(:,:,1:T_IS), 3);
            V = zeros(k,k,T);
            loglik_vec = zeros(T,1);
            V(:,:,1) = V0;
            for j=1:T
                j_V          = V(:,:,j);
                j_RC         = RC(:,:,j);
                
                check_ind = rcond(j_V);
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                % update cov
                if j<T
                    V(:,:,j+1) = (1-beta-alpha)*V0 + beta * j_V + alpha * j_RC;
                end
                
                if j>T_IS
                    j_invV       = inv(j_V);
                    j_logdetV    = log(det(j_V));
                    j_logdetRC   = log(det(j_RC));
                    
                    loglik_vec(j,1) = - nu/2 * j_logdetV   +  (nu-k-1)/2 * j_logdetRC  - 0.5 * nu * trace(j_invV * j_RC);
                end
                
            end
            
            loglik_vec= loglik_vec -(k*nu)/2 * log(2)  + (k*nu)/2*log(nu) -  LogLik.lnmultigamma(k, nu/2 );
            
            V(:,:,1:T_IS) = [];
            loglik_vec(1:T_IS) = [];
            
        end
        
        
        function [loglik_RC,V,Flag] = CAW_Riesz_I_fcst(k,T, params,RC,T_IS);
            
            % assign params
            alpha   = params(1);
            beta    = params(2);
            nu_vec      = params(3:end);
            Flag = 0;
            
            % score objects: score and GAS
            V0  = mean(RC(:,:,1:T_IS), 3);
            V = zeros(k,k,T);
            loglik_RC = zeros(T,1);
            V(:,:,1) = V0;
            for j=1:T
                j_V          = V(:,:,j);
                j_RC         = RC(:,:,j);
                
                check_ind = rcond(j_V);
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                % update cov
                if j<T
                    V(:,:,j+1) = (1-beta-alpha)*V0 + beta * j_V + alpha * j_RC;
                end
                
                if j>T_IS
                    j_RC         = RC(:,:,j);
                    ln_PWD_j_RC  = LogLik.ln_PWD(j_RC,0.5*(nu_vec-k-1));
                    
                    L_tilde      = chol(j_V)';
                    L            = L_tilde * diag(1./sqrt(nu_vec));
                    Sigma_j      = L * L';
                    
                    ln_PWD_j_Sigma = LogLik.ln_PWD(Sigma_j,0.5*nu_vec);
                    j_inv_Sigma       = Sigma_j\eye(k);
                    
                    loglik_RC(j,1) = -ln_PWD_j_Sigma  + ln_PWD_j_RC - 0.5 * trace(j_inv_Sigma * j_RC);
                end
                
            end
            
            loglik_RC = loglik_RC - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
            
            V(:,:,1:T_IS) = [];
            loglik_RC(1:T_IS) = [];
            
        end
        
        
        function [loglik_RC,V,Flag] = CAW_Inv_Riesz_I_fcst(k,T, params,RC,T_IS);
            
            % assign params
            alpha   = params(1);
            beta    = params(2);
            nu_vec      = params(3:end);
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            Flag = 0;
            
            % score objects: score and GAS
            V0  = mean(RC(:,:,1:T_IS), 3);
            V = zeros(k,k,T);
            loglik_RC = zeros(T,1);
            V(:,:,1) = V0;
            for j=1:T
                j_V          = V(:,:,j);
                j_RC         = RC(:,:,j);
                
                check_ind = rcond(j_V);
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                % update cov
                if j<T
                    V(:,:,j+1) = (1-beta-alpha)*V0 + beta * j_V + alpha * j_RC;
                end
                
                if j>T_IS
                    j_RC       = RC(:,:,j);
                    j_inv_RC   = inv(j_RC);
                    ln_PWD_j_RC  = LogLik.ln_PWD(j_inv_RC,0.5*(nu_vec+k+1));
                    
                    j_Sigma_0        = j_V;
                    Sigma_inv_j_0    = j_Sigma_0\eye(k);
                    
                    L_0_trans = chol(Sigma_inv_j_0);
                    L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                    inv_Sigma_j = L_trans' * L_trans;
                    j_Sigma = inv(inv_Sigma_j);
                    
                    ln_PWD_Sigma_j = LogLik.ln_PWD(inv_Sigma_j,-0.5*nu_vec);
                    
                    %j_RC        = RC(:,:,j);
                    %j_inv_RC    = inv(j_RC);
                    
                    loglik_RC(j,1) = ln_PWD_Sigma_j  + ln_PWD_j_RC - 0.5 * trace(j_Sigma  * j_inv_RC);
                end
                
            end
            
            loglik_RC = loglik_RC - (nu_vec'*ones(k,1))/2 * log(2)  -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
            
            V(:,:,1:T_IS) = [];
            loglik_RC(1:T_IS) = [];
            
        end
        
        
        
        function [loglik_RC,V,Flag] = CAW_FRiesz_fcst(k,T, params,RC,T_IS)
            
            alpha   = params(1);
            beta    = params(2);
            nu1_vec   = params(3:3+k-1);
            nu2_vec = params(3+k:3+2*k-1);
            gamma_t = 0.5*(k+1) - [k:-1:1]';
            
            mu_vec = Admin.Compute_1st_moment_FRiesz(nu1_vec,nu2_vec);
            
            Flag = 0;
            
            % score objects: score and GAS
            V0  = mean(RC(:,:,1:T_IS), 3);
            V = zeros(k,k,T);
            loglik_RC = zeros(T,1);
            V(:,:,1) = V0;
            for j=1:T
                j_V          = V(:,:,j);
                j_RC         = RC(:,:,j);
                
                check_ind = rcond(j_V);
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                % update cov
                if j<T
                    V(:,:,j+1) = (1-beta-alpha)*V0 + beta * j_V + alpha * j_RC;
                end
                
                if j>T_IS
                    ln_PWD_j_RC  = LogLik.ln_PWD(j_RC,0.5*(nu1_vec-k-1));
                    
                    j_Sigma_0        = V(:,:,j);
                    L_trans_0        = chol(j_Sigma_0);
                    L_trans          = diag(1./sqrt(mu_vec))*L_trans_0;
                    j_Sigma          = L_trans' * L_trans;
                    
                    j_X = j_RC + j_Sigma;
                    
                    ln_PWD_X_j = LogLik.ln_PWD(j_X,-0.5*(nu2_vec+nu1_vec)+gamma_t);
                    ln_PWD_Sigma_j = LogLik.ln_PWD(j_Sigma,-0.5*nu2_vec+gamma_t);
                    
                    loglik_RC(j,1) = -ln_PWD_Sigma_j  + ln_PWD_j_RC + ln_PWD_X_j;
                end
                
            end
            
            loglik_RC = loglik_RC + LogLik.lnmultigamma_nu_vec(k, (nu2_vec+nu1_vec)/2 ) ...
                -  LogLik.lnmultigamma_nu_vec(k, nu2_vec/2 ) - LogLik.lnmultigamma_nu_vec(k, nu1_vec/2 );
            
            V(:,:,1:T_IS) = [];
            loglik_RC(1:T_IS) = [];
            
            
        end
        
        
        function [loglik_vec,V,Flag] = CAIW_fcst(k,T, params,RC,T_IS);
            
            % assign params
            alpha   = params(1);
            beta    = params(2);
            nu      = params(3);
            Flag = 0;
            
            % score objects: score and GAS
            V0  = mean(RC(:,:,1:T_IS), 3);
            V = zeros(k,k,T);
            loglik_vec = zeros(T,1);
            V(:,:,1) = V0;
            constant = (nu-k-1);
            
            for j=1:T
                j_V          = V(:,:,j);
                j_RC         = RC(:,:,j);
                
                check_ind = rcond(j_V);
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                % update cov
                if j<T
                    V(:,:,j+1) = (1-beta-alpha)*V0 + beta * j_V + alpha * j_RC;
                end
                
                if j>T_IS
                    j_invRC      = inv(j_RC);
                    j_logdetV    = log(det(j_V));
                    j_logdetRC   = log(det(j_RC));
                    
                    loglik_vec(j) = (nu/2) * j_logdetV -0.5*(nu+k+1)*j_logdetRC -constant/2 * trace(j_invRC*j_V);
                end
                
            end
            
            loglik_vec = loglik_vec + 0.5*nu*k*log(constant/2) - LogLik.lnmultigamma(k, 0.5*nu);
            
            V(:,:,1:T_IS) = [];
            loglik_vec(1:T_IS) = [];
            
        end
        
        
        
        
        function [V,Flag]  = IW_HAR(k, T,params,RC)
            beta_1 = params(1);
            beta_2 = params(2);
            beta_3 = params(3);
            
            Flag = 0;
            
            l_2 = 12;
            l_3 = 60;
            RC_bar = mean(RC,3);
            B0 = (1-sum(params))*RC_bar;
            V = zeros(k,k,T);
            V(:,:,1) = RC_bar;
            for i = 1:T-1
                
                check_ind = rcond(V(:,:,i));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                
                if i<l_2
                    V(:,:,i+1) = B0 + beta_1 * RC(:,:,i) + beta_2 * mean(RC(:,:,1:i),3) ...
                        + beta_3 * mean(RC(:,:,1:i),3);
                elseif i<l_3
                    
                    V(:,:,i+1) = B0 + beta_1 * RC(:,:,i) + beta_2 * mean(RC(:,:,i-l_2+1:i),3) ...
                        + beta_3 * mean(RC(:,:,1:i),3);
                else
                    V(:,:,i+1) = B0 + beta_1 * RC(:,:,i) + beta_2 * mean(RC(:,:,i-l_2+1:i),3)...
                        + beta_3 * mean(RC(:,:,i-l_3+1:i),3);
                end
            end
        end
        
        
        function [V,Flag]  = IW_HAR_CT_fcst(k, T,params,RC,T_IS)
            beta_1 = params(1);
            beta_2 = params(2);
            beta_3 = params(3);
            
            Flag = 0;
            
            l_2 = 12;
            l_3 = 60;
            RC_bar = mean(RC(:,:,1:T_IS),3);
            B0 = (1-sum(params))*RC_bar;
            V = zeros(k,k,T);
            V(:,:,1) = RC_bar;
            for i = 1:T-1
                
                check_ind = rcond(V(:,:,i));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                
                if i<l_2
                    V(:,:,i+1) = B0 + beta_1 * RC(:,:,i) + beta_2 * mean(RC(:,:,1:i),3) ...
                        + beta_3 * mean(RC(:,:,1:i),3);
                elseif i<l_3
                    
                    V(:,:,i+1) = B0 + beta_1 * RC(:,:,i) + beta_2 * mean(RC(:,:,i-l_2+1:i),3) ...
                        + beta_3 * mean(RC(:,:,1:i),3);
                else
                    V(:,:,i+1) = B0 + beta_1 * RC(:,:,i) + beta_2 * mean(RC(:,:,i-l_2+1:i),3)...
                        + beta_3 * mean(RC(:,:,i-l_3+1:i),3);
                end
            end
            
            V(:,:,1:T_IS) = [];
            
        end
        
        function [V,Flag]  = BEKK_FULL_CT( k, T,params, y_mat);
            
            alpha_vec = params(1:k+1);
            beta_vec = params(k+2:end);
            
            [i_LT] = Admin.UnitLT(k);
            A_mat = diag(alpha_vec(1:k)) + alpha_vec(end) * i_LT;
            B_mat = diag(beta_vec(1:k)) + beta_vec(end) * i_LT;
            
            %alpha = params(1);
            %beta  = params(2);
            
            Vbar  = cov(y_mat);
            Omega = Vbar - A_mat * Vbar * A_mat' - B_mat * Vbar * B_mat';
            V     = zeros(k,k,T);
            V(:,:,1) =  Vbar;
            Flag     = 0;
            for j = 1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                
                if j<T
                    %V(:,:,j+1) =  Omega + alpha*r2(:,:,j) + beta*V(:,:,j);
                    V(:,:,j+1) =  Omega + A_mat*(y_mat(j,:)'*y_mat(j,:))*A_mat' + B_mat*V(:,:,j)*B_mat';
                end
                
            end
            
            
            
        end
        
        
        function [V,Flag]  = BEKK_RES_CT( k, T,params, y_mat);
            
            alpha_vec = params(1:k+1);
            beta_vec = params(k+2:end);
            
            [i_LT] = Admin.UnitLT(k);
            A_mat = diag(alpha_vec(1:k)) + alpha_vec(end) * i_LT;
            B_mat = beta_vec(1)*eye(k) + beta_vec(2) * i_LT;
            
            %alpha = params(1);
            %beta  = params(2);
            
            Vbar  = cov(y_mat);
            Omega = Vbar - A_mat * Vbar * A_mat' - B_mat * Vbar * B_mat';
            V     = zeros(k,k,T);
            V(:,:,1) =  Vbar;
            Flag     = 0;
            for j = 1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                
                if j<T
                    %V(:,:,j+1) =  Omega + alpha*r2(:,:,j) + beta*V(:,:,j);
                    V(:,:,j+1) =  Omega + A_mat*(y_mat(j,:)'*y_mat(j,:))*A_mat' + B_mat*V(:,:,j)*B_mat';
                end
                
            end
            
            
            
        end
        
        function [V,Flag]  = BEKK_DIAG_CT( k, T,params, y_mat);
            
            alpha_vec = params(1:k);
            beta = params(k+1);
            
            A_mat = diag(alpha_vec(1:k));
            %B_mat = beta_vec(1)*eye(k) + beta_vec(2) * i_LT;
            
            %alpha = params(1);
            %beta  = params(2);
            
            Vbar  = cov(y_mat);
            Omega = Vbar - A_mat * Vbar * A_mat - beta*Vbar;
            V     = zeros(k,k,T);
            V(:,:,1) =  Vbar;
            Flag     = 0;
            for j = 1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                
                if j<T
                    %V(:,:,j+1) =  Omega + alpha*r2(:,:,j) + beta*V(:,:,j);
                    V(:,:,j+1) =  Omega + A_mat*(y_mat(j,:)'*y_mat(j,:))*A_mat + beta*V(:,:,j);
                end
                
            end
            
            
            
        end
        
        function [V,Flag]  = HEAVY_DIAG_CT( k, T,params, y_mat,RC_mat);
            
            alpha_vec = params(1:k);
            beta = params(k+1);
            
            A_mat = diag(alpha_vec(1:k));
            %B_mat = beta_vec(1)*eye(k) + beta_vec(2) * i_LT;
            
            RC_bar = mean(RC_mat,3);
            
            %alpha = params(1);
            %beta  = params(2);
            
            Vbar  = cov(y_mat);
            Omega = Vbar - A_mat * RC_bar * A_mat - beta*Vbar;
            V     = zeros(k,k,T);
            V(:,:,1) =  Vbar;
            Flag     = 0;
            for j = 1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                
                if j<T
                    %V(:,:,j+1) =  Omega + alpha*r2(:,:,j) + beta*V(:,:,j);
                    V(:,:,j+1) =  Omega + (A_mat*RC_mat(:,:,j))*A_mat + beta*V(:,:,j);
                end
                
            end
            
            
            
        end
        
        
        function [V,Flag]  = BEKK_SPILL_CT( k, T,params, y_mat);
            
            alpha_vec = params(1:k);
            beta = params(k+1);
            
            A_mat = alpha_vec*alpha_vec';
            %B_mat = beta_vec(1)*eye(k) + beta_vec(2) * i_LT;
            
            %alpha = params(1);
            %beta  = params(2);
            
            Vbar  = cov(y_mat);
            Omega = Vbar - A_mat * Vbar * A_mat - beta*Vbar;
            V     = zeros(k,k,T);
            V(:,:,1) =  Vbar;
            Flag     = 0;
            for j = 1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                
                if j<T
                    %V(:,:,j+1) =  Omega + alpha*r2(:,:,j) + beta*V(:,:,j);
                    V(:,:,j+1) =  Omega + A_mat*(y_mat(j,:)'*y_mat(j,:))*A_mat + beta*V(:,:,j);
                end
                
            end
            
            
            
        end
        
        function [V_mat,Flag]  = BEKK_CT( k, T,params, y_mat);
            
            alpha = params(1);
            beta  = params(2);
            
            Vbar  = cov(y_mat);
            Omega = (1-alpha-beta)*Vbar;
            V_mat = zeros(k,k,T);
            V_mat(:,:,1) = Vbar;
            
            Flag     = 0;
            for j = 1:T
                
                check_ind = rcond(V_mat(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                
                if j<T
                    %V(:,:,j+1) =  Omega + alpha*r2(:,:,j) + beta*V(:,:,j);
                    V_mat(:,:,j+1) =  Omega + alpha*(y_mat(j,:)'*y_mat(j,:)) + beta*V_mat(:,:,j);
                end
                
            end
            
            
            
        end
        
        
        function [V,Flag]  = HEAVY_CT( k, T,params, y_mat,RC_mat);
            
            alpha = params(1);
            beta  = params(2);
            
            Vbar  = cov(y_mat);
            %Omega = (1-alpha-beta)*Vbar;
            
            Omega = (1-beta)*Vbar - alpha*mean(RC_mat,3);
            V     = zeros(k,k,T);
            V(:,:,1) =  Vbar;
            Flag     = 0;
            for j = 1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                if j<T
                    V(:,:,j+1) =  Omega + alpha*RC_mat(:,:,j) + beta*V(:,:,j);
                end
                
            end
            
            
            
        end
                                     
        
        function [V_mat_vech,Flag,V_mat_vech_IS,V_mat_vech_TOT]  = HEAVY_CT_fcst( k, T,params, y_mat,RC_mat,T_IS);
            
            alpha = params(1);
            beta  = params(2);
            
            Vbar  = cov(y_mat(1:T_IS,:));
            %Omega = (1-alpha-beta)*Vbar;
            
            Omega = (1-beta)*Vbar - alpha*mean(RC_mat(:,:,1:T_IS),3);
            
            nr_vech = k*(k+1)/2;
            V_old =  Vbar;
            V_mat_vech     = zeros(T,nr_vech);
            Flag     = 0;
            
            for j = 1:T
                
                check_ind = rcond(V_old);
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                %if j>T_IS
                V_mat_vech(j,:) = Admin.Sym2Vech(k,V_old);
                %end
                
                if j<T
                    V_new =  Omega + alpha*RC_mat(:,:,j)+ beta*V_old;
                    V_old = V_new;
                end
                
            end
            
            V_mat_vech_IS = V_mat_vech(1:T_IS,:);
            V_mat_vech_TOT = V_mat_vech;
            V_mat_vech(1:T_IS,:) = [];
            
        end                        
        
        
        function [V_mat_vech,Flag,V_mat_vech_IS,V_mat_vech_TOT]  = HEAVY_LT_spec_fcst(k, T,params,RC_mat,T_IS);
            
            alpha = params(1);
            beta  = params(2);
            
            I_k = eye(k);
            
            %Vbar  = cov(y_mat(1:T_IS,:));
            V_bar = mean(RC_mat(:,:,1:T_IS),3);
            chol_H_LT_bar = chol(V_bar)';
            inv_chol_H_LT_bar = chol_H_LT_bar\I_k;
                                       
            %Omega = (1-beta)*Vbar - alpha*mean(RC_mat(:,:,1:T_IS),3);
            
            Omega_ST = (1-alpha-beta)*I_k;
            
            nr_vech = k*(k+1)/2;
            V_old =  V_bar;
            V_mat_vech     = zeros(T,nr_vech);
            Flag     = 0;
                              
            
            for j = 1:T
                
                H_TOT = chol_H_LT_bar  *V_old * chol_H_LT_bar';
                
                check_ind = rcond(H_TOT);
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                %if j>T_IS
                V_mat_vech(j,:) = Admin.Sym2Vech(k,H_TOT);
                %end
                
                if j<T
                    V_new = Omega_ST + alpha*inv_chol_H_LT_bar*(RC_mat(:,:,j)*inv_chol_H_LT_bar') + beta*V_old;
                    V_old = V_new;
                end
                
            end
            
            V_mat_vech_IS = V_mat_vech(1:T_IS,:);
            V_mat_vech_TOT = V_mat_vech;
            V_mat_vech(1:T_IS,:) = [];
            
        end                        
        
        
        
        function [V_mat_vech,Flag]  = BEKK_CT_fcst( k, T,params, y_mat,T_IS);
            
            alpha = params(1);
            beta  = params(2);
            
            Vbar  = cov(y_mat(1:T_IS,:));
            Omega = (1-alpha-beta)*Vbar;
            
            nr_vech = k*(k+1)/2;
            V_old =  Vbar;
            V_mat_vech     = zeros(T-T_IS,nr_vech);
            
            Flag     = 0;
            for j = 1:T
                
                check_ind = rcond(V_old);
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                if j>T_IS
                    V_mat_vech(j-T_IS,:) = Admin.Sym2Vech(k,V_old);
                end
                
                
                if j<T
                    %V(:,:,j+1) =  Omega + alpha*r2(:,:,j) + beta*V(:,:,j);
                    V_new =  Omega + alpha*(y_mat(j,:)'*y_mat(j,:)) + beta*V_old;
                    V_old = V_new;
                end
                
            end
            
            
            
        end
        
        function [V_CTC_mat_vech,c_vec,Flag]=Mult_T_tv_c_given_V_DAY_fcst(k, T,params, y_mat,V_DAY_vech,T_IS,c_1);
            
            c_vec = zeros(T,1);
            c_vec(1) = c_1;
            
            Flag = 0;
            
            omega = params(1);
            alfa = params(2);
            beta = params(3);
            
            nu_t = params(end);
            nr_k = k*(k+1)/2;
            V_CTC_mat_vech = zeros(T,nr_k);
            
            
            for j =1 :T
                V_CTC_t  = c_vec(j)*Admin.Vech2Sym(k,V_DAY_vech(j,:));
                
                V_CTC_mat_vech(j,:) = c_vec(j)*V_DAY_vech(j,:);
                
                if isnan(rcond(V_CTC_t))
                    Flag = 1;
                    break
                end
                
                inv_Sigma_t = V_CTC_t\eye(k);
                y_V_inv_y = y_mat(j,:)*inv_Sigma_t*y_mat(j,:)';
                
                if j<T
                    w_t = (nu_t+k)/(nu_t-2 + y_V_inv_y);
                    nabla_t = 0.5 * ((w_t/c_vec(j)) *  y_V_inv_y - k/c_vec(j));
                    s_t = c_vec(j)^2 * nabla_t;
                    c_vec(j+1) = omega + alfa * s_t + beta * c_vec(j);
                end
                
              
                
            end
            
           
            
            V_CTC_mat_vech(1:T_IS,:) = [];
            c_vec(1:T_IS) = [];
            
            
            
        end
        
        
        function [V_CTC_mat_vech,c_mat,Flag]=Mult_T_tv_c_vec_given_V_DAY_fcst(k, T,params, y_mat,V_DAY_vech,T_IS,c_1);
            
            c_mat = zeros(T,k);
            c_mat(1,:) = c_1;
            
            Flag = 0;
            
            omega_vec = params(1:k)';
            alfa = params(k+1);
            beta = params(k+2);
            
            nu_t = params(end);
            
            nr_k_vech = k*(k+1)/2;
            V_CTC_mat_vech = zeros(T,nr_k_vech);
            
            
            for j =1 :T
                
                c_vec_j = c_mat(j,:);
                C_mat_j = diag(c_vec_j);
                V_CTC_t  = C_mat_j*Admin.Vech2Sym(k,V_DAY_vech(j,:)) * C_mat_j;
                
                V_CTC_mat_vech(j,:) = Admin.Sym2Vech(k,V_CTC_t);
                
                if isnan(rcond(V_CTC_t))
                    Flag = 1;
                    break
                end
                
                inv_Sigma_t = V_CTC_t\eye(k);
                y_V_inv_y = y_mat(j,:)*inv_Sigma_t*y_mat(j,:)';
                
                
                if j<T
                    w_t = (nu_t+k)/(nu_t-2 + y_V_inv_y);
                                        
                    inv_C_t = diag(1./ c_vec_j);
                    A_t = inv_C_t *  (y_mat(j,:)'*y_mat(j,:))*inv_Sigma_t;
                    A_t_trans = A_t';
                    
                    nabla_t   = 0.5 * w_t * diag(A_t + A_t_trans)' - 1./c_vec_j;
                    s_t = c_mat(j,:).^2 .* nabla_t;
                    
                    c_mat(j+1,:) = omega_vec + alfa * s_t + beta * c_mat(j,:);
                end
                
                
            end
                                    
            V_CTC_mat_vech(1:T_IS,:) = [];
            c_mat(1:T_IS,:) = [];
                                    
        end
        
        
        function [V_CTC_mat_vech,c_vec,Flag]=Mult_Tr_tv_c_given_V_DAY_fcst(k, T,params, y_mat,V_DAY_vech,T_IS,c_1);
            
            c_vec = zeros(T,1);
            c_vec(1) = c_1;
            
            Flag = 0;
            
            omega = params(1);
            alfa = params(2);
            beta = params(3);
            
            nu_vec = params(4:end);
           
            nr_k_vech = k*(k+1)/2;
            V_CTC_mat_vech = zeros(T,nr_k_vech);
                                    
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
            
            
            for j =1 :T
                
                    V_DAY_t = Admin.Vech2Sym(k,V_DAY_vech(j,:));
                                        
                    V_CTC_t  = c_vec(j)* V_DAY_t;                
                    V_CTC_mat_vech(j,:) =  c_vec(j)*V_DAY_vech(j,:);                                            
                    Sigma_inv_j_0    = V_CTC_t\eye(k);
                    
                    [L_0_trans,ind_chol] = chol(Sigma_inv_j_0);
                    if ind_chol>0
                        Flag =1;
                        break;
                    end
                    
                    L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                    inv_Sigma_j = L_trans' * L_trans;
                    
                    [~,ind_chol] = chol(inv_Sigma_j);
                    if ind_chol>0
                        Flag =1;
                        break;
                    end
                                                                  
                    j_Sigma = inv(inv_Sigma_j);                    
                    j_yy       = y_mat(j,:)'*y_mat(j,:);
                    j_yy_Sigma = j_yy + j_Sigma;
                                                                                
                    if rcond(j_yy_Sigma)<1e-15
                        Flag =1;
                        break
                    end                                        
                                                                                                                      
                    j_yy_Sigma_inv = j_yy_Sigma\eye(k);
                    
                    [~,ind_chol] = chol(j_yy_Sigma_inv);
                    if ind_chol>0
                        Flag =1;
                        break
                    end
                                        
                    if j<T
                       L_inv_S  = chol(j_yy_Sigma_inv);
                       U_t      = inv(L_inv_S);                       
                       diag_U_t_vec = diag(U_t);
                                                                 
                       s_t_1 = 0.5*sum(nu_vec)/c_vec(j);
                       
                       dau_ln_PWD_dau_diag_U = -(nu_vec'+1)./diag_U_t_vec';                                                                                                                              
                       
                       Z = kron(I_k,U_t) + kron(U_t,I_k)*K_k;                       
                       temp_score_2                 = inv_D_prime_D_D_prime * Z*D_tilde_k;                                                                    
                       dau_vech_U_prime_dau_vech_S  = temp_score_2\I_k_vech;                   
                       dau_vech_S_dau_vec_S         = inv_D_prime_D_D_prime;                                              
                       dau_vec_S_dau_c              = j_Sigma(:)/c_vec(j);
                                                                     
                       s_t = c_vec(j)^2 * (s_t_1 + dau_ln_PWD_dau_diag_U*dau_diag_U_dau_vech_U_prime*dau_vech_U_prime_dau_vech_S...
                           *dau_vech_S_dau_vec_S*dau_vec_S_dau_c);

                     % s_vec(j) = s_t;
                       
                       c_vec(j+1) = omega + alfa * s_t + beta * c_vec(j);
                    end                                                
                                    
            end
            
            V_CTC_mat_vech(1:T_IS,:) = [];
            c_vec(1:T_IS) = [];
        
        end
        
        
        
        function [V_CTC_mat_vech,c_mat,Flag]=Mult_Tr_tv_c_vec_given_V_DAY_fcst(k, T,params, y_mat,V_DAY_vech,T_IS,c_1);
            
            c_mat = zeros(T,k);
            c_mat(1,:) = c_1;
            
            Flag = 0;
            
            omega = params(1:k)';
            alfa = params(k+1);
            beta = params(k+2);            
            nu_vec = params(k+3:end);
           
            nr_k_vech = k*(k+1)/2;
            V_CTC_mat_vech = zeros(T,nr_k_vech);
                                    
            mu_vec = Admin.Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            K_k = Admin.commutation(k,k);
            D_k = Admin.duplication(k);
            D_tilde_k = Admin.duplication_gen(k);
            I_k = eye(k);
            I_k_vech = eye(k*(k+1)/2);
            
            inv_D_prime_D_D_prime    = inv(D_k'*D_k) *D_k';                   
            K_k_D_tilde = K_k*D_tilde_k;
            S_D_tilde = Admin.Compute_S_d_matrix(1:k)';            
            dau_diag_U_dau_vech_U_prime = S_D_tilde*K_k_D_tilde;
            dau_vec_C_dau_diag_C = S_D_tilde';   
            
            for j =1 :T
                
                    V_DAY_t = Admin.Vech2Sym(k,V_DAY_vech(j,:));
                    Sigma_inv_j_0    = V_DAY_t\I_k;
                    
                    C_mat = diag(c_mat(j,:));
                    
                    V_CTC_t = C_mat * V_DAY_t * C_mat;
                    V_CTC_mat_vech(j,:) = Admin.Sym2Vech(k,V_CTC_t);
             
                    L_0_trans = chol(Sigma_inv_j_0);                                      
                    L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                    inv_Sigma_j = L_trans' * L_trans;
                    
                    [~,ind_chol] = chol(inv_Sigma_j);
                    if ind_chol>0
                        Flag = 1;
                        break;
                    end
                              
                    j_Sigma_DAY = inv_Sigma_j\I_k;                    
                    j_Sigma_CTC = C_mat* j_Sigma_DAY * C_mat; 
                    
                    if rcond(j_Sigma_CTC)<1e-15
                        Flag = 1;
                        break;
                    end                                                                              
                                                         
                    j_yy       = y_mat(j,:)'*y_mat(j,:);                               
                    j_yy_Sigma = j_yy + j_Sigma_CTC;
                    
                    if rcond(j_yy_Sigma)<1e-15
                        Flag = 1;
                        break;
                    end                                                                                                                                                              
                                                            
                    j_yy_Sigma_inv = j_yy_Sigma\eye(k);
                    
                    [~,ind_chol] = chol(j_yy_Sigma_inv);
                    if ind_chol>0
                        Flag = 1;
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

                       %s_mat(j,:) = s_t;
                       
                       c_mat(j+1,:) = omega + alfa * s_t + beta * c_mat(j,:);
                    end                                                                  
                                    
            end
            
            V_CTC_mat_vech(1:T_IS,:) = [];
            c_mat(1:T_IS,:) = [];
        
        end
        
        
        function port_return_mat_sim  = BEKK_CT_SIM(k,h_step,params,w_vec,Vbar,Vt1,z_cell_N,ind_t_dist,N_sim);
            
            alpha = params(1);
            beta  = params(2);
            
            if ind_t_dist==1
                nu = params(3);
                inv_gamma_vec = sqrt((nu-2)./ (2.*randg(nu./2,N_sim,1)));
            end
            
            Omega = (1-alpha-beta)*Vbar;
            
            z_t_sim = zeros(N_sim,k);
            port_return_mat_sim = zeros(N_sim,h_step);
            Vnew = zeros(k,k,N_sim);
            
            for j = 1:h_step
                
                if j ==1
                    if ind_t_dist == 1
                        z_t_sim     =  inv_gamma_vec*ones(1,k).*(z_cell_N{j}*chol(Vt1));
                    else
                        z_t_sim     =  z_cell_N{j}*chol(Vt1);
                    end
                    
                    port_return_mat_sim(:,j) = z_t_sim*w_vec;
                    
                    for o = 1: N_sim
                        Vnew(:,:,o) = Omega + alpha * (z_t_sim(o,:)'*z_t_sim(o,:)) + beta *Vt1;
                    end
                    
                else
                    if ind_t_dist ==1
                        for o = 1: N_sim
                            z_t_sim(o,:) =inv_gamma_vec(o)*(z_cell_N{j}(o,:)*chol(Vold(:,:,o)));
                            Vnew(:,:,o) = Omega + alpha * (z_t_sim(o,:)'*z_t_sim(o,:)) + beta *Vold(:,:,o);
                        end
                        
                        port_return_mat_sim(:,j) = z_t_sim*w_vec;
                    else
                        for o = 1: N_sim
                            z_t_sim(o,:) = z_cell_N{j}(o,:)*chol(Vold(:,:,o));
                            Vnew(:,:,o) = Omega + alpha * (z_t_sim(o,:)'*z_t_sim(o,:)) + beta *Vold(:,:,o);
                        end
                        
                        port_return_mat_sim(:,j) = z_t_sim*w_vec;
                    end
                end
                
                
                Vold = Vnew;
                
            end
            
            
        end
        
        
        function port_return_mat_sim  = BEKK_CT_TRIESZ_SIM(k,h_step,params,w_vec,Vbar,Vt1,z_cell_N,X_iRiesz_cell,N_sim);
            
            alpha = params(1);
            beta  = params(2);
            
            Omega = (1-alpha-beta)*Vbar;
            
            z_t_sim = zeros(N_sim,k);
            port_return_mat_sim = zeros(N_sim,h_step);
            Vnew = zeros(k,k,N_sim);
            
            for j = 1:h_step
                
                if j ==1
                    
                    V_inv_t = Vt1\eye(k);
                    L_inv_V_t = chol(V_inv_t)';
                    L_V_t = L_inv_V_t\eye(k);
                    
                    
                    for o = 1:N_sim
                        W_ps_o = L_V_t'*X_iRiesz_cell{j}(:,:,o)*L_V_t;
                        z_t_sim(o,:) = z_cell_N{j}(o,:)*chol(W_ps_o);
                    end
                    
                    
                    port_return_mat_sim(:,j) = z_t_sim*w_vec;
                    
                    for o = 1: N_sim
                        Vnew(:,:,o) = Omega + alpha * (z_t_sim(o,:)'*z_t_sim(o,:)) + beta *Vt1;
                    end
                    
                else
                    for o = 1: N_sim
                        
                        V_inv_t = Vold(:,:,o)\eye(k);
                        L_inv_V_t = chol(V_inv_t)';
                        L_V_t = L_inv_V_t\eye(k);
                        
                        W_ps_i = L_V_t'*X_iRiesz_cell{j}(:,:,o)*L_V_t;
                        z_t_sim(o,:) = z_cell_N{j}(o,:)*chol(W_ps_i);
                        
                        Vnew(:,:,o) = Omega + alpha * (z_t_sim(o,:)'*z_t_sim(o,:)) + beta *Vold(:,:,o);
                    end
                    
                    port_return_mat_sim(:,j) = z_t_sim*w_vec;
                    
                end
                
                
                Vold = Vnew;
                
            end
            
            
        end
        
        
        function port_return_mat_sim  = BEKK_CT_DIAG_TRIESZ_SIM(k,h_step,params,w_vec,Vbar,Vt1,z_cell_N,X_iRiesz_cell,N_sim);
            
            A_mat = diag(params(1:k));
            beta  = params(k+1);
            
            Omega = (1-beta)*Vbar - A_mat * Vbar * A_mat;
            
            z_t_sim = zeros(N_sim,k);
            port_return_mat_sim = zeros(N_sim,h_step);
            Vnew = zeros(k,k,N_sim);
            
            for j = 1:h_step
                
                if j ==1
                    
                    V_inv_t = V_t1\eye(k);
                    L_inv_V_t = chol(V_inv_t)';
                    L_V_t = L_inv_V_t\eye(k);
                    
                    
                    for o = 1:N_sim
                        W_ps_o = L_V_t'*X_iRiesz_cell{j}(:,:,o)*L_V_t;
                        z_t_sim(o,:) = z_cell_N{j}(o,:)*chol(W_ps_o);
                    end
                    
                    
                    port_return_mat_sim(:,j) = z_t_sim*w_vec;
                    
                    for o = 1: N_sim
                        Vnew(:,:,o) = Omega +  A_mat *  (z_t_sim(o,:)'*z_t_sim(o,:))*A_mat + beta *Vt1;
                    end
                    
                else
                    for o = 1: N_sim
                        
                        V_inv_t = Vold(:,:,o)\eye(k);
                        L_inv_V_t = chol(V_inv_t)';
                        L_V_t = L_inv_V_t\eye(k);
                        
                        W_ps_i = L_V_t'*X_iRiesz_cell{j}(:,:,o)*L_V_t;
                        z_t_sim(o,:) = z_cell_N{j}(o,:)*chol(W_ps_i);
                        
                        Vnew(:,:,o) = Omega + A_mat * (z_t_sim(o,:)'*z_t_sim(o,:))* A_mat + beta *Vold(:,:,o);
                    end
                    
                    port_return_mat_sim(:,j) = z_t_sim*w_vec;
                    
                end
                
                
                Vold = Vnew;
                
            end
            
            
        end
        
        function port_return_mat_sim  = BEKK_CT_DIAG_SIM(k,h_step,params,w_vec,Vbar,Vt1,z_cell_N,ind_t_dist,N_sim);
            
            alpha_vec = params(1:k);
            beta  = params(k+1);
            A_mat = diag(alpha_vec);
            
            if ind_t_dist==1
                nu = params(3);
                inv_gamma_vec = sqrt((nu-2)./ (2.*randg(nu./2,N_sim,1)));
            end
            
            Omega = (1-beta)*Vbar - A_mat * Vbar * A_mat;
            
            z_t_sim = zeros(N_sim,k);
            port_return_mat_sim = zeros(N_sim,h_step);
            Vnew = zeros(k,k,N_sim);
            
            for j = 1:h_step
                
                if j ==1
                    if ind_t_dist == 1
                        z_t_sim     =  inv_gamma_vec*ones(1,k).*(z_cell_N{j}*chol(Vt1));
                    else
                        z_t_sim     =  z_cell_N{j}*chol(Vt1);
                    end
                    
                    port_return_mat_sim(:,j) = z_t_sim*w_vec;
                    
                    for o = 1: N_sim
                        Vnew(:,:,o) = Omega + A_mat  * (z_t_sim(o,:)'*z_t_sim(o,:))*A_mat  + beta *Vt1;
                    end
                    
                else
                    if ind_t_dist ==1
                        for o = 1: N_sim
                            z_t_sim(o,:) =inv_gamma_vec(o)*(z_cell_N{j}(o,:)*chol(Vold(:,:,o)));
                            Vnew(:,:,o) = Omega + A_mat * (z_t_sim(o,:)'*z_t_sim(o,:))*A_mat + beta *Vold(:,:,o);
                        end
                        
                        port_return_mat_sim(:,j) = z_t_sim*w_vec;
                    else
                        for o = 1: N_sim
                            z_t_sim(o,:) = z_cell_N{j}(o,:)*chol(Vold(:,:,o));
                            Vnew(:,:,o) = Omega + A_mat * (z_t_sim(o,:)'*z_t_sim(o,:))*A_mat + beta *Vold(:,:,o);
                        end
                        
                        port_return_mat_sim(:,j) = z_t_sim*w_vec;
                    end
                end
                
                
                Vold = Vnew;
                
            end
            
            
        end
        
        
        
        
        function [V,Flag]  = BEKK_DIAG_CT_fcst( k, T,params, return_mat,T_IS);
            
            alpha_vec = params(1:k);
            beta  = params(k+1);
            
            A_mat = diag(alpha_vec);
            
            Vbar  = cov(return_mat(1:T_IS,:));
            Omega = (1-beta)*Vbar - A_mat*Vbar*A_mat;
            V     = zeros(k,k,T);
            V(:,:,1) =  Vbar;
            Flag     = 0;
            for j = 1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                
                if j<T
                    V(:,:,j+1) =  Omega + A_mat*(return_mat(j,:)'*return_mat(j,:))*A_mat + beta*V(:,:,j);
                end
                
            end
            
            V(:,:,1:T_IS) = [];
            
        end
        
        
        function [V,Flag]  = HEAVY_DIAG_CT_fcst( k, T,params, return_mat,RC_mat,T_IS);
            
            alpha_vec = params(1:k);
            beta  = params(k+1);
            
            A_mat = diag(alpha_vec);
            
            Vbar  = cov(return_mat(1:T_IS,:));
            RC_bar= mean(RC_mat(:,:,1:T_IS),3);
            
            Omega = (1-beta)*Vbar - A_mat*RC_bar*A_mat;
            V     = zeros(k,k,T);
            V(:,:,1) =  Vbar;
            Flag     = 0;
            for j = 1:T
                
                check_ind = rcond(V(:,:,j));
                if check_ind < 1e-15
                    Flag = 1;
                    break
                end
                
                
                if j<T
                    V(:,:,j+1) =  Omega + A_mat*(RC_mat(:,:,j)*A_mat) + beta*V(:,:,j);
                end
                
            end
            
            V(:,:,1:T_IS) = [];
            
        end
        
        
        function [Rt_mat,loglike_vec,C_q_mat,Ind_tail_mat,Qbar_cDCC,Qt_mat] = cDCC_Copula(T,params, u_mat,ind_t_dist,T_IS,ind_csl,q_vec,M)
            
            if ind_t_dist ==1
                nu = params(end);
                Stdresid_0 = tinv(u_mat,nu);
                Stdresid = sqrt((nu-2)/nu) * Stdresid_0;
                clear Stdresid_0
            else
                Stdresid = norminv(u_mat);
            end
            
            k= size(u_mat,2);
            Qbar_DCC = cov(Stdresid(1:T_IS,:));
            
            Rt_mat = zeros(k,k,T-T_IS);
            Qt_mat = zeros(k,k,T-T_IS);
            
            Qt_current = Qbar_DCC;
            alpha = params(1);
            beta =  params(2);
            loglike_vec = zeros(T,1);
            
            C_q_mat = [];
            Ind_tail_mat = [];
            
            if ind_csl ==1
                C_q_mat = nan*ones(T,length(q_vec));
                Ind_tail_mat = zeros(T,length(q_vec));
            end
            
            q_ii_mat = zeros(T_IS,k);
            q_ii_mat(1,:) = diag(Qbar_DCC);
            Qbar_new  = zeros(k,k);
            for m = 1:T_IS
                if m<T_IS
                    q_ii_mat(m+1,:) = (1-alpha-beta) + alpha*q_ii_mat(m,:).*Stdresid(m,:).^2 + beta*q_ii_mat(m,:);
                end
                Q_square_root_star_t = diag(sqrt(q_ii_mat(m,:)));
                Qbar_new = Qbar_new +  Q_square_root_star_t * Stdresid(m,:)'* Stdresid(m,:) * Q_square_root_star_t;
            end
            
            Qbar_cDCC = Qbar_new/T_IS;
            
            for j = 1:T
                
                Rt_DCC= Qt_current./(sqrt(diag(Qt_current))*sqrt(diag(Qt_current))');
                
                if j<T
                    temp = sqrt(diag(Qt_current));
                    Q_star_t = diag(temp(1:end));
                    Qt_new = Qbar_cDCC*(1 - alpha-beta) + alpha*(Q_star_t*Stdresid(j,:)'*Stdresid(j,:)*Q_star_t) + beta*Qt_current;
                end
                
                if j>T_IS
                    Rt_mat(:,:,j-T_IS) = Rt_DCC;
                    Qt_mat(:,:,j-T_IS) = Qt_current;
                end
                
                if j>T_IS
                    x_R_inv_x = Stdresid(j,:)*inv(Rt_DCC)*Stdresid(j,:)';
                    log_det_Rt   = log(det(Rt_DCC));
                    if ind_t_dist ==1
                        % tdist
                        loglike_vec(j,1) = -0.5 *(nu+k)*log(1 + x_R_inv_x/(nu-2)) - 0.5* log_det_Rt;
                    else
                        loglike_vec(j,1) = -0.5 * x_R_inv_x - 0.5*log_det_Rt;
                    end
                    
                    if ind_csl==1
                        
                        for m = 1:length(q_vec)
                            q_i = q_vec(m)*ones(k,1);
                            if sum(Stdresid(j,:)'<=q_i)==k
                                Ind_tail_mat(j,m) = 1;
                            end
                            % now compute C(q)
                            if ind_t_dist ==1
                                C_q_mat(j,m) = qsimvt(M,nu,Rt_DCC,-inf*ones(k,1),q_i);
                            else
                                C_q_mat(j,m) = qsimvn(M,Rt_DCC,-inf*ones(k,1),q_i);
                            end
                        end
                    end
                end
                
                Qt_current = Qt_new;
            end
            
            if ind_t_dist==1
                loglike_vec = loglike_vec + gammaln(0.5*(nu+k)) + (k-1)*gammaln(nu/2) - k * gammaln(0.5*(nu+1)) ...
                    + 0.5*(nu+1)*sum(log(ones(T,k)+(1/(nu-2))*Stdresid.^2),2);
            else
                loglike_vec = loglike_vec + 0.5*sum(Stdresid.^2,2);
            end
            
            
            loglike_vec(1:T_IS) = [];
            if ind_csl==1
                C_q_mat(1:T_IS,:) = [];
                Ind_tail_mat(1:T_IS,:) = [];
            end
            
            
        end
        
        
        
        function port_return_sim = cDCC_t_Copula_PLUS_MARG_SIM(k,h_step,params_DCC,params_GARCH,Qt_current,Qbar_cDCC,h_vec,z_cell_N,w_vec,N_sim)
            % marginals: t_aggr
            
            
            nu_vec_i = params_GARCH(5:end)';
            nu_Cop_t = params_DCC(end);
            inv_gamma_vec = sqrt(nu_Cop_t ./ (2.*randg(nu_Cop_t./2,N_sim,1)));
            
            
            alpha = params_DCC(1);
            beta =  params_DCC(2);
            Qt_new = zeros(k,k,N_sim);
            
            port_return_sim = zeros(N_sim,h_step);
            
            for j = 1:h_step
                
                if j ==1
                    Rt_DCC= Qt_current./(sqrt(diag(Qt_current))*sqrt(diag(Qt_current))');
                    z_t_sim_c     =  inv_gamma_vec*ones(1,k).*(z_cell_N{j}*chol(Rt_DCC));
                    
                    temp = sqrt(diag(Qt_current));
                    Q_star_t = diag(temp(1:end));
                    
                    for o = 1:N_sim
                        Qt_new(:,:,o) = Qbar_cDCC*(1 - alpha-beta) + alpha*(Q_star_t*z_t_sim_c(o,:)'*z_t_sim_c(o,:)*Q_star_t) + beta*Qt_current;
                    end
                    
                    u_sim       = tcdf(z_t_sim_c,nu_Cop_t);
                    z_t_sim     = ones(N_sim,1)*sqrt((nu_vec_i-2)./nu_vec_i).*tinv(u_sim,ones(N_sim,1)*nu_vec_i);
                    
                    D_mat = diag(sqrt(h_vec));
                    return_mat_sim_j = z_t_sim * D_mat;
                    port_return_sim(:,j) = return_mat_sim_j*w_vec;
                    
                    h_mat_new = params_GARCH(2) + params_GARCH(3)*z_t_sim.^2 + params_GARCH(4)*(ones(N_sim,1)*h_vec);
                    
                    
                    
                else
                    
                    for o = 1:N_sim
                        Rt_DCC_o           = Qt_current(:,:,o)./(sqrt(diag(Qt_current(:,:,o)))*sqrt(diag(Qt_current(:,:,o)))');
                        z_t_sim_c(o,:)     =  inv_gamma_vec(o)*(z_cell_N{j}(o,:)*chol(Rt_DCC_o));
                        Qt_new(:,:,o)      = Qbar_cDCC*(1 - alpha-beta) + alpha*(Q_star_t*z_t_sim_c(o,:)'*z_t_sim_c(o,:)*Q_star_t) + beta*Qt_current(:,:,o);
                        
                    end
                    
                    u_sim       = tcdf(z_t_sim_c,nu_Cop_t);
                    z_t_sim     = ones(N_sim,1)*sqrt((nu_vec_i-2)./nu_vec_i).*tinv(u_sim,ones(N_sim,1)*nu_vec_i);
                    
                    return_mat_sim_j     = sqrt(h_mat_old).*z_t_sim;
                    port_return_sim(:,j) = return_mat_sim_j*w_vec;
                    
                    h_mat_new = params_GARCH(2) + params_GARCH(3)*z_t_sim.^2 + params_GARCH(4)*h_mat_old;
                    
                end
                
                h_mat_old  = h_mat_new;
                Qt_current = Qt_new;
            end
            
            
        end
        
        
        
        function [Rt_mat,loglik_r2,Qbar_cDCC] = cDCC_TRIESZ_fcst(k,T,params,Stdresid,ind_Rt,T_IS)
            
            Qbar_DCC = cov(Stdresid(1:T_IS,:));
            D_mat = diag(sqrt(diag(Qbar_DCC)));
            
            if ind_Rt==1
                Rt_mat = zeros(k,k,T-T_IS);
            else
                Rt_mat =[];
            end
            
            
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
                
                if j>T_IS
                    Rt_mat(:,:,j-T_IS) = Rt_DCC;
                end
                
                Ht_DCC = D_mat*Rt_DCC*D_mat;
                %                 [~,ind_chol] = chol(Ht_DCC);
                %                 if ind_chol>0
                %                     loglik_r2 = NaN;
                %                     break;
                %                 end
                %
                
                H_inv_j_0    = Ht_DCC\eye(k);
                
                [L_0_trans,ind_chol] = chol(H_inv_j_0);
                %                 if ind_chol>0
                %                     loglik_r2 = NaN;
                %                     break;
                %                 end
                
                L_trans = diag(sqrt(mu_vec)) * L_0_trans;
                inv_H_j = L_trans' * L_trans;
                
                %                [~,ind_chol] = chol(inv_H_j);
                %                 if ind_chol>0
                %                     loglik_r2 = NaN;
                %                     break;
                %                 end
                
                
                j_H = inv(inv_H_j);
                j_yy       = Stdresid(j,:)'*Stdresid(j,:);
                j_yy_H = j_yy + j_H;
                
                ln_PWD_H = LogLik.ln_PWD(inv_H_j,0.5*nu_vec);
                j_yy_H_inv = j_yy_H\eye(k);
                
                %                [~,ind_chol] = chol(j_yy_H_inv);
                %                 if ind_chol>0
                %                     loglik_r2= NaN;
                %                     break
                %                 end
                
                ln_PWD_j_yy_H  = LogLik.ln_PWD(j_yy_H_inv,0.5*(nu_vec+1));
                loglik_r2(j) = -ln_PWD_H  + ln_PWD_j_yy_H;
                
                Qt_current = Qt_new;
                
            end
            
            loglik_r2 = loglik_r2 - (k/2)*log(pi) + LogLik.lnmultigamma_nu_vec(k, (nu_vec+1)/2 )   -  LogLik.lnmultigamma_nu_vec(k, nu_vec/2 );
            
            loglik_r2(1:T_IS) = [];
            
        end
        
        
        
        
        function [V] = EWMA_Uni(T, params,RV)
            
            V0  = mean(RV,1);
            b_EWMA = params;
            V = zeros(T,1);
            V(:,1) = V0;
            for j = 1:T
                if j<T
                    V(j+1,1) = b_EWMA * V(j,1) + (1-b_EWMA)*RV(j);
                end
            end
            
            
            
        end
        
        
        function [V] = EWMA_fcst(k, T, params,RC,T_IS)
            
            V0  = mean(RC(:,:,1:T_IS), 3);
            b_EWMA = params;
            V = zeros(k,k,T);
            V(:,:,1) = V0;
            for j = 1:T
                if j<T
                    V(:,:,j+1) = b_EWMA * V(:,:,j) + (1-b_EWMA)*RC(:,:,j);
                end
            end
            V(:,:,1:T_IS) =[];
            
            
        end
        
        
        function [h_vec,loglik_r2,u_vec,z_vec] =  GARCH_t_fcst(T,params,ret_vec,T_IS);
            
            % garch (verschillende opties mogelijk) model met scheve t verdeling
            c     = params(1);
            omega = params(end-3);
            alpha = params(end-2);
            beta  = params(end-1);
            
            nu_t = params(end);
            eps_vec = ret_vec - c;
            
            
            eps2_vec =eps_vec.^2;
            h_vec = zeros(T,1);
            h_vec(1) = var(ret_vec(1:T_IS));
            for i = 1:T+1
                if i<T
                    h_vec(i+1) = omega  + alpha * eps2_vec(i) + beta*h_vec(i);
                end
            end
            % h_last = h_vec(end);
            
            z_vec = eps_vec./sqrt(h_vec);
            ksi_vec = sqrt(nu_t/(nu_t-2))*z_vec;
            u_vec = tcdf(ksi_vec,nu_t);
            
            loglik_r2 = -0.5 * log(h_vec) - 0.5*(nu_t+1)*log(1 + eps2_vec./((nu_t-2)*h_vec));
            loglik_r2 = loglik_r2 + gammaln(0.5*(nu_t+1)) - gammaln(0.5*nu_t) - 0.5 * log(pi*(nu_t-2));
            
            h_vec(1:T_IS) = [];
            loglik_r2(1:T_IS) = [];
            
        end
        
        function [h_mat_fcst,loglik_mat_fcst,u_mat] =  GARCH_AGG_t_fcst(T,params,return_mat,T_IS);
            
            k = size(return_mat,2);
            h_mat_fcst = zeros(T-T_IS,k);
            u_mat = zeros(T,k);
            loglik_mat_fcst = zeros(T-T_IS,k);
            
            for i = 1:k
                params_i = params([1:4 4+i]);
                [h_vec_f,loglik_r2_f,u_vec_i]=Filter_New.GARCH_t_fcst(T,params_i,return_mat(:,i),T_IS);
                
                h_mat_fcst(:,i)= h_vec_f;
                loglik_mat_fcst(:,i)= loglik_r2_f;
                u_mat(:,i) = u_vec_i;
            end
            
            %loglik_vec_aggr_fcst = sum(loglik_mat_fcst,2);
            
        end
        
        
        
        
        
        
    end % end methods
end % end class
