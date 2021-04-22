classdef Filter_New_I % in Filter.m
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
        
                
        %% GAS filter with CT = Covariance Targeting
        function [ V, s,Flag] = GAS_CT( k, T, params, r2, RC )
                                    
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
            [s, f, V] = Filter_New_I.Store(k, T, V0);
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
                d1    =  1*w_j * r2(:,:,j) - 1*V_j;
                d2    = 1*nu_f1*((nu_f1 +nu_f2)/(nu_f2-k-1)*RC_j*inv(eye(k)+ constant*invV*RC_j) - V_j) ;
                
                nabla = (d1 + d2)/(1+nu_f1);        
                s(:,:,j) = nabla;
                
                if j<T
                    f(:,:,j+1) = (1-beta)*fbar + beta * f(:,:,j) + alpha * s(:,:,j);
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
            [s,~,V] = Filter_New_I.Store(k, T, V0);
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
        
        
        
        
        
    end % end methods
end % end class
