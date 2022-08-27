
classdef Admin % in Admin.m
    methods (Static = true)
        
        
        function H_TOT_mat= Compute_H_mat_given_corr_and_vol(Rt_mat,h_it_mat)
            
            % maak 3wegmatrices voor schatten
            [T,k] = size(h_it_mat);
            H_TOT_mat = zeros(k,k,T);
            sqrt_h_it_mat = sqrt(h_it_mat);
            
            for j = 1:T
                H_TOT_mat(:,:,j) = diag(sqrt_h_it_mat(j,:)) * Rt_mat(:,:,j) * diag(sqrt_h_it_mat(j,:));
            end
            
        end
        
        
        function RCOV_mat= Compute_RK_mat(RM_mat_vec,k)
            
            % maak 3wegmatrices voor schatten
            T_end = size(RM_mat_vec,1);
            RCOV_mat = zeros(k,k,T_end);
            
            for j = 1:T_end
                RCOV_mat(:,:,j) = Admin.Vec2Sym(k, RM_mat_vec(j,:));
            end
            
        end
        
        function RCOV_mat= Compute_RK_mat_given_vech(RM_mat_vech,k)
            
            % maak 3wegmatrices voor schatten
            T_end = size(RM_mat_vech,1);
            RCOV_mat = zeros(k,k,T_end);
            
            for j = 1:T_end
                RCOV_mat(:,:,j) = Admin.Vech2Sym(k, RM_mat_vech(j,:));
            end
            
        end
        
        
        function S_d_sp= Compute_S_d_matrix(asset_group_vec)
            nr_groups = max(asset_group_vec);
            n = length(asset_group_vec);
            
            S_d = zeros(n^2,nr_groups);
            
            for i = 1:n
                S_d(i + (i-1)*n,asset_group_vec(i)) = 1;
            end
            S_d_sp = sparse(S_d);
            
        end
        
        
        
        function  [R,Var_mat] = Compute_implicit_Corr_matrices(V)
            
            T = size(V,3);
            k = size(V,2);
            
            R = zeros(k,k,T);
            Var_mat = zeros(T,k);
            for i = 1:T
                V_i = V(:,:,i);
                [sigma_vec, R_i] = cov2corr(V_i);
                Var_mat(i,:) = sigma_vec.^2;
                R(:,:,i) = R_i;
                
            end
            
        end
        
        function mu_vec = Compute_1st_moment_Inv_Riesz_I(nu_vec);
            
            k = length(nu_vec);
            mu_vec = zeros(k,1);
            for j = 1:k-1
                i_vec = j+1:k;
                mu_vec(j,1) = 1/(nu_vec(j)-(j+1)) * prod( (nu_vec(j+1:k)' - i_vec)./(nu_vec(j+1:k)'- (i_vec + 1)));
            end
            mu_vec(end,1) = 1/(nu_vec(end)-(k+1));
            %E_X = diag(mu_vec);
            
        end
        
        function mu_vec = Compute_1st_moment_Inv_Riesz_II(nu_vec);
            
            k = length(nu_vec);
            mu_vec = zeros(k,1);
            mu_vec(1,1) = 1/(nu_vec(1)-(k+1));
            
            for i = 2:k
                j_vec_0 = 1:i-1;
                j_vec = k - j_vec_0 + 1;
                mu_vec(i,1) = 1/(nu_vec(i)-(k-i+2)) * prod( (nu_vec(1:i-1)' - j_vec)./(nu_vec(1:i-1)'- (j_vec + 1)));
            end
        end
        
        
        function mu_vec = Compute_1st_moment_FRiesz_I(nu1_vec,nu2_vec);
            % NOTE: this is the moments of an FRIESZ-I obtained by Theorem
            % 12 of the paper!
            % X - Riesz(I), LL' = Y = Riesz(I)
            % Z = Linv X Linv'
            
            k = length(nu1_vec);
            mu_vec = zeros(k,1);
            mu_vec(1) = nu1_vec(1)/(nu2_vec(1)-2);
            for j = 2:k
                mu_vec(j,1) = 1/(nu2_vec(j)-(j+1)) *  (nu1_vec(j) + sum(mu_vec(1:j-1)));
            end
            
        end
        
        
        function mu_vec = Compute_1st_moment_FRiesz_II(nu1_vec,nu2_vec);
            % NOTE: this is the moments of an FRIESZ-II obtained by Theorem
            % X of the paper!
            % X - Riesz(II), UU' = Y = Riesz(I)
            % Z = Uinv X Uinv'
            
            % NEW
            k = length(nu1_vec);
            mu_vec = zeros(k,1);
            mu_vec(k) = nu1_vec(k)/(nu2_vec(k)-2);
            for i = k:-1:2
                mu_vec(i-1,1) = 1/(nu2_vec(i-1)-(2+k-(i-1))) *  (nu1_vec(i-1) + sum(mu_vec(k:-1:i)));
            end
            
        end
        
        %%
        function [uLT] = UnitLT(k)
            % vult lower trianguler met enen.
            uLT= zeros(k,k);
            
            for j=1:k
                for i=j+1:k
                    uLT(i,j) =1;
                end
            end
            
        end
        
        
        %%
        function [Cov] = CovOneFactor(k, sigma, rho)
            % bouwt cov matrix met 1 sigma en equicorrelation
            %Cov = zeros(k, k);
            Vol = zeros(k, k);
            Corr= eye(k, k);
            
            for i=1:k
                Vol(i,i)= sigma;
            end
            
            
            for i=1:k
                for j=(i+1):k
                    Corr(i,j)= rho;
                    Corr(j,i)= Corr(i,j);
                end
            end
            
            Cov= Vol*Corr*Vol;
        end
        
        
        function RCOV_vech_tot = Convert_vec_to_vech(RCOV_mat,k)
            
            date_vec = RCOV_mat(:,1);
            RCOV_mat(:,1) = [];
            
            T = size(RCOV_mat,1);
            RCOV_vech = zeros(T,k*(k+1)/2);
            
            for t = 1:T
                A = Admin.Vec2Sym(k, RCOV_mat(t,:));
                [vechA] = Admin.Sym2Vech(k, A);
                
                RCOV_vech(t,:) = vechA;
                
            end
            
            RCOV_vech_tot = [date_vec RCOV_vech];
        end
        
        function RCOV_mat= Convert_vec_to_cov(RCOV_mat_vec,k)                        
                        
            T = size(RCOV_mat_vec,1);
            RCOV_mat = zeros(k,k,T);
            
            for t = 1:T
                A = Admin.Vec2Sym(k, RCOV_mat_vec(t,:));
                RCOV_mat(:,:,t) = A;                
            end
        end
        
        function RCOV_mat= Convert_vech_to_cov(RCOV_mat_vech,k)                                    
            
            T = size(RCOV_mat_vech,1);
            RCOV_mat = zeros(k,k,T);
            
            for t = 1:T
                A = Admin.Vech2Sym(k, RCOV_mat_vech(t,:));
                RCOV_mat(:,:,t) = A;
                
            end
        end
        
        
        %%
        function [ A ] = Vech2Sym(k, vechA )
            % this function puts vech of A into A
            
            A= zeros(k,k);
            l=1;
            for j=1:k
                for i=j:k
                    A(i,j)= vechA(l);
                    A(j,i)= A(i,j);
                    l=l+1;
                end
            end
            
        end
        
        function [ A ] = Vec2Sym(k, vecA )
            % this function puts vech of A into A
            
            A= zeros(k,k);
            l=1;
            for j=1:k
                for i=1:k
                    A(i,j)= vecA(l);
                    l=l+1;
                end
            end
            
        end
        
        
        function [ A ] = VecCorr2Sym(k, vecA )
            % this function puts vec of LT mat of A of into A
            
            A= zeros(k,k);
            l=1;
            for j=1:k
                for i=j+1:k
                    A(i,j)= vecA(l);
                    l=l+1;
                end
            end
            
        end
        
        
        
        %%
        function [ vechA ] = Sym2Vech(k, A )
            % this function puts A into vech of A
            
            vechA= zeros(k*(k+1)/2, 1);
            l=1;
            for j=1:k
                for i=j:k
                    vechA(l) = A(i,j);
                    l=l+1;
                end
            end
            
        end
        
        %%
        function [ vechltA ] = LowerTr2Vech(k, ltA )
            % this function puts lower triangular A into  vech lt A
            
            vechltA= zeros(k*(k+1)/2, 1);
            l=1;
            for j=1:k
                for i=j:k
                    vechltA(l)= ltA(i,j);
                    l=l+1;
                end
            end
            
        end
        
        function [UtA] =Sym2UT(k,A);
            
            UtA = zeros(k*(k-1)/2,1);
            l = 1;
            for j = 1:k
                UtA(l:l+k-j-1,1) = A(j,j+1:end);
                l =l+k-j;
            end
            
        end
        
        function [A] =Ut2Sym(k,UtA);
            % with ones on the diagonal
            A = eye(k);
            l = 1;
            for j = 1:k
                A(j,j+1:end) = UtA(l:l+k-j-1,1);
                A(j+1:end,j) =  A(j,j+1:end);
                l =l+k-j;
            end
            
        end
        
        
        
        
        %%
        function [ ltA ] = Vech2LowerTr(k, vechltA )
            % this function puts vector of vech of lt A into lt A
            
            ltA= zeros(k,k);
            l=1;
            for j=1:k
                for i=j:k
                    ltA(i,j)= vechltA(l);
                    l=l+1;
                end
            end
            
        end
        
        %%
        function [K_sp] = commutation(m,  n)
            
            I= eye(m*n, m*n);
            K= zeros(m*n, m*n);
            
            for j=1:n
                for i=1:m
                    K((i-1)*n + j,:)= I((j-1)*m + i, :);
                end
            end
            
            K_sp = sparse(K);
            
        end
        
        
        %%
        function [ D ] = duplication(k)
            % creates duplication matrix
            
            I= eye(k*(k+1)/2);
            D= zeros( k^2, k*(k+1)/2 );
            
            for j=1:k*(k+1)/2
                for i=1:k
                    if (i>=j)
                        D((j-1)*k + i ,:)= I((j-1)*(k - j/2) +i ,:);
                        D((i-1)*k + j ,:)= I((j-1)*(k - j/2) +i ,:);
                    end
                end
            end
            
        end
        
        %%
        function [ L ] = elimination(k)
            % creates elimination matrix
            
            D= Admin.duplication(k);
            L= inv(D' * D) * D';
            
        end
        
        
        %%
        function [DL] = duplication_gen(k)
            
            D = Admin.duplication(k);
            LT= ones(k, k);
            
            for i=1:k
                for j=i+1:k
                    LT(i,j)= 0;
                end
            end
            
            LT= reshape(LT, k^2, 1) * ones(1, k*(k+1)/2 );
            DL= LT .* D;
            
        end
        
        %%
        function [LL] = elimination_gen(k)
            
            DL= Admin.duplication_gen(k);
            LL= DL';
            
        end
        
        
        %%
        function [sqroot_matA] = sqroot_mat(A)
            
            [V, D]= eig(A);
            sqroot_matA= V * D.^(0.5) * inv(V);
            
        end
        
        %%
        function [sqroot_matA] = sqroot_symmat(A)
            
            [V, D]= eig(A);
            sqroot_matA= V * D.^(0.5) * V';
            
        end
        
        %%
        function [sqroot_matA] = sqroot_invsymmat(A)
            
            [V, D]= eig(A); %TODO: eig robust for eigen > 0
            sqroot_matA= V* (diag(diag(D).^(-1/2))) * V';
            
        end
    end % end methods
end % end class