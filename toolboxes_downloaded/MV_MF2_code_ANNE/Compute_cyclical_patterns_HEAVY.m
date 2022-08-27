function Compute_cyclical_patterns_HEAVY(V_HEAVY,return_mat,RCOV_mat,ind_HF_only);


T = size(RCOV_mat,3);
k = size(RCOV_mat,1);
u_tilde_mat = zeros(k,k,T);
u_tilde_mat_HF = zeros(k,k,T);
for i = 1:T
   
    V_i = V_HEAVY(:,:,i);
    chol_V = chol(V_i)';
    inv_chol_V = inv(chol_V);
    if ind_HF_only==1
        u_tilde_mat_HF(:,:,i) = inv_chol_V * RCOV_mat(:,:,i) * inv_chol_V';
    else
        
        return_outer_i = return_mat(i,:)'*return_mat(i,:);
        u_tilde_i = inv_chol_V * return_outer_i * inv_chol_V';
        u_tilde_mat(:,:,i) = u_tilde_i;
        u_tilde_mat_HF(:,:,i) = inv_chol_V * RCOV_mat(:,:,i) * inv_chol_V';
    end
    
    
end

apie =1;

T_w = 63;

 nr_off  = k*(k-1)/2;
 k_plot = ceil(sqrt(nr_off));
 k_plot_diag = ceil(sqrt(k));

if ind_HF_only==1
    
    
    teller = 1;
    for t = T_w : T_w : T
        u_cycle_t_HF = mean(u_tilde_mat_HF(:,:,t-T_w+1:t),3);
        u_diag_mat_HF(teller,:) = diag(u_cycle_t_HF);
        u_off_diag_mat_HF(teller,:) = Admin.Sym2UT(k,u_cycle_t_HF);
        teller= teller + 1;
    end
    figure;
    
    for i = 1:k
        subplot(k_plot_diag ,k_plot_diag ,i)
        autocorr(u_diag_mat_HF(:,i));
    end
    
    figure;
   
    for i = 1:nr_off
        subplot(k_plot,k_plot,i)
        autocorr(u_off_diag_mat_HF(:,i));
    end
    
else
    teller = 1;
    for t = T_w : T_w : T
        u_cycle_t = mean(u_tilde_mat(:,:,t-T_w+1:t),3);
        u_diag_mat(teller,:) = diag(u_cycle_t);
        u_off_diag_mat(teller,:) = Admin.Sym2UT(k,u_cycle_t);
        teller= teller + 1;
    end
    figure;
    for i = 1:k
       subplot(k_plot_diag ,k_plot_diag ,i)
        autocorr(u_diag_mat(:,i));
    end
    
    figure;
    
    for i = 1:nr_off
        subplot(k_plot,k_plot,i)
        autocorr(u_off_diag_mat(:,i));
    end
    
    teller = 1;
    for t = T_w : T_w : T
        u_cycle_t_HF = mean(u_tilde_mat_HF(:,:,t-T_w+1:t),3);
        u_diag_mat_HF(teller,:) = diag(u_cycle_t_HF);
        u_off_diag_mat_HF(teller,:) = Admin.Sym2UT(k,u_cycle_t_HF);
        teller= teller + 1;
    end
    
    figure;
    for i = 1:k
       subplot(k_plot_diag ,k_plot_diag ,i)
        autocorr(u_diag_mat_HF(:,i));
    end
    
    figure;    
    for i = 1:nr_off
        subplot(k_plot,k_plot,i)
        autocorr(u_off_diag_mat_HF(:,i));
    end
end












