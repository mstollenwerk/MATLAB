function HAC_var_mat = Construct_HAC_se_NEW(resids_vec,X_mat,L);

[T k]       = size(X_mat);

if isempty(L)
    L           = ceil(T^(1/5));
end

%construct V_hat
V_cross  = zeros(k,k);
X_mat_new = ones(k,1)*resids_vec'.*X_mat';
V_hat = X_mat_new*X_mat_new';

for i  = 1:L
    
    w_i = 1- i/L;
    %w_i = 1- i/(L+1);
    % NOTE: 1- i/L is bartlett kernel, used by matlab ifyou use 'hac'
    
    %w_resid_vec_i = w_i*resid_vec(i+1:end).*resid_vec(1:end-i);    
    X_eps_i = (ones(k,1)*resids_vec(i+1:end)'.*X_mat(i+1:end,:)')*(X_mat(1:end-i,:).*(resids_vec(1:end-i)*ones(1,k)));           
    V_cross   = V_cross + w_i *(X_eps_i + X_eps_i');
end
%V_cross/T
V_tot_mat = (V_hat + V_cross)/T;

inv_XX    = inv((1/T)*X_mat'*X_mat);
HAC_var_mat = (1/T) * inv_XX*V_tot_mat*inv_XX;
