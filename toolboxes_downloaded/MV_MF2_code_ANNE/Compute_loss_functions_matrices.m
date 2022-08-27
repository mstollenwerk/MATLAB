function [mean_loss,loss_mat_FR,loss_mat_QLIK] = Compute_loss_functions_matrices(R_mat,R_mat_cell_HAT);

T = size(R_mat,3);
k = size(R_mat,1);

nr_models = length(R_mat_cell_HAT);

loss_mat_FR = zeros(T,nr_models);
loss_mat_QLIK = zeros(T,nr_models);
I_k = eye(k);

for i =1:T
    R_mat_i = R_mat(:,:,i);
    
    for j = 1: nr_models
        
        R_mat_HAT_i = R_mat_cell_HAT{j}(:,:,i);
        R_mat_HAT_inv_i = R_mat_HAT_i\I_k;
        
        diff_i = R_mat_i - R_mat_HAT_i;
        vech_i = Admin.Sym2Vech(k, diff_i);
        loss_mat_FR(i,j) = vech_i'*vech_i;
        loss_mat_QLIK(i,j) = log(det(R_mat_HAT_i)) + trace(R_mat_HAT_inv_i*R_mat_i);
        
    end
    
end
mean_loss = [mean(loss_mat_FR); mean(loss_mat_QLIK)];