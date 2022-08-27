function [sum_DM_vec,mean_sigma,sigma_real_port_vec ] = Do_GMV_exercise(H_mat_cell,Sigma,return_mat)

nr_models = length(H_mat_cell);
T = size(H_mat_cell{1},3);
k = size(H_mat_cell{1},1);
ones_vec = ones(k,1);
%port_return_mat = zeros(T,nr_models);
ex_post_sigma_mat = zeros(T,nr_models);
ex_post_var_mat = zeros(T,nr_models);

for j = 1:nr_models
    w_cell{j} =  zeros(T,k);
end

for i = 1:T
    
    for j = 1:nr_models
        H_j = H_mat_cell{j};
                
        if size(H_j,2)>k
            H_ji = Admin.Vech2Sym(k,H_j(i,:));
        else
            H_ji = H_j(:,:,i);
        end
        
        
        inv_H_ji = inv(H_ji);
        w_ji = (inv_H_ji * ones_vec)/(ones_vec'*inv_H_ji*ones_vec);
        
        w_cell{j}(i,:) = w_ji;
        port_return_mat(i,j) = return_mat(i,:) *  w_ji;
        ex_post_sigma_mat(i,j) = sqrt(w_ji'*Sigma(:,:,i) * w_ji);
        %ex_post_var_mat(i,j) = w_ji'*Sigma(:,:,i) * w_ji;
    end
    
end

mean_ep_sigma_vec = mean(ex_post_sigma_mat);
sigma_real_port_vec  = std(port_return_mat);

nr_model_vec = 1:nr_models;
mod_i = 2;

sum_DM_vec = Do_DM_HAC_test_given_matrix(ex_post_sigma_mat,nr_model_vec ,mod_i);
mean_sigma = [mean_ep_sigma_vec];
