function  [sum_DM_mat,mean_loss] = Do_DM_HAC_test_given_matrix(loss_mat,model_ind_vec,mod_i)

mean_loss = mean(loss_mat);
nr_models = length(model_ind_vec);
model_ind_vec(mod_i) = [];
sum_DM_mat = zeros(2,nr_models-1);
T = size(loss_mat,1);
for m = 1:nr_models-1
     mod_j =  model_ind_vec(m);
 
     d_vec = loss_mat(:,mod_i) - loss_mat(:,mod_j) ;
     result        = olk(d_vec,ones(T,1));
     d_mean        = result.beta;
     HAC_var_mat   = Construct_HAC_se_NEW(result.resid,ones(T,1),[]);
     tstat_HAC_sigma_p  = d_mean/sqrt(HAC_var_mat);
     %sum_DM_vec = [sum_DM_vec tstat_HAC_sigma_p d_mean];

     sum_DM_mat(:,m) = [tstat_HAC_sigma_p d_mean];
end