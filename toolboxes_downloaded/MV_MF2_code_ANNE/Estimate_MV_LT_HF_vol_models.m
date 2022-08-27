function [Sum_mat,para_cell] = Estimate_MV_LT_HF_vol_models(return_mat,RC_mat,ind_1_step_out)

% JULY 2022
% Estimates various MV models including RV/RV data
% MV HEAVY (with N dist); 2 versions
% MF2 DCC (DCC without RCORR)
% MF2 DCC (DCC with LT RCORR)
% MF2 realized DCC 
% MF2 realized DCC with LT RCORR 

% INPUT: 
%  T x k return_mat
%  vec of k * k RCOV mat
%  ind_1_step out = 1 if you want to discard the 1-step estimation of the
%  MF2-DCC model (with returns only), otherwise 0.

options = optimset('fmincon');
options = optimset(options, 'Display', 'iter-detailed');
options = optimset(options, 'Diagnostics', 'off');
options = optimset(options, 'Algorithm', 'sqp');
options.MaxIter = 30000;
options.MaxFunEvals = 30000;
import LogLik_MV.*

[T,k] = size(return_mat);
T_IS = T;
T_w = 63;

theta_N = [0.05 0.85]';
Lb_N = [zeros(2,1) + eps];
Rb_N = [1-eps 1-eps]';

% select  data, compute RV and RCORR out of
return_mat_i  = return_mat(1:T_IS,:);
RC_mat_i = RC_mat(:,:,1:T_IS);
[RCORR_mat_i,RV_mat_i] = Admin.Compute_implicit_Corr_matrices(RC_mat_i);

% estimate 2 heavy models
% f_HEAVY_N =  @(params)LogLik_MV.LogLik_MV_HEAVY_N(k,T_IS,T_w,params,return_mat_i,RC_mat_i,1);
% [theta_HEAVY_N_hat,fval_HEAVY_N] = fmincon(f_HEAVY_N,theta_N,[],[],[],[],Lb_N,Rb_N,[],options);
% f_HEAVY_N2 =  @(params)LogLik_MV.LogLik_MV_HEAVY_LT_spec_N(k,T_IS,T_w,params,return_mat_i,RC_mat_i);
% [theta_HEAVY_N2_hat,fval_HEAVY_N2] = fmincon(f_HEAVY_N2,theta_N,[],[],[],[],Lb_N,Rb_N,[],options);
% 
% para_cell{1} = theta_HEAVY_N_hat;
% para_cell{2} = theta_HEAVY_N2_hat;
% 
% [~,logl_vec_HEAVY_N,V_HEAVY_N] = f_HEAVY_N(theta_HEAVY_N_hat);
% [~,logl_vec_HEAVY_N2,V_HEAVY_N2] = f_HEAVY_N2(theta_HEAVY_N2_hat);

%----------------------------------------------------------------------------------------------



% OLD check autocorrelatoins of the ''V_tilde'' of the MF2 model
%Compute_cyclical_patterns_HEAVY(V_HEAVY_N2,return_mat_i,RC_mat_i,0);


% OLD MF2 with CCC
% theta_MF_N_diag_CCC = [0.2*ones(k,1); theta_HEAVY_N2_hat(1:2); 0.02; 0.8];
% Lb_MF_N = [zeros(k+4,1) + eps;];
% Rb_MF_N = [30*ones(k,1); ones(4,1)-eps];
%
% A = [zeros(1,k) 1 1 0 0; zeros(1,k) 0 0 1 1];
% b = ones(2,1)-eps;
%
% f_MF_HEAVY_N_diag_CCC =  @(params)LogLik_MV.LogLik_MV_MF2_HEAVY_diag_CCC_N(k,T,T_w,params,return_mat_i,RC_mat_i);
% [theta_HEAVY_MF_N_CCC_hat,fval_HEAVY_MF_N_CCC] = fmincon(f_MF_HEAVY_N_diag_CCC,theta_MF_N_diag_CCC,A,b,[],[],Lb_MF_N,Rb_MF_N,[],options);
% [~,~,H_TOT_mat_CCC,H_ST_mat_CCC, Z_t_mat_CCC,V_t_mat_CCC] = f_MF_HEAVY_N_diag_CCC(theta_HEAVY_MF_N_CCC_hat);

% theta_DCC = [0.02 0.8]';
% Lb_DCC = zeros(2,1)+eps;
% Rb_DCC = ones(2,1)-eps;
% A_DCC =[1 1]; b_DCC = 1-eps;
%
% f_DCC = @(params)LogLik_MV.LogLik_cDCC_PSEUDO_N(k,T,params,u_mat);
% [theta_DCC_hat,fval_DCC] = fmincon(f_DCC,theta_DCC,A_DCC,b_DCC,[],[],Lb_DCC,Rb_DCC,[],options);
%
% [~,Rt_mat_DCC] = f_DCC(theta_DCC_hat);

%------------------------------------------------------------------------------------------------
% MF2 DCC (1 step estimation)
% if ind_1_step_out==0
%     
%     theta_MF_N_diag_DCC = [0.2*ones(k,1); theta_HEAVY_N2_hat(1:2);0.02; 0.8;0.01; 0.9];
%     Lb_MF_N = [zeros(k+6,1) + eps;];
%     Rb_MF_N = [30*ones(k,1); ones(6,1)-eps];
%     A = [zeros(1,k) 1 1 0 0 0 0; zeros(1,k) 0 0 1 1 0 0; zeros(1,k+4) 1 1];
%     b = ones(3,1)-eps;
%     
%     T_w = 63;
%     f_MF_HEAVY_N_diag_DCC =  @(params)LogLik_MV.LogLik_MV_MF2_HEAVY_diag_DCC_N(k,T_IS,T_w,params,return_mat_i,RC_mat_i);
%     [theta_HEAVY_MF_N_DCC_hat,fval_HEAVY_MF_N_DCC] = fmincon(f_MF_HEAVY_N_diag_DCC,theta_MF_N_diag_DCC,A,b,[],[],Lb_MF_N,Rb_MF_N,[],options);
%     [~,~,H_TOT_mat_DCC,H_ST_mat_DCC,Z_t_mat_DCC,R_t_mat_DCC,V_t_mat_DCC] = f_MF_HEAVY_N_diag_DCC(theta_HEAVY_MF_N_DCC_hat);
%     
%     para_cell{3} = theta_HEAVY_MF_N_DCC_hat;
% end


% 2 step estimation; 
%STEP 1 VOL estimation
% Para_start_2 = [0.02 0.05 0.90 0.02 0.02 0.96]';
% Lb_2 = [-0.5; zeros(5,1)+eps];
% Rb_2= [0.4 ones(1,2)-eps 5 ones(1,2)-eps]';
% 
% % estimate volatilities asset by asset, to get a clue if it works
% for j =1 : k
%      f_GARCH_MF_N = @(params)LogLik_UV.Loglike_GARCH_RV_MF_N(T_IS,params,return_mat_i(:,j),RV_mat_i(:,j),T_w);
%      [theta_2,LLF_2,Flag_2] = fmincon(f_GARCH_MF_N,Para_start_2,[],[],[],[],Lb_2,Rb_2,[],options);                    
%      theta_2_mat(j,:) = [theta_2' LLF_2];  
% end
    

% estimate vol models by pooling alfa and beta.
theta_MF_HEAVY_N_VOL = [0.05*ones(k,1); 0.2; 0.6; 0.2; 0.6];
Lb_MF_N_VOL = [zeros(k+3,1) + eps; 0.5];
Rb_MF_N_VOL = [0.5*ones(k,1); ones(2,1)-eps; 0.4; 1-eps];
A = [zeros(1,k) 1 1 0 0; zeros(1,k) 0 0 1 1];
b = ones(2,1)-eps;

f_MF_HEAVY_N_VOL =  @(params)LogLik_MV.LogLik_MV_MF2_HEAVY_VOL_PART(k,T_IS,T_w,params,return_mat_i,RV_mat_i);
[theta_HEAVY_MF_N_VOL_hat,fval_HEAVY_MF_N_VOL] = fmincon(f_MF_HEAVY_N_VOL ,theta_MF_HEAVY_N_VOL,A,b,[],[],Lb_MF_N_VOL,Rb_MF_N_VOL,[],options);
[~,h_it_mat,tau_mat,h_it_mat_ST] = f_MF_HEAVY_N_VOL(theta_HEAVY_MF_N_VOL_hat);
u_it_mat = return_mat_i./sqrt(h_it_mat);


% STEP 2 Several DCC models
theta_DCC = [0.02 0.8]';
theta_DCC_RC = [0.15 0.6]';
Lb_DCC = zeros(2,1)+eps;
Rb_DCC = ones(2,1)-eps;
A_DCC =[1 1]; b_DCC = 1-eps;

% 2.1 DCC with returns only
% f_DCC = @(params)LogLik_MV.LogLik_cDCC_PSEUDO_N(k,T_IS,T_w,params,u_it_mat);
% [theta_DCC_hat,fval_DCC_2s] = fmincon(f_DCC,theta_DCC,A_DCC,b_DCC,[],[],Lb_DCC,Rb_DCC,[],options);
% [~,Rt_mat_DCC] = f_DCC(theta_DCC_hat);
% para_cell{4} = [theta_HEAVY_MF_N_VOL_hat; theta_DCC_hat];

% 2.2 DCC with RC (in long term)
% f_DCC_RC = @(params)LogLik_MV.LogLik_cDCC_PSEUDO_N_RC(k,T_IS,T_w,params,u_it_mat,RCORR_mat_i);
% [theta_DCC_hat_RC,fval_DCC_RC] = fmincon(f_DCC_RC,theta_DCC_RC,A_DCC,b_DCC,[],[],Lb_DCC,Rb_DCC,[],options);
% [~,Rt_mat_DCC_RC,Qbar_DCC_t] = f_DCC_RC(theta_DCC_hat_RC);
% para_cell{5} = [theta_HEAVY_MF_N_VOL_hat; theta_DCC_hat_RC];

% 2.3: realzied DCC(geen Q stuff!)
f_DCC_RC2 = @(params)LogLik_MV.LogLik_REALIZED_DCC_PSEUDO_N(k,T_IS,T_w,params,u_it_mat,RCORR_mat_i);
[theta_DCC_hat_RC2,fval_DCC_RC2] = fmincon(f_DCC_RC2,theta_DCC_RC,A_DCC,b_DCC,[],[],Lb_DCC,Rb_DCC,[],options);
[~,Rt_mat_DCC_RC2] = f_DCC_RC2(theta_DCC_hat_RC);
para_cell{6} = [theta_HEAVY_MF_N_VOL_hat; theta_DCC_hat_RC2];

% 2.4: Realized DCC plus Longterm
f_DCC_RC2_LT = @(params)LogLik_MV.LogLik_REALIZED_DCC_PSEUDO_N2(k,T_IS,T_w,params,u_it_mat,RCORR_mat_i);
[theta_DCC_hat_RC2_LT,fval_DCC_RC_LT] = fmincon(f_DCC_RC2_LT,theta_DCC_hat_RC2,A_DCC,b_DCC,[],[],Lb_DCC,Rb_DCC,[],options);
[~,Rt_mat_DCC_RC2_LT,Qbar_DCC_RC2_LT] = f_DCC_RC2_LT(theta_DCC_hat_RC2_LT);
para_cell{7} = [theta_HEAVY_MF_N_VOL_hat; theta_DCC_hat_RC2_LT];

% compute final fitted covariance matrices.
% H_TOT_mat_DCC_2st= Admin.Compute_H_mat_given_corr_and_vol(Rt_mat_DCC,h_it_mat);
% H_TOT_mat_DCC_RC= Admin.Compute_H_mat_given_corr_and_vol(Rt_mat_DCC_RC,h_it_mat);
% H_TOT_mat_DCC_RC2= Admin.Compute_H_mat_given_corr_and_vol(Rt_mat_DCC_RC2,h_it_mat);
% H_TOT_mat_DCC_RC_LT= Admin.Compute_H_mat_given_corr_and_vol(Rt_mat_DCC_RC2_LT,h_it_mat);

% compute loglik
% fval_DCC_2s = fval_HEAVY_MF_N_VOL + fval_DCC_2s;
% fval_DCC_RC = fval_HEAVY_MF_N_VOL + fval_DCC_RC;
fval_DCC_RC2 = fval_HEAVY_MF_N_VOL + fval_DCC_RC2;
% fval_DCC_RC_LT = fval_HEAVY_MF_N_VOL + fval_DCC_RC_LT;

% put fitted cov matrices into a cell
% H_mat_cell_HAT{1}  = V_HEAVY_N;
% H_mat_cell_HAT{2}  = V_HEAVY_N2;
% if ind_1_step_out ==1
%     H_mat_cell_HAT{3}  = H_TOT_mat_DCC_2st;
%     H_mat_cell_HAT{4}  = H_TOT_mat_DCC_RC;
%     H_mat_cell_HAT{5}  = H_TOT_mat_DCC_RC2;
%     H_mat_cell_HAT{6}  = H_TOT_mat_DCC_RC_LT;
%     loglik_vec = -[fval_HEAVY_N fval_HEAVY_N2 fval_DCC_2s fval_DCC_RC fval_DCC_RC2 fval_DCC_RC_LT];
% else
%     loglik_vec = -[fval_HEAVY_N fval_HEAVY_N2 fval_HEAVY_MF_N_DCC fval_DCC_2s fval_DCC_RC fval_DCC_RC2 fval_DCC_RC_LT];
%     H_mat_cell_HAT{3}  = H_TOT_mat_DCC;
%     H_mat_cell_HAT{4}  = H_TOT_mat_DCC_2st;
%     H_mat_cell_HAT{5}  = H_TOT_mat_DCC_RC;
%     H_mat_cell_HAT{6}  = H_TOT_mat_DCC_RC2;
%     H_mat_cell_HAT{7}  = H_TOT_mat_DCC_RC_LT;
% end

% compute loss functions and play the GMV portfolio
% [mean_loss,loss_mat_FR,loss_mat_QLIK] = Compute_loss_functions_matrices(RC_mat_i,H_mat_cell_HAT);
% [sum_DM_vec,mean_sigma,sigma_real_port_vec ] = Do_GMV_exercise(H_mat_cell_HAT,RC_mat_i,return_mat);
% Sum_mat = [loglik_vec; mean_loss; mean_sigma; sigma_real_port_vec];


