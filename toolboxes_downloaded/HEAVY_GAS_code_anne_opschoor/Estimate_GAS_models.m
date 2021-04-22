function [Para_cell,Stat_mat] = Estimate_GAS_models(return_mat,RM_mat_vec,ind_model_vec);

%ind_model_vec
%1: GAS
%2: GAS HAR

Para_cell = cell(2,1);
Stat_mat = zeros(2,3);
%LLF LLF_1 LLF_2 


% define the class with all admin utils, DGP and filters
import Admin.*
import Filter_New_I.*

options = optimset('fmincon');
options = optimset(options, 'Display', 'off');
options = optimset(options, 'Diagnostics', 'off');
options = optimset(options, 'Algorithm', 'sqp');

[T,k] = size(return_mat);
nr_models = length(ind_model_vec);

% our model
params0  = [0.08 0.9 15 25 25]';
Lb0 = [0.001 0.0001 5.01 k+1+0.01 k+1+0.01]';
Rb0 = [3 0.9990 1000 1000 1000]';

params_H  = [0.08 0.3 0.3 0.3 15 25 25]';
Lb_H = [0.001 0.0001*ones(1,3) 5.01 k+1+0.01 k+1+0.01]';
Rb_H = [3 0.999*ones(1,3) 1000 1000 1000]';
A_H  = [0 1 1 1 0 0 0];
b = 1 - eps;

% maak 3wegmatrices voor schatten
RK_mat = zeros(k,k,T);
r2_mat = zeros(k,k,T);
for j = 1:T
    RK_mat(:,:,j) = Admin.Vec2Sym(k, RM_mat_vec(j,:));    
    r2_mat(:,:,j) = return_mat(j,:)'*return_mat(j,:);
end

% Schat modellen (insample)
f_GAS_tF     = @(params)LogLik_I.LogLik_GAS_tF_CT(k, T, params, r2_mat, RK_mat);
f_GAS_tF_HAR = @(params)LogLik_I.LogLik_GAS_tF_HAR_CT(k, T, params, r2_mat, RK_mat);

for i = 1 :nr_models
    
    ind_model_vec_i = ind_model_vec(i);
    
    
    switch ind_model_vec_i
        
        case 1
            
            [theta_JLO, LLF_JLO] = fmincon(f_GAS_tF, params0,[],[],[],[],Lb0,Rb0,[],options);
            %Hessian = hessian_2sided(f_GAS_tF,theta_JLO);
            %se_vec = sqrt(diag(inv(Hessian)));
            %Para_cell{1} = [theta_JLO se_vec theta_JLO./se_vec];
            
            [Vfil_JLO, ffil_JLO]    = Filter_New_I.GAS_CT(k, T, theta_JLO, r2_mat, RK_mat);            
            [~, LLF_T, LLF_F] = f_GAS_tF(theta_JLO);
            
            Stat_mat(1,:) = [-LLF_JLO LLF_T LLF_F];
            
            
        case 2
            [theta_JLO_H, LLF_JLO_H] = fmincon(f_GAS_tF_HAR, params_H,A_H,b,[],[],Lb_H,Rb_H,[],options);
            
            %Hessian      = hessian_2sided(f_GAS_tF_HAR,theta_JLO_H);
            %se_vec       = sqrt(diag(inv(Hessian)));
            %Para_cell{2} = [theta_JLO_H se_vec theta_JLO_H./se_vec];
            
            Vfil_JLO_H               = Filter_New_I.GAS_HAR_CT(k, T, theta_JLO_H,r2_mat,RK_mat);            
            [~, LLF_T_H, LLF_F_H]    = f_GAS_tF_HAR(theta_JLO_H);
            
            Stat_mat(2,:) = [-LLF_JLO_H LLF_T_H LLF_F_H];
            
            
            
            
            
    end
    
end










