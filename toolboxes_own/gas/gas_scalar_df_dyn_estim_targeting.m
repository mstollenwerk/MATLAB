function [eparam, tstats, logL, fit, fcst, optimoutput] = ...
	gas_scalar_df_dyn_estim_targeting( R, dist, x0, varargin )
%
%
% Michael Stollenwerk
% michael.stollenwerk@live.com
% 05.12.2022

[p,~,T] = size(R);
p_ = p*(p+1)/2;

options_ = [{'noperm'}, varargin];
[eparam.gasHARmean, tstats.gasHARmean, logL.gasHARmean, ...
 fit.gasHARmean, fcst.gasHARmean, optimoutput.gasHARmean] = ...
    gas_scalar_har_mean_estim_targeting( ...
        R, dist, x0, options_{:} ...
    );

%% Optimization
obj_fun = @(paramDFs) obj_fun_wrapper(eparam.gasHARmean, paramDFs, R, dist);
x0_dyn_df = [0.0001, .6, .2, .19];
% A = [ 0 1 1 1 ];
% b = 1;
lb = zeros(1,4);
ub = [inf .999 .999 .999];
if any(strcmp( dist, {'tWish','itWish','F','tRiesz','itRiesz2','FRiesz','iFRiesz2'} ))
    x0_dyn_df = [x0_dyn_df, x0_dyn_df];
%     A = [ 0 1 1 1 0 0 0 0;
%           0 0 0 0 0 1 1 1 ];
%     b = [.99;.99];
    lb = zeros(1,8);
    ub = [inf .999 .999 .999 inf .999 .999 .999];
end

eparam_dyn_df = ...
    my_fmincon(...
        @(param_dyn_df) obj_fun(param_dyn_df), ...
        x0_dyn_df', ...
        [],[],[],[], ...
        lb,ub,[], ...
        varargin{:} ...
    );
% [eparam_dyn_df] = ...
%     fminsearch(...
%         @(param_dyn_df) obj_fun(abs(param_dyn_df)), ...
%         x0_dyn_df', ...
%         optimset('Display','iter','PlotFcns',@optimplotx) ...
%     );

[ ~, ~, ~, ~, eparam_dyn_df ] = obj_fun( abs(eparam_dyn_df) );

[eparam.gas, tstats.gas, logL.gas, fit.gas, fcst.gas, optimoutput.gas] = ...
    gas_scalar_estim_targeting_( R, dist, eparam_dyn_df.all(p_+1:end), options_{:} );

end

function [ nLogL, logLcontr, dyn, S, param_out ] = ...
    obj_fun_wrapper(param_gasHARmean, paramDFs, R, dist) 
    
    if size(paramDFs,1)<size(paramDFs,2)
        paramDFs = paramDFs';
    end
        
    if isfield(param_gasHARmean,'n')
        persistence_n = sum(paramDFs(2:4));
        if persistence_n >= .999
            nLogL = inf;
            return
        end
        intrcpt_n = param_gasHARmean.n*(1-persistence_n);
        if isfield(param_gasHARmean,'nu')
            persistence_nu = sum(paramDFs(6:8));
            if persistence_nu >= 1
                nLogL = inf;
                return
            end
            intrcpt_nu = param_gasHARmean.nu*(1-persistence_nu);
            paramDFs_all = [intrcpt_n;paramDFs(1:4);intrcpt_nu;paramDFs(5:8)];
        else
            paramDFs_all = [intrcpt_n;paramDFs(1:4)];
        end
    elseif isfield(param_gasHARmean,'nu')
        persistence_nu = sum(paramDFs(2:4));
        if persistence_nu >= .999
            nLogL = inf;
            return
        end
        intrcpt_nu = param_gasHARmean.nu*(1-persistence_nu);
        paramDFs_all = [intrcpt_nu;paramDFs(1:4)];
    end
    
    s1Sig = param_gasHARmean.Sig.score1;
    s2Sig = param_gasHARmean.Sig.score2;
    gDSig = param_gasHARmean.Sig.garchD;
    gWSig = param_gasHARmean.Sig.garchW;
    gMSig = param_gasHARmean.Sig.garchM;
    meanSig = mean(R,3);
    vechcholIntrcptSig = vechchol( (1-gDSig-gWSig-gMSig)*meanSig );  

    param_all = [vechcholIntrcptSig; s1Sig; s2Sig; gDSig; gWSig; gMSig; paramDFs_all];
    
    [ nLogL, logLcontr, dyn, S, param_out ] = gas_scalar_likeRec_( ...
        param_all, ...
        R, ...
        dist ...
    );

end