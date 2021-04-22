function [X,optimoutput] = my_fmincon_robust(fun,X,LB,UB,subgroup_idx,sensitivity_,varargin)
%MY_fminunc_robust Uncontstraint optimization over all parameters first and
%then over all parameter subgroups secified in subgroup_idx. Keep doing
%this until function value is improved leass than sensitivity_.

subgroup_idx = [0; subgroup_idx];
improvement_ = inf;
nlogL = inf;
ii = 1;
if isempty(LB)
    LB = inf(size(X,1),1);
end
if isempty(UB)
    UB = inf(size(X,1),1);
end
while improvement_ > sensitivity_

    [X_new,optimoutput{ii,1}] = my_fmincon(fun,X,[],[],[],[],LB,UB,[],varargin{:});
    
    jj = 1;
    for kk = 1:length(subgroup_idx)-1
        idx_ = subgroup_idx(kk)+1:subgroup_idx(kk+1);
        if ~isempty(idx_)
            indices{jj} = idx_;
            jj = jj+1;
        end
    end
    
    for jj = 1:numel(indices)
         
        fun_jj = @(x) wrapper_function(x,fun,X_new,indices{jj});

        [X_new(indices{jj}),optimoutput{ii,jj+1}] = ...
                my_fmincon( ...
                    fun_jj, ...
                    X_new(indices{jj}),[],[],[],[], ...
                    LB(indices{jj}), ...
                    UB(indices{jj}),[], ...
                    varargin{:} ...
                );

    end

    improvement_ = nlogL - optimoutput{ii,end}.history.fval(end);
    if improvement_ > sensitivity_
        nlogL = optimoutput{ii,end}.history.fval(end);
        X = X_new;
    else
        optimoutput = optimoutput(1:end-1,:);
    end
    ii = ii + 1;
end

end

function y = wrapper_function(x,fun,x0,free_indices)
    x_fun = x0;
    x_fun(free_indices) = x;
    y = fun(x_fun);    
end