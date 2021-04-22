function [x0,optimoutput] = my_fminunc_partial_old(fun,x0,free_indices,varargin)

numelx0 = numel(x0);
for ii = 1:length(free_indices)
    fixed_indices = setdiff(1:numelx0, free_indices{ii});
    [Aeq, beq] = equality_restrictions(fixed_indices,x0);
    [x0,optimoutput{ii}] = ...
            my_fmincon(fun,x0,[],[],Aeq,beq,[],[],[],varargin{:});
end

end

