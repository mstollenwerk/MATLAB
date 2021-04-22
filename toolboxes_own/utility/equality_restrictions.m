function [Aeq, beq] = equality_restrictions(fixed_indices,param)

I = eye(length(param));
Aeq = I(fixed_indices,:);
beq = param(fixed_indices);

end

