function [x0,optimoutput] = my_fminunc_partial(fun,x0,free_indices,varargin)

for ii = 1:length(free_indices)
    
    if ~isempty(free_indices{ii})
        fun_ii = @(x) wrapper_function(x,fun,x0,free_indices{ii});
        [x0(free_indices{ii}),optimoutput{ii}] = ...
                my_fminunc(fun_ii,x0(free_indices{ii}),varargin{:});
    end
    
end

end

function y = wrapper_function(x,fun,x0,free_indices)
    x_fun = x0;
    x_fun(free_indices) = x;
    y = fun(x_fun);    
end