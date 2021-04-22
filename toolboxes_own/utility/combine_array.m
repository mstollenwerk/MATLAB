function [S] = combine_array(varargin)
%COMBINE_ARRAY Combine multiple arrays with NaN and nonNaN entries to one.
%   If at any index 

nans = cellfun(@isnan,varargin,'UniformOutput',false);
sizes = cellfun(@size,nans,'UniformOutput',false);

% Do all arrays have the same size?
if ~isequal(sizes{:})
    error('All inputs must have same size.')
end

% Is one index covered by two or more of the inputs?
if any(sum(~cell2mat(nans),2)>1)
    error('An index is covered by more than one input array.')
end

S = NaN(size(nans{1}));
for ii = 1:numel(nans)
  S(~nans{ii}) = varargin{ii}(~nans{ii});
end

end

