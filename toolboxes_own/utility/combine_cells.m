function [y] = combine_cells(varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% y = cat(2, varargin{:});
% y = y(~cellfun('isempty', y));

tmp = cat(3,varargin{:});
[~,idp] = sort(cellfun('isempty',tmp),3);
sz = size(tmp);
[idr,idc] = ndgrid(1:sz(1),1:sz(2));
idx = sub2ind(sz,idr,idc,idp(:,:,1));
y = tmp(idx);
end

