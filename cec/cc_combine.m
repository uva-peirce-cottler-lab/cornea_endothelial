function cc_out = cc_combine(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

cc_out = varargin{1};

cc_out.NumObjects = sum(cellfun(@(x) x.NumObjects, varargin));
% keyboard
temp_idx_list = cellfun(@(x) x.PixelIdxList(:),varargin,'UniformOutput',0);
cc_out.PixelIdxList = vertcat(temp_idx_list{:});


end

