function ncc = cc_subset(cc,idx)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ncc = cc;
ncc.PixelIdxList = cc.PixelIdxList(idx);
ncc.NumObjects = sum(idx);
end

