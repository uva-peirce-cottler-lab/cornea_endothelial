function cc_merged = merge_cc_fragments_disk(cc_frags,rad)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% keyboard
% Remove united CCs from b+w frag image
bw_frags = labelmatrix(cc_frags)>0;

nuc_score = imfilter(bw_frags,fspecial('disk',rad),'symmetric');


% bwco
cc = bwconncomp(nuc_score>0.8);
lbl = labelmatrix(cc);

for n=1:cc.NumObjects
    bv_merged_cc_cell = cellfun(@(x) any(lbl(x)==n),cc_frags.PixelIdxList);
    
    temp = cc_frags.PixelIdxList(bv_merged_cc_cell);
    
    
    PixIdxList{n} =  unique(vertcat(vertcat(temp{:}, cc.PixelIdxList{n})));
end


% cc_merged =cc;
cc_merged.PixelIdxList = PixIdxList;

