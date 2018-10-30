function cc_merged = merge_cc_fragments(cc_frags, area_range)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Remove united CCs from b+w frag image
bw_frags = labelmatrix(cc_frags)>0;
% Paint current fragment and then clear at end of each loop
bw_curr_frag = false(cc_frags.ImageSize);
keyboard

lind_bw = 1:prod(cc_frags.ImageSize);
lind_ccs = 1:cc_frags.NumObjects;
for n=1:cc_frags
    bw_curr_frag(cc_frags.PixelIdxList{n})=1;
    bw_outer_border = imdilate(bw_curr_frag,strel('disk',2,0))&~bw_curr_frag;
    lind_border_coords = lind_bw(bw_outer_border);
    
    % Find any CCs that border current CC
    nint = cellfun(@(x) numel(intersect(lind_border_coords,x)),cc_frags.PixelIdxList);
    
    % Get linear index of adjacent CCs
    lind_adj_ccs = lind_ccs(nint>0);
    % Get pixel area
    npix_adj_ccs = cellfun(@(x) numel(x), cc_frags.PixelIdxList(lind_adj_ccs));
     
    % Sort based on area of each CC
    [sorted_npix_adj_ccs,ix] = sort(npix_adj_ccs, 'ascend');
    sorted_lind_adj_ccs = lind_adj_ccs(ix);
    
    % Calculate how much area would be added by iteratively merging each CC
    cc_added_areas = numel(cc_frags.PixelIdxList{n}) + cumsum(sorted_npix_adj_ccs);
    % Find index of added list that would make combined area in acceptable
    % range
    bv = cc_added_areas > area_range(1) & cc_added_areas < area_range(2); 
    index_bv = 1:numel(bv);
    lind_of_bv = index_bv(bv);
    
    % Take first area that puts nucleus in acceptable range
    last_lind_cc_to_add = lind_of_bv(1);
    
    linds_cc_to_add = last_lind_to_add
    
    % For each bordering CC, iteratively add CCs to current until area
    % requirement is met
    
     % Ad to cc_merge
     
     % Remove CC from PixelIdxList
     
     % Reomve from bw_frags
     
     % clear bw_frags
     
end



end

