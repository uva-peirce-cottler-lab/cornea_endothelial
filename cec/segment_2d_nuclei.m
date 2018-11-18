function cc_passed_single_nuc = segment_2d_nuclei(z_proj_thresh2);


    %    medfilt(xyz_gs_proc,[5 5],'symmetric')
    cc_seg = bwconncomp(z_proj_thresh2);
    
    % Split CCs bvased on area, discard too small, keep just right, and save
    % too large to split
    nuc_low = 200;
    nuc_high = 550;
    area_cc = cellfun(@(x) numel(x), cc_seg.PixelIdxList);
    
    cc_lt_nuc = cc_subset(cc_seg,area_cc<200);
    
    cc_nuc = cc_subset(cc_seg,area_cc>=200 & area_cc<550);
    
    cc_mt_nuc = cc_subset(cc_seg,area_cc>=550);
    bw_mt_nuc = labelmatrix(cc_mt_nuc);

    % Perform euclidean distance inside and outside of nucleii, combine for
    % segment distance elevation map
    sm_nuc_border = imclose(imopen(bw_mt_nuc,strel('disk',2)),strel('disk',2));
    inner_elev = bwdist(~sm_nuc_border);
    outer_elev = bwdist(sm_nuc_border);
    seg_elev = outer_elev-inner_elev;
    %    inner_ed_elev= bwdist(inner_seg_dist>0);
    
    % Perform watershed on segment distance elevation map
    ws = watershed(imopen(seg_elev,strel('disk',5)));
    
    % Use wastershed boundaries to split up segmentaion
    nuc_sep = bw_mt_nuc;
    %    nuc_sep( ~ws)=0;
    nuc_sep(~ws)=0;
    
    % Find nucleii of acceptable size
    cc_splitted = bwconncomp(nuc_sep);
    cc_splitted_area = cellfun(@(x) numel(x),cc_splitted.PixelIdxList);
    % Additional nucleii from splitting connected componenets
    cc_splitted_single_nuc = cc_subset(cc_splitted, cc_splitted_area >= nuc_low...
        & cc_splitted_area <= nuc_high);
    % Try tp combine smalled connected components iteratively to get nucleii
    % within specified size range
    cc_splitted_lt_nuc = cc_subset(cc_splitted, cc_splitted_area < nuc_low);
    
    
    cc_united_nuc_frags = merge_cc_fragments_disk(cc_splitted_lt_nuc, 8);
    cc_united_nuc_frags_area = cellfun(@(x) numel(x), cc_united_nuc_frags.PixelIdxList);
    
    cc_united_single_nuc = cc_subset(cc_united_nuc_frags, cc_united_nuc_frags_area >= nuc_low &  ...
        cc_united_nuc_frags_area <= nuc_high);
    
%     keyboard
    cc_passed_single_nuc = cc_combine(cc_nuc,cc_splitted_single_nuc,cc_united_single_nuc);
    