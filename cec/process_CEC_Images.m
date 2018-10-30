function process_CEC_Images()

% Load base path for box sync
proj_path = getappdata(0,'proj_path');
mkdir([proj_path '/temp/']);
if isempty(dir([proj_path '/temp/base_path.mat']))
    base_path = uigetdir('Select Box Folder');
    save([proj_path '/temp/base_path.mat'],'base_path');
else
    load([proj_path '/temp/base_path.mat']);
end
in_path = [base_path '\16. CEC Project\Exp_292 CEC My11-TOM Uninduced_TMX'];
out_path = [in_path '\out'];

tif_names_sarr = dir([in_path '/orig*.tif']);


%% Create surface zproj image
for n=1:numel(tif_names_sarr);
    
    fprintf('Loading: %s\n',tif_names_sarr(n).name);
     if ~isempty(dir([out_path '/z_proj_' tif_names_sarr(n).name 'f'])); continue; end
    
    % Load Image
    st = imfinfo([in_path '/' tif_names_sarr(n).name]);
    z_tot = str2double(regexp(st(1).ImageDescription,'slices=(\d*)','tokens','once'));
    ch_tot = str2double(regexp(st(1).ImageDescription,'channels=(\d*)','tokens','once'));
    
    xychz_img = zeros(st(1).Height, st(1).Width,3,z_tot,'uint8');
    k=1;
    for z=1:z_tot
        for ch = 1:ch_tot
            fprintf('Z: %.f, Ch %.f\n', z,ch)
            xychz_img(:,:,ch,z) = ...
                imread([in_path '/' tif_names_sarr(n).name],k);
            k=k+1;
        end
    end
    img_dim = size(xychz_img);
    
    %    keyboard
    if ~isempty(regexp( tif_names_sarr(n).name,'Exp_283.2','once')');
        %        keyboard
        xychz_img(:,:,2,:)=xychz_img(:,:,3,:);
        xychz_img(:,:,3,:)=0;
    end
    if ~isempty(regexp( tif_names_sarr(n).name,'Exp_290','once')');
        %        keyboard
        xychz_img(:,:,2,:)=xychz_img(:,:,3,:);
        xychz_img(:,:,3,:)=0;
    end
    if ~isempty(regexp( tif_names_sarr(n).name,'Exp_283.3','once')');
        %        keyboard
        xychz_img(:,:,2,:)=xychz_img(:,:,3,:);
        xychz_img(:,:,3,:)=0;
    end
    
    
    
    % Seclect second channels for DAPI cell nucleii
    xyz_gs = squeeze(xychz_img(:,:,2,:));
    % Slight blurring of raw input image data
    xyz_bl_gs = imfilter(xyz_gs,fspecial('disk',1));
    
    % Function to enhance contrast
    %     ndadjust = @(x) uint8( double(x-min(x(:))).*double(intmax(class(x)))./double((max(x)-min(x))));
    %     xyz_gs_proc = uint8(subimage_map(xyz_bl_gs, [1024 1024],ndadjust));
    xyz_gs_proc = xyz_bl_gs;
    
    % Threshold 3d image z slice by z slice
    out_args = zslice_map(xyz_gs_proc, @CEC_zslice_thresh);
    bw_nucleii_xyz = out_args{1};
    min_xyz = out_args{2};
    max_xyz = out_args{3};
    clear out_args;
    
    
    % keyboard
    %     ratio_over_roir
    
    % Locate extent of 3d image to z project for a surface z projection
    elev_xy = seg3d_2_elev_img(bw_nucleii_xyz,1);
    z_proj_thresh = elev_img_to_zsurf_proj(bw_nucleii_xyz, elev_xy);
    
    ratio_xyz = bsxfun(@times, max_xyz./min_xyz, z_proj_thresh);
    
    % Find z layer where the greatest contrast between forground and
    % background is found
    [~,z_surf] = max(ratio_xyz,[],3);
    
    %    keyboard
    
    % Inverse z layer so we can run a max filter (for top z layer)
    % TODO need to do this image by image
    z_surf_inv = (z_tot - z_surf).*z_proj_thresh;
    z_surf_inv_max = double(ordfilt2(z_surf_inv,49.^2,ones(49,49)));
    z_surf_max = (z_tot - z_surf_inv_max);
    
    % Makr function given an elevation image and a 3D zstack, max projection
    % the portions of 3d zstack that is determined from elevation image
    % Perform max intensity proj based on z_surf_max level specified for
    % each mixel
    [z_proj_thresh2, z_fill] = elev_img_to_zsurf_proj(bw_nucleii_xyz, z_surf_max);
    
    % Logical indexing with cumsum >0 will yield a z project that only
    % captures the top surface and a little below for each image
    z_proj_ch1 = max( uint8(z_fill) .*squeeze(xychz_img(:,:,1,:)),[],3);
    z_proj_ch2 = max( uint8(z_fill) .*squeeze(xyz_gs_proc),[],3);
    z_proj_ch1_thresh = z_proj_ch1>35;
    
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
    
    %    keyboard
    
    % Remove segments smaller than 250 area, directly keep those in middle
    % size range, and split up those in large range.
    %     RGB = zeros(img_dim(1),img_dim(2),3,'uint8');
    %    RGB(:,:,1) = 256*(labelmatrix(cc_mt_nuc)>0);
    %    RGB(:,:,2) = 256*(labelmatrix(cc_single_nuc)>0);
    %    RGB(:,:,3) = 256*(labelmatrix(cc_lt_nuc)>0);
    %    RGB(:,:,1) = RGB(:,:,1) + RGB(:,:,3);
    %    RGB(:,:,2)= RGB(:,:,2) + RGB(:,:,3);
    
    
    
    
    % Remove spots not large enough
    %    inner_seg_dist = bwdist(~z_proj_thresh2);
    %    inner_seg_dist(inner_seg_dist<7)=0;
    
    
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
    
    

    %Find centroid of each connected component
    ctrs = regionprops(cc_passed_single_nuc, 'centroid');
    % Sample pixels at each centroid, determine if + or - for marker
    rc_ctr_ind = [round(arrayfun(@(x) x.Centroid(1,2),ctrs)),...
        round(arrayfun(@(x) x.Centroid(1,1),ctrs))];
    lind_nucleii_ctrs = sub2ind(img_dim(1:2), rc_ctr_ind(:,1), rc_ctr_ind(:,2));
    nucleii_ctrs=false(img_dim(1:2)); nucleii_ctrs(lind_nucleii_ctrs)=1;
    
    % Determeine if nucleii is marked
    is_nucleii_marked = cellfun(@(x) mean(z_proj_ch1_thresh(x)),cc_passed_single_nuc.PixelIdxList)>.75;
    

    % write output coordinates and classification
    csvwrite([out_path '/' tif_names_sarr(n).name 'coords.csv'],[rc_ctr_ind is_nucleii_marked])
    
    
    
    %    Export surface z projected image
    z_proj = zeros(img_dim(1:3),'uint8');
    z_proj(:,:,1)=z_proj_ch1;
    z_proj(:,:,3)=z_proj_ch2;
    z_proj(:,:,2)=z_proj_thresh2*100 + nucleii_ctrs*100;
    imwrite(z_proj, [out_path '/z_proj_' tif_names_sarr(n).name 'f']);
    
    % Export threshlded image
    rgb_thresh =zeros(img_dim(1:3),'uint8');
    rgb_thresh(:,:,1) = z_proj_ch1_thresh*256;
    rgb = label2rgb(labelmatrix(cc_passed_single_nuc));
    rgb_thresh(:,:,2:3) =  rgb(:,:,2:3);
    imwrite(rgb_thresh, [out_path '/thresh_' tif_names_sarr(n).name 'f']);
end



end