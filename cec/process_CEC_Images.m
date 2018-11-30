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
in_path = [base_path '\16. CEC Project\input_images'];
out_path = [in_path '\out'];

tif_names_sarr = dir([in_path '/orig*.tif']);


%% Create surface zproj image
for n=1:numel(tif_names_sarr);
    
    fprintf('Loading: %s\n',tif_names_sarr(n).name);
%      if ~isempty(dir([out_path '/z_proj_' tif_names_sarr(n).name 'f'])); continue; end
    
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
    
    % Switch channels according to experiment, they were taken with
    % different channel configurations
    if ~isempty(regexp( tif_names_sarr(n).name,'Exp_283.2','once')');
        xychz_img(:,:,2,:)=xychz_img(:,:,3,:);
        xychz_img(:,:,3,:)=0;
    end
    if ~isempty(regexp( tif_names_sarr(n).name,'Exp_290','once')');
        xychz_img(:,:,2,:)=xychz_img(:,:,3,:);
        xychz_img(:,:,3,:)=0;
    end
    if ~isempty(regexp( tif_names_sarr(n).name,'Exp_283.3','once')');
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
    
    % Locate extent of 3d image to z project for a surface z projection
    [bridged_elev_xy, unbridged_elev_xy] = seg3d_2_elev_img(bw_nucleii_xyz,1);
    [z_proj_thresh_1, z_proj_fill_1] = elev_img_to_zsurf_proj(bw_nucleii_xyz, bridged_elev_xy);
    
    z_proj_gs_xyz = xyz_gs;
    z_proj_gs_xyz(~z_proj_fill_1)=0;
    z_proj_gs_1 = max(z_proj_gs_xyz,[],3);
    
    ratio_xyz = bsxfun(@times, max_xyz./min_xyz, z_proj_thresh_1);
    
    % Find z layer where the greatest contrast between forground and
    % background is found
    [~,z_surf] = max(ratio_xyz,[],3);
    
    %    keyboard
    
    % Inverse z layer so we can run a max filter (for top z layer)
    % TODO need to do this image by image
    z_surf_inv = (z_tot - z_surf).*z_proj_thresh_1;
    z_surf_inv_max = double(ordfilt2(z_surf_inv,49.^2,ones(49,49)));
    z_surf_max = (z_tot - z_surf_inv_max);
    
    % Make function given an elevation image and a 3D zstack, max projection
    % the portions of 3d zstack that is determined from elevation image
    % Perform max intensity proj based on z_surf_max level specified for
    % each mixel
    [z_proj_thresh2, z_proj_fill2] = elev_img_to_zsurf_proj(bw_nucleii_xyz, z_surf_max);
    
    
    % Logical indexing with cumsum >0 will yield a z project that only
    % captures the top surface and a little below for each image
    z_proj_ch1 = max( uint8(z_proj_fill2) .*squeeze(xychz_img(:,:,1,:)),[],3);
    z_proj_ch2 = max( uint8(z_proj_fill2) .*squeeze(xyz_gs_proc),[],3);
    z_proj_ch1_thresh = z_proj_ch1>35;
    
    
    cc_passed_single_nuc = segment_2d_nuclei(z_proj_thresh2);

    %     c_sub = 1024*3+256+(1:350);
    % r_sub = 1024*3+256+400+(1:350);
    
    %Find centroid of each connected component
    ctrs = regionprops(cc_passed_single_nuc, 'centroid');
    % Sample pixels at each centroid, determine if + or - for marker
    rc_ctr_ind = [round(arrayfun(@(x) x.Centroid(1,2),ctrs)),...
        round(arrayfun(@(x) x.Centroid(1,1),ctrs))];
    lind_nucleii_ctrs = sub2ind(img_dim(1:2), rc_ctr_ind(:,1), rc_ctr_ind(:,2));
    nucleii_ctrs=false(img_dim(1:2)); nucleii_ctrs(lind_nucleii_ctrs)=1;
    
    % Determeine if nucleii is marked
    is_nucleii_marked = cellfun(@(x) mean(z_proj_ch1_thresh(x)),cc_passed_single_nuc.PixelIdxList)>.75;
    
    % Write output coordinates and classification
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
    rgb_thresh(:,:,3) = z_proj_thresh2;
    rgb = label2rgb(labelmatrix(cc_passed_single_nuc));
    rgb_thresh(:,:,2) =  rgb(:,:,2);
    imwrite(rgb_thresh, [out_path '/thresh_' tif_names_sarr(n).name 'f']);



    % Export basic validation images
     validation_images_from_tile(z_proj,rgb_thresh>0, {tif_names_sarr(:).name}', n);
        
     
     processing_figures
end



end