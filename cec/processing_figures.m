% function processing_figures


suff = [num2str(n) '_'];
out_path = [getappdata(0,'proj_path') '\temp\process_images'];
mkdir(out_path);

% n=11

elev_uint8=@(x) uint8(50+x.* 200./size(xyz_gs,3));

% Sub image index from greater image tile
c_sub = 1024*3+256+(1:350);
r_sub = 1024*3+256+400+(1:350);

% Greyscale raw dapi image
gs_dapi_img = max(xyz_gs(c_sub,r_sub,:),[],3);
imwrite(gs_dapi_img,[out_path '/' suff 'gs_dapi_img.tif']);

% Initial 3D thresh
first_thresh = max(bw_nucleii_xyz(c_sub,r_sub,:),[],3) ;
imwrite(first_thresh,[out_path '/' suff 'initial_3d_thresh.tif']);

% Initial unbridged elevtation image
initial_unbridged_elev_xy = unbridged_elev_xy(c_sub,r_sub,:);
imwrite(elev_uint8(initial_unbridged_elev_xy),[out_path '/' suff 'initial_elev_xy_unbridged.tif']);

% Initial bridged Elevation image
initial_bridged_elev_xy = bridged_elev_xy(c_sub,r_sub,:);
imwrite(elev_uint8(initial_bridged_elev_xy),[out_path '/' suff 'initial_elev_xy_bridged.tif']);

% Initial Surface Z projection threshold
initial_zproj_surf_thresh = max(z_proj_thresh_1(c_sub,r_sub,:),[],3) ;
imwrite(initial_zproj_surf_thresh,[out_path '/' suff 'initial_zproj_surf_thresh.tif']);

% Initial surface z projection of thresh
% z_proj_fig = max(z_proj(c_sub,r_sub,:),[],3) ;
% imwrite(z_proj_fig,[out_path '/initial_z_surf_proj.tif']);

% Initial surface z projection of grayscale image
initial_zproj_surf_gs = max(z_proj_gs_1(c_sub,r_sub,:),[],3);
imwrite(initial_zproj_surf_gs,[out_path '/' suff 'initial_zproj_surf_gs.tif']);


%% Second pass
% Second unbriged elevation image
second_unbridged_elev_xy = z_surf(c_sub,r_sub,:);
imwrite(elev_uint8(second_unbridged_elev_xy),[out_path '/' suff 'second_unbridged_elev_xy.tif']);

% Second bridged elevation image
second_bridged_elev_xy = z_surf_max(c_sub,r_sub,:);
imwrite(elev_uint8(second_bridged_elev_xy),[out_path '/' suff 'second_bridged_elev_xy.tif']);

% Second z_projection thresh
second_z_proj_thresh = z_proj_thresh2(c_sub,r_sub,:);
imwrite(second_z_proj_thresh,[out_path '/' suff 'second_zproj_surf_thresh.tif']);

% Second z_projection grayscale
second_zproj_surf_gs = z_proj_ch2(c_sub,r_sub,:);
imwrite(second_zproj_surf_gs,[out_path '/' suff 'second_zproj_surf_gs.tif']);

% Final segmented nucleotides
cc_passed_single_nuc.PixelIdxList = ...
    cc_passed_single_nuc.PixelIdxList(randperm(cc_passed_single_nuc.NumObjects));
RGB = label2rgb(labelmatrix(cc_passed_single_nuc));
final_nuc_seg = RGB(c_sub,r_sub,:);
imwrite(final_nuc_seg,[out_path '/' suff 'final_nuc_seg.tif']);


%Second z surface projection rgb
second_z_surf_r = max( uint8(z_proj_fill2) .*squeeze(xychz_img(:,:,1,:)),[],3);
second_z_surf_b = imadjust(max( uint8(z_proj_fill2) .*squeeze(xychz_img(:,:,2,:)),[],3));
second_z_surf_rgb = cat(3, second_z_surf_r,...
    zeros(img_dim(1:2),'uint8'),second_z_surf_b);
second_z_surf_proj_rgb = second_z_surf_rgb(c_sub,r_sub,:);
imwrite(second_z_surf_proj_rgb,[out_path '/' suff 'second_z_surf_proj_rgb.tif']);

% RFP z projection threshold
second_z_proj_thresh_cell_marker = z_proj_ch1_thresh(c_sub,r_sub,:);
imwrite(second_z_proj_thresh_cell_marker,[out_path '/' suff 'second_z_proj_thresh_cell_marker.tif']);
return

% Nuclei with/without rfp labeling
% Recalculate marked stat of 
ctrs = regionprops(cc_passed_single_nuc, 'centroid');
% Sample pixels at each centroid, determine if + or - for marker
rc_ctr_ind = [round(arrayfun(@(x) x.Centroid(1,2),ctrs)),...
    round(arrayfun(@(x) x.Centroid(1,1),ctrs))];
lind_nucleii_ctrs = sub2ind(img_dim(1:2), rc_ctr_ind(:,1), rc_ctr_ind(:,2));
nucleii_ctrs=false(img_dim(1:2)); nucleii_ctrs(lind_nucleii_ctrs)=1;
% Determeine if nucleii is marked
is_nucleii_marked = cellfun(@(x) mean(z_proj_ch1_thresh(x)),cc_passed_single_nuc.PixelIdxList)>.75;
cc_marked=cc_passed_single_nuc; cc_marked.NumObjects = sum(is_nucleii_marked);
cc_marked.PixelIdxList = cc_passed_single_nuc.PixelIdxList(is_nucleii_marked);
cc_unmarked=cc_passed_single_nuc; cc_unmarked.NumObjects = sum(~is_nucleii_marked);
cc_unmarked.PixelIdxList = cc_passed_single_nuc.PixelIdxList(~is_nucleii_marked);
nuc_mark_state_gs_tile = uint8(255.*(labelmatrix(cc_marked)>0) + 100.*(labelmatrix(cc_unmarked)>0));
nuc_mark_state_gs = nuc_mark_state_gs_tile(c_sub,r_sub,:);
imwrite(nuc_mark_state_gs,[out_path '/' suff 'nuc_mark_state_gs.tif']);

% Color bar for elevation image
wd = 20; n_slices = size(xyz_gs,3);
delta = floor(350./n_slices);
 color_bar = zeros(n_slices*delta,wd,'uint8');
 for z=1:n_slices
     color_bar(1+delta*(z-1):delta*z,:)=50+z*200./n_slices;
 end
 imwrite(color_bar,[out_path '/' suff 'color_bar.tif']);

 
 % 3D isosurface
 bw_3d = cat(3,bwareaopen(bw_nucleii_xyz(c_sub,r_sub,:),20),false(350,350));
 isosurface(bw_3d ,1/2);
 
% Plot of ratio values foreground/background 
 ratio_xyz_vector = ratio_xyz(c_sub,r_sub,:);
 ratio_xyz_col = reshape(ratio_xyz_vector,[size(ratio_xyz_vector,1).*...
     size(ratio_xyz_vector,2) size(ratio_xyz_vector,3)]);
nonzero_ratio  = ratio_xyz_col(sum(ratio_xyz_col>0,2)>0,:);
 plot(nonzero_ratio(1:100,:)')
 xr = xlim; yr = ylim;
 axis([1 xr(2) yr])
 
 
 