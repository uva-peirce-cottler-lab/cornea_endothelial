function [fs_proc_elev_xy,fs_flipped_elev_xy] = seg3d_2_elev_img(xyz_img,DownSample_Image)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%     keyboard
img_dim = size(xyz_img);


% Downsample image if requested
if DownSample_Image
    rs_xyz_img = zeros(img_dim.*[0.25 0.25 1]);
    for z = 1:size(xyz_img,3)
        rs_xyz_img(:,:,z) = imresize(bwareaopen(xyz_img(:,:,z),10),1/4);
    end
else
    rs_xyz_img = xyz_img;
end

% Go through each image, find vertical index at each white
elev_xyz = uint8(bsxfun(@times, permute(1:size(rs_xyz_img,3),[3 1 2]),...
    ones(size(rs_xyz_img))) .*rs_xyz_img);
% XY pixels where no thresholded pixels (and elevation values) found
no_elev_xy = sum(elev_xyz,3)==0;

% Elevation is zslice of true pixel closest to beginning of zstack
elev_xy = min(elev_xyz+uint8((elev_xyz==0).*img_dim(3)), [],3);
elev_xy(no_elev_xy)=0;

% Flip elevation so that small zslice is a big # that survives dilation
flipped_elev_xy = (img_dim(3) - elev_xy);
flipped_elev_xy(no_elev_xy)=0;

% Perform additional filtering on elevation image on an image by image
% basis (sinve elevation can change abrubtly between images)
%    keyboard
% 150

%    proc_flipped_elev_xy = imdilate(flipped_elev_xy,strel('disk',100,0));

proc_flipped_elev_xy = ordfilt2(flipped_elev_xy,100^2,ones(100,100),'symmetric');

% Unflip elevation index
proc_elev_xy = (img_dim(3) - proc_flipped_elev_xy);



% Convert elevation image to fullsize
fs_proc_elev_xy = imresize(uint8(proc_elev_xy),img_dim([1 2]));
fs_flipped_elev_xy = imresize(uint8(flipped_elev_xy),img_dim([1 2]));

% z_surf_proj = elev_img_to_surf_proj(xyz_img, fs_proc_elev_xy);
%     keyboard

end

