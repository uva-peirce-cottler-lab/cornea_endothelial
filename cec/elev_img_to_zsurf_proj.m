function [z_surf_proj, z_fill] = elev_img_to_zsurf_proj(xyz_img, elev_xy)
%UNTITLED Summary of this function goes here
%   Takes a 2D elevation image (pixel value is zslice index and performs a
%   pixel by pixel variable maximum intensity z projection starting at the
%   beginning of the zstack on going up to the specified z index

img_dim = size(xyz_img);

% Make row coord image
   col_img = bsxfun(@times, 1:size(elev_xy,2),...
       ones(size(elev_xy,1),size(elev_xy,2)));
   % Make col coord image
   row_img = bsxfun(@times, (1:size(elev_xy,1))',...
       ones(size(elev_xy,1),size(elev_xy,2)));
   
   % Take all nonzero pixels from elevation, row, and col images, they will
   % line up
   z_ind = elev_xy(elev_xy>0);
   c_ind = col_img(elev_xy>0);
   r_ind = row_img(elev_xy>0);
   
   %Convert all images to triplits of coordinates, subtract 2 to elevation
   %ebcause we want to capture the top 2-3 images for a surface z projects
   l_ind = sub2ind(size(xyz_img),r_ind, c_ind,...
       min([max([z_ind+1,ones(size(z_ind))],[],2) img_dim(3)*ones(size(z_ind))],[],2));
   
   
   % Mark all coordinates in new binary image
   z_extent = false(size(xyz_img));
   z_extent(l_ind)=1;
   
   % Compute cumsum with flipped image
   z_fill = cumsum(z_extent,3,'reverse')>0;
   
   z_surf_proj = max(double(z_fill) .*xyz_img,[],3)>0;

end

