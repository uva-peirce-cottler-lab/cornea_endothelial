



c_sub = 1024*3+256+(1:350);
r_sub = 1024*3+256+400+(1:350);

raw_img = max(xyz_gs(c_sub,r_sub,:),[],3);

first_thresh = max(bw_nucleii_xyz(c_sub,r_sub,:),[],3) ;

first_elev_xy = max(elev_xy(c_sub,r_sub,:),[],3) ;

first_z_proj_thresh = max(z_proj_thresh(c_sub,r_sub,:),[],3) ;


z_proj_fig = max(z_proj(c_sub,r_sub,:),[],3) ;

rgb_thresh_fig = max( rgb_thresh(c_sub,r_sub,:),[],3) ;

 
imshow(raw_img)
imshow(first_thresh)
imshow(first_elev_xy)
imshow(first_z_proj_thresh)

imshow(z_proj_fig)
imtool( first_elev_xy)
