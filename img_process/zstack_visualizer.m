

 [raw_xychz_img, meta] = img_open(['C:\Users\bac\Box Sync\13. CBurn\Exp_274 Low Mag Image\3_cornea.lsm']);
 
 fv = isosurface(squeeze(raw_xychz_img(:,:,1,:)),0) ;
 
 
 p = patch(isosurface(x,y,z,v,-3));
isonormals(x,y,z,v,p)
