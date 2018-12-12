

img_path = 'C:\Users\bruce\Box Sync\16. CEC Project\Richard 3D CellNucleii VIsual\60_zstacj_cornea_cec_dapi3.tif';
info = imfinfo(img_path);

img = zeros(info(1).Width,info(1).Height,3,ceil(length(info)/2),'uint8');

num_slices = length(info);
for k = 1:2:num_slices-1
  img(:,:,3,ceil(k/2)) = imread(img_path,'Info',info,'Index',k);
  img(:,:,1,ceil(k/2)) = imread(img_path,'Info',info,'Index',k+1);
end

img2 = img(200:200+712, 200:200+712,:,:);
size(img2)

img3 = permute(img2,[1 4 3 2]);

size(img3)

x_proj = imresize( max(img3,[],4), [size(img2,1) size(img2,1)]);

imwrite(x_proj,'x_proj.tif')