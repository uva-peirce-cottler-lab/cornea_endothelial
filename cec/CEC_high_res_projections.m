




% Load base path for box sync
proj_path = getappdata(0,'proj_path');
mkdir([proj_path '/temp/']);
if isempty(dir([proj_path '/temp/base_path.mat']))
    base_path = uigetdir('Select Box Folder');
    save([proj_path '/temp/base_path.mat'],'base_path');
else
    load([proj_path '/temp/base_path.mat']);
end
in_path = [base_path '\16. CEC Project\Richard 3D CellNucleii VIsual'];

% tif_names_sarr = dir([in_path '/']);
tif_name = '60x_dapi_cornea_layers_all.lsm';


 [raw_xychz_img, meta] = img_open([in_path '/' tif_name]);
 
 % X axis project to show layers in image
 xyz_img = squeeze(raw_xychz_img); 
 yzx = permute(xyz_img,[2 3 1]);
 x_proj = imrotate(max(yzx(:,:,1:256),[],3),90);
 x_proj_rgb = zeros(size(x_proj),'uint8');
 x_proj_rgb(:,:,3) = x_proj;
 x_proj_out  = imresize(x_proj_rgb, [1024 1024]);
 imwrite(x_proj_out,[in_path '/' tif_name 'x_proj.tif'])

 % Pseudo color image to a specified brightest target value
 % use LUT?
 B = intlut(A, LUT)
bottom_layer_img = imread([in_path '/bot_z_proj.tif']);

mid_layer_img = imread([in_path '/mid_z_proj.tif']);

top_layer_img = imread([in_path '/top_z_proj.tif']);


help
