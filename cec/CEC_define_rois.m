

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



for n=1:numel(tif_names_sarr); 
    
 
   
   % Load Image
   img_path = [out_path '/z_proj_' tif_names_sarr(n).name 'f'];
   fprintf('Loading: %s\n',tif_names_sarr(n).name);
%    st = imfinfo(img_path)
   

       % Draw polygon around cornea
   mask =[]; mask_coords=[];
   if ~isempty(dir([in_path '/' tif_names_sarr(n).name '.mat']));
       continue;
       load([in_path '/' tif_names_sarr(n).name '.mat']);
   end
   
   % load image
   xych_img = imread(img_path);
   
   % st.mask
   % st.mask_coords
   
   % Display image and prompt user for polgyon
   hf = imshow(xych_img);
%   keyboard
    h=impoly(gca,mask_coords);
    k=0;
    while k~=1
        k = ~ishandle(hf); 
        if ~k 
            mask = h.createMask();
            mask_coords = h.getPosition;
        end
        pause(.25); 
    end
       
       % Get mask and coordinates, save to disk
       
       
       save([out_path '/' tif_names_sarr(n).name '.mat'],'mask','mask_coords');
%        imwrite(bw_cornea, [in_path '/' base_name '_bw.tif']);
       close(gcf);
end
