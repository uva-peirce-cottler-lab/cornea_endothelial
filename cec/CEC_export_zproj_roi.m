

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
    
   fprintf('Loading: %s\n',tif_names_sarr(n).name);
  
   % Load zmax image, load roi
   z_proj = imread([out_path '/z_proj_' tif_names_sarr(n).name 'f']);
   st = load([out_path '/' tif_names_sarr(n).name '.mat']);
   
   % Zero out section of zmax
   z_proj(cat(3, ~st.mask, ~st.mask, ~st.mask)) = 0;
   
   imwrite(z_proj,[out_path '/roi_' tif_names_sarr(n).name 'f']);
       
%        close(gcf);
end
