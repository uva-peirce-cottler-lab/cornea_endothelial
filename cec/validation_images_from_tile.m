function validation_images_from_tile(rgb_gs,rgb_thresh, img_names, ind)
%UNTITLED2 Summary of this function goes here


% Load base path for box sync
proj_path = getappdata(0,'proj_path');
mkdir([proj_path '/temp/']);
if isempty(dir([proj_path '/temp/base_path.mat']))
    base_path = uigetdir('Select Box Folder');
    save([proj_path '/temp/base_path.mat'],'base_path');
else
    load([proj_path '/temp/base_path.mat']);
end
out_path = [base_path '\16. CEC Project\validation_images'];

mkdir(out_path);

% Sub image index from greater image tile
c_sub = 1024*3+150+(1:512);
r_sub = 1024*3+150+(1:512);

% keyboard
cc_cropped_nuclei = segment_2d_nuclei(bwareaopen(rgb_thresh(r_sub,c_sub,3),25));



img_dim = cc_cropped_nuclei.ImageSize;

% keyboard


nuclei_thresh = labelmatrix(cc_cropped_nuclei)>0;
cell_thresh = rgb_thresh(r_sub,c_sub,1);

%   Detailed explanation goes here
%  Recalculate marked stat of 
ctrs = regionprops(cc_cropped_nuclei, 'centroid');
% Sample pixels at each centroid, determine if + or - for marker
rc_ctr_ind = [round(arrayfun(@(x) x.Centroid(1,2),ctrs)),...
    round(arrayfun(@(x) x.Centroid(1,1),ctrs))];
lind_nucleii_ctrs = sub2ind(img_dim(1:2), rc_ctr_ind(:,1), rc_ctr_ind(:,2));
nucleii_ctrs=false(img_dim(1:2)); nucleii_ctrs(lind_nucleii_ctrs)=1;
% Determeine if nucleii is marked
is_nucleii_marked = cellfun(@(x) mean(cell_thresh(x)),cc_cropped_nuclei.PixelIdxList)>.75;
cc_marked=cc_cropped_nuclei; cc_marked.NumObjects = sum(is_nucleii_marked);
cc_marked.PixelIdxList = cc_cropped_nuclei.PixelIdxList(is_nucleii_marked);
cc_unmarked=cc_cropped_nuclei; cc_unmarked.NumObjects = sum(~is_nucleii_marked);
cc_unmarked.PixelIdxList = cc_cropped_nuclei.PixelIdxList(~is_nucleii_marked);
nuc_mark_state_gs = uint8(255.*(labelmatrix(cc_marked)>0) + 100.*(labelmatrix(cc_unmarked)>0));


% nuc_mark_state_gs = nuc_mark_state_gs(c_sub,r_sub,:);

csv_path = [out_path '/img_tbl.csv'];

if ~isempty(dir(csv_path))
    tbl = readtable(csv_path);
else
    tbl = table();
    tbl.img_name = img_names;
    tbl.auto_count_rfp_minus = zeros(numel(img_names),1);
    tbl.auto_count_rfp_plus = zeros(numel(img_names),1);
    tbl.manual_count_rfp_minus = zeros(numel(img_names),1);
    tbl.manual_count_rfp_plus = zeros(numel(img_names),1);
    
end

% Quantify particular entry and update table
tbl.auto_count_rfp_minus(ind) = sum(~is_nucleii_marked);
tbl.auto_count_rfp_plus(ind) = sum(is_nucleii_marked);

% Update csv table
writetable(tbl, csv_path);

% write image to disk
imwrite(rgb_gs(r_sub, c_sub,:),[out_path '/' img_names{ind} '_valid.tif']);

% keyboard
end

