% function [outputArg1,outputArg2] = CEC_compile_results(inputArg1,inputArg2)


% Load base path for box sync
proj_path = getappdata(0,'proj_path');
mkdir([proj_path '/temp/']);
if isempty(dir([proj_path '/temp/base_path.mat']))
    base_path = uigetdir('Select Box Folder');
    save([proj_path '/temp/base_path.mat'],'base_path');
else
    load([proj_path '/temp/base_path.mat']);
end
in_path = [base_path '\16. CEC Project\input_images'];
out_path = [in_path '\out'];

tif_names_sarr = dir([in_path '/orig*.tif']);
% tif_names = {tif_names_sarr(:).name};

% ex_str = "orig_TMX Midterm_EyeDrop_Ind,Exp_283.2 My11TOM Eyedrop D21 CEC Tile,1L_cornea";
% tif_names(cellfun(@(x) ~isempty(regexp(x,ex_str,'once')), tif_names))'


% parse names to save master csv file
% nuc_tbl = table('VariableNames',{'r_coord','c_coord','class','norm_rad_from_outer'});

meta_tbl = table();

meta_tbl.tx_type = cellfun(@(x) regexp(x,'[Tt]erm_([^_*]*)_','tokens','once'),{tif_names_sarr(:).name}');
meta_tbl.term =  cellfun(@(x) regexp(x,'\s(.*)[Tt]erm','tokens','once'),{tif_names_sarr(:).name}');
meta_tbl.exp = cellfun(@(x) str2double(regexp(x,'Exp_([0-9\.]*)','tokens','once')),{tif_names_sarr(:).name}');
meta_tbl.mouse_num = cellfun(@(x) str2double(regexp(x,'(\d*)[LR]_','tokens','once')),{tif_names_sarr(:).name}');
meta_tbl.eye_ind = cellfun(@(x) isempty(regexp(x,'\d*R_','once')),{tif_names_sarr(:).name}')+1;

for n=1:numel(tif_names_sarr)
    val =  str2double(regexp(tif_names_sarr(n).name,'D(\d*)\s','once','tokens'));
    if isempty(val);
        val =  str2double(regexp(tif_names_sarr(n).name,'(\d*)HR','once','tokens'));
    end
    if ~isempty(regexp(tif_names_sarr(n).name,'Late_IP_TMX','once'))
        val =  16*7;
    end
    if ~isempty(regexp(tif_names_sarr(n).name,'STD_TMX','once'))
        val =  4*7;
    end
    
    meta_tbl.day(n) = val;
end

% temp=cellfun(@(x) str2double(regexp(x,'D(\d*)\s','once','tokens')),{tif_names_sarr(:).name}','UniformOutput',0);
 

for n=1:numel(tif_names_sarr)
   %load mask
   st = load([out_path '/' tif_names_sarr(n).name '.mat']);
   % load csv file 
   rc_class = csvread([out_path '/' tif_names_sarr(n).name 'coords.csv']);
   
   
%    culled_rc_marked= rc_marked(~st.mask(rc_marked))=[];
   
   % calculate normlaized distance from edge of polygoned ROI
   CH = bwconvhull(st.mask);
   ed = bwdist(~CH);
   norm_ed = ed./max(ed(:));
   
   % determine distance for each nucleii center
   lind = sub2ind(size(st.mask),rc_class(:,1), rc_class(:,2));
   culled_rc_marked= rc_class(CH(lind'),:);

   norm_rad_from_outer = norm_ed(sub2ind(size(norm_ed),culled_rc_marked(:,1), culled_rc_marked(:,2)));

%    out_matrix{n} = horzcat(rc_class,norm_rad_from_outer);  
   
   meta_tbl.tot_cells(n) = size(rc_class,1);
   meta_tbl.marked_cells_tot(n) = sum(rc_class(:,3));
   meta_tbl.frac_marked_cells(n) = meta_tbl.marked_cells_tot(n)./meta_tbl.tot_cells(n);
   meta_tbl.unmarked_cells_tot(n) = sum(~rc_class(:,3));
   ixq =  norm_rad_from_outer>=0 & norm_rad_from_outer<1/5;
   meta_tbl.frac_mark_cells_q1(n) = sum(rc_class(ixq,3))./sum(ixq);
   
   ixq =  norm_rad_from_outer>=1/5 & norm_rad_from_outer<2/5;
   meta_tbl.frac_mark_cells_q2(n) = sum(rc_class(ixq,3))./sum(ixq);
   
   ixq = norm_rad_from_outer>=2/5 & norm_rad_from_outer<3/5;
   meta_tbl.frac_mark_cells_q3(n) = sum(rc_class(ixq,3))./sum(ixq);
   
   ixq = norm_rad_from_outer>=3/5 & norm_rad_from_outer<4/5; 
   meta_tbl.frac_mark_cells_q4(n) = sum(rc_class(ixq,3))./sum(ixq); 
   
   ixq = norm_rad_from_outer>=4/5 & norm_rad_from_outer<5/5;  
   meta_tbl.frac_mark_cells_q5(n) = sum(rc_class(ixq,3))./sum(ixq); 
   
    
end

fig_pos = [500 500 200 200];
save([in_path '/CEC_results.mat']);

load([in_path '/CEC_results.mat']);
keyboard
cell_count_field = 'marked_cells_tot';


%% IP Induction Experiment
hf = figure('Units','Pixels');
ix1=strcmp(meta_tbl.tx_type,'IP') & ...
    strcmp(meta_tbl.term,'Long') & meta_tbl.day==28;
ix2=strcmp(meta_tbl.tx_type,'IP') & ...
    strcmp(meta_tbl.term,'Long') & meta_tbl.day==112;

plot_cell_count(meta_tbl.(cell_count_field)(ix1),...
    meta_tbl.(cell_count_field)(ix2),'RFP+ Cell Count',{'6','16'})
set(gcf,'position',[200 200 210 180]);

plot_radial_density(0:.25:1,[meta_tbl.frac_mark_cells_q1(ix1) meta_tbl.frac_mark_cells_q2(ix1) ...
    meta_tbl.frac_mark_cells_q3(ix1) meta_tbl.frac_mark_cells_q4(ix1) ...
    meta_tbl.frac_mark_cells_q5(ix1)], ...
    'Relative Rad. to Center','Fraction of RFP+ Cells','r')
plot_radial_density(0:.25:1,[meta_tbl.frac_mark_cells_q1(ix2) meta_tbl.frac_mark_cells_q2(ix2) ...
    meta_tbl.frac_mark_cells_q3(ix2) meta_tbl.frac_mark_cells_q4(ix2) ...
    meta_tbl.frac_mark_cells_q5(ix2)], ...
    'Relative Rad. to Center','Fraction of RFP+ Cells','b')



%% Eyedrop Induction short Experiment
hf = figure('Units','Pixels');
ix1 = strcmp(meta_tbl.tx_type,'EyeDrop') & ...
    strcmp(meta_tbl.term,'Short') & meta_tbl.day==24;
ix2 = strcmp(meta_tbl.tx_type,'EyeDrop') & ...
    strcmp(meta_tbl.term,'Short') & meta_tbl.day==48;

plot_cell_count(meta_tbl.(cell_count_field)(ix1),...
    meta_tbl.(cell_count_field)(ix2),'RFP+ Cell Count')


plot_radial_density(0:.25:1,[meta_tbl.frac_mark_cells_q1(ix1) meta_tbl.frac_mark_cells_q2(ix1) ...
    meta_tbl.frac_mark_cells_q3(ix1) meta_tbl.frac_mark_cells_q4(ix1) ...
    meta_tbl.frac_mark_cells_q5(ix1)], ...
    'Relative Rad. to Center','Fraction of RFP+ Cells','r')
plot_radial_density(0:.25:1,[meta_tbl.frac_mark_cells_q1(ix2) meta_tbl.frac_mark_cells_q2(ix2) ...
    meta_tbl.frac_mark_cells_q3(ix2) meta_tbl.frac_mark_cells_q4(ix2) ...
    meta_tbl.frac_mark_cells_q5(ix2)], ...
    'Relative Rad. to Center','Fraction of RFP+ Cells','r')


%% Eyedrop Induction mid Experiment
hf = figure('Units','Pixels');
ix1 = strcmp(meta_tbl.tx_type,'EyeDrop') & ...
    strcmp(meta_tbl.term,'Mid') & meta_tbl.day==2;
ix2 = strcmp(meta_tbl.tx_type,'EyeDrop') & ...
    strcmp(meta_tbl.term,'Mid') & meta_tbl.day==21;

plot_cell_count(meta_tbl.(cell_count_field)(ix1),...
    meta_tbl.(cell_count_field)(ix2),'RFP+ Cell Count',{'2','21'})
set(gcf,'position',[200 200 210 180]);

plot_radial_density(0:.25:1,[meta_tbl.frac_mark_cells_q1(ix1) meta_tbl.frac_mark_cells_q2(ix1) ...
    meta_tbl.frac_mark_cells_q3(ix1) meta_tbl.frac_mark_cells_q4(ix1) ...
    meta_tbl.frac_mark_cells_q5(ix1)], ...
    'Relative Rad. to Center','Fraction of RFP+ Cells','r')
plot_radial_density(0:.25:1,[meta_tbl.frac_mark_cells_q1(ix2) meta_tbl.frac_mark_cells_q2(ix2) ...
    meta_tbl.frac_mark_cells_q3(ix2) meta_tbl.frac_mark_cells_q4(ix2) ...
    meta_tbl.frac_mark_cells_q5(ix2)], ...
    'Relative Rad. to Center','Fraction of RFP+ Cells','b')
