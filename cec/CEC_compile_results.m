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


save([in_path '/CEC_results.mat']);

load([in_path '/CEC_results.mat']);
keyboard

fig_pos = [500 500 250 200];

% IP Induction Experiment
hf = figure('Units','Pixels');
ix1=strcmp(meta_tbl.tx_type,'IP') & ...
    strcmp(meta_tbl.term,'Long') & meta_tbl.day==28;
ix2=strcmp(meta_tbl.tx_type,'IP') & ...
    strcmp(meta_tbl.term,'Long') & meta_tbl.day==112;
y1=meta_tbl.tot_cells(ix1);
y2=meta_tbl.tot_cells(ix2);
plot(1+(rand(size(y1,1),1)-.5)/10,y1,'ro', 'MarkerSize',4); hold on
add_errorbar(gcf, 1, y1, 0.4,'ABSOLUTE_WIDTH',1)
plot(2+(rand(size(y2,1),1)-.5)/10,y2,'bo', 'MarkerSize',4)
add_errorbar(gcf, 2, y2, 0.4,'ABSOLUTE_WIDTH',1)
hold off
ylabel('Total Cells')
set(gca,'xtick',[1:2],'xticklabel',{'6','16'})
set(gca,'Box', 'off', 'TickDir' , 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'FontName', 'Arial', 'YMinorTick', 'on', ...
    'YGrid', 'on', 'XGrid', 'off', 'XColor', [.3 .3 .3], ...
    'YColor', [.3 .3 .3],  'LineWidth', 1);
set(gcf,'position',fig_pos)
axis([.5 2.5 ylim])
% plot(repmat(1:5,[5 1]),
[h, p] = ttest2(y1,y2);

hf = figure('Units','Pixels');
y1 = [meta_tbl.frac_mark_cells_q1(ix1) meta_tbl.frac_mark_cells_q2(ix1) ...
    meta_tbl.frac_mark_cells_q3(ix1) meta_tbl.frac_mark_cells_q4(ix1) ...
    meta_tbl.frac_mark_cells_q5(ix1)];
y2 = [meta_tbl.frac_mark_cells_q1(ix2) meta_tbl.frac_mark_cells_q2(ix2) ...
    meta_tbl.frac_mark_cells_q3(ix2) meta_tbl.frac_mark_cells_q4(ix2) ...
    meta_tbl.frac_mark_cells_q5(ix2)];
errorbar(0:.25:1,mean(y1),std(y1),'r-') 
hold on
errorbar(0:.25:1,mean(y2),std(y2),'b-') 
ylabel('Fraction of TOM+ Cells')
xlabel('Relative Rad. to Center')
% set(gca,'xtick',[1:5],'xticklabel',{'0','16w'})
set(gca,'Box', 'off', 'TickDir' , 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'FontName', 'Arial', 'YMinorTick', 'on', ...
    'YGrid', 'on', 'XGrid', 'off', 'XColor', [.3 .3 .3], ...
    'YColor', [.3 .3 .3],  'LineWidth', 1);
set(gcf,'position',fig_pos)
axis([-.1 1.1 ylim])
[p,tbl,stats] =anova2(vertcat(y1, y2(1:size(y1,1),:)),8);
multcompare(stats)



% Eyedrop Induction short Experiment
hf = figure('Units','Pixels');
ix1 = strcmp(meta_tbl.tx_type,'EyeDrop') & ...
    strcmp(meta_tbl.term,'Short') & meta_tbl.day==24;
ix2 = strcmp(meta_tbl.tx_type,'EyeDrop') & ...
    strcmp(meta_tbl.term,'Short') & meta_tbl.day==48;
y1=meta_tbl.tot_cells(ix1);
y2=meta_tbl.tot_cells(ix2);
plot(1+(rand(size(y1,1),1)-.5)/5,y1,'ro'); hold on
add_errorbar(gcf, 1, y1, 0.4,'ABSOLUTE_WIDTH',1)
plot(2+(rand(size(y2,1),1)-.5)/5,y2,'bo')
add_errorbar(gcf, 2, y2, 0.4,'ABSOLUTE_WIDTH',1)
hold off
ylabel('Total Cells')
set(gca,'xtick',[1:2],'xticklabel',{'24','48'})
set(gca,'Box', 'off', 'TickDir' , 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'FontName', 'Arial', 'YMinorTick', 'on', ...
    'YGrid', 'on', 'XGrid', 'off', 'XColor', [.3 .3 .3], ...
    'YColor', [.3 .3 .3],  'LineWidth', 1);
set(gcf,'position',fig_pos)
axis([.5 2.5 ylim])
[h, p] = ttest2(y1,y2);

hf = figure('Units','Pixels');
y1 = [meta_tbl.frac_mark_cells_q1(ix1) meta_tbl.frac_mark_cells_q2(ix1) ...
    meta_tbl.frac_mark_cells_q3(ix1) meta_tbl.frac_mark_cells_q4(ix1) ...
    meta_tbl.frac_mark_cells_q5(ix1)];
y2 = [meta_tbl.frac_mark_cells_q1(ix2) meta_tbl.frac_mark_cells_q2(ix2) ...
    meta_tbl.frac_mark_cells_q3(ix2) meta_tbl.frac_mark_cells_q4(ix2) ...
    meta_tbl.frac_mark_cells_q5(ix2)];
errorbar(0:.25:1,mean(y1),std(y1),'r') 
hold on
errorbar(0:.25:1,mean(y2),std(y2),'b') 
ylabel('Fraction of TOM+ Cells')
xlabel('Relative Rad. to Center')
% set(gca,'xtick',[1:5],'xticklabel',{'0','16w'})
set(gca,'Box', 'off', 'TickDir' , 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'FontName', 'Arial', 'YMinorTick', 'on', ...
    'YGrid', 'on', 'XGrid', 'off', 'XColor', [.3 .3 .3], ...
    'YColor', [.3 .3 .3],  'LineWidth', 1);
set(gcf,'position',fig_pos)
axis([-.1 1.1 ylim])
[p,tbl,stats] =anova2(vertcat(y1(1:size(y2,1),:), y2),size(y2,1));
multcompare(stats)


% Eyedrop Induction mid Experiment
hf = figure('Units','Pixels');
ix1 = strcmp(meta_tbl.tx_type,'EyeDrop') & ...
    strcmp(meta_tbl.term,'Mid') & meta_tbl.day==2;
ix2 = strcmp(meta_tbl.tx_type,'EyeDrop') & ...
    strcmp(meta_tbl.term,'Mid') & meta_tbl.day==21;
ix3 = strcmp(meta_tbl.tx_type,'EyeDrop') & ...
    strcmp(meta_tbl.term,'Mid') & meta_tbl.day==42;

y1=meta_tbl.tot_cells(ix1);
y2=meta_tbl.tot_cells(ix2);
y3=meta_tbl.tot_cells(ix3);
plot(1+(rand(size(y1,1),1)-.5)/5,y1,'ro'); hold on
add_errorbar(gcf, 1, y1, 0.4,'ABSOLUTE_WIDTH',1)
plot(2+(rand(size(y2,1),1)-.5)/5,y2,'bo')
add_errorbar(gcf, 2, y2, 0.4,'ABSOLUTE_WIDTH',1)
plot(3+(rand(size(y3,1),1)-.5)/5,y3,'bo')
add_errorbar(gcf, 3, y3, 0.4,'ABSOLUTE_WIDTH',1)
hold off
ylabel('Total Cells')
set(gca,'xtick',[1:3],'xticklabel',{'2','21','42'})
set(gca,'Box', 'off', 'TickDir' , 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'FontName', 'Arial', 'YMinorTick', 'on', ...
    'YGrid', 'on', 'XGrid', 'off', 'XColor', [.3 .3 .3], ...
    'YColor', [.3 .3 .3],  'LineWidth', 1);
set(gcf,'position',fig_pos)
axis([.5 3.5 ylim])
% [h, p] = anova(y1,y2);

hf = figure('Units','Pixels');
y1 = [meta_tbl.frac_mark_cells_q1(ix1) meta_tbl.frac_mark_cells_q2(ix1) ...
    meta_tbl.frac_mark_cells_q3(ix1) meta_tbl.frac_mark_cells_q4(ix1) ...
    meta_tbl.frac_mark_cells_q5(ix1)];
y2 = [meta_tbl.frac_mark_cells_q1(ix2) meta_tbl.frac_mark_cells_q2(ix2) ...
    meta_tbl.frac_mark_cells_q3(ix2) meta_tbl.frac_mark_cells_q4(ix2) ...
    meta_tbl.frac_mark_cells_q5(ix2)];
y3 = [meta_tbl.frac_mark_cells_q1(ix3) meta_tbl.frac_mark_cells_q2(ix3) ...
    meta_tbl.frac_mark_cells_q3(ix3) meta_tbl.frac_mark_cells_q4(ix3) ...
    meta_tbl.frac_mark_cells_q5(ix3)];
errorbar(0:.25:1,mean(y1),std(y1),'r') 
hold on
errorbar(0:.25:1,mean(y2),std(y2),'b') 
errorbar(0:.25:1,mean(y3),std(y3),'g') 
ylabel('Fraction of TOM+ Cells')
xlabel('Relative Rad. to Center')
% set(gca,'xtick',[1:5],'xticklabel',{'0','16w'})
set(gca,'Box', 'off', 'TickDir' , 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'FontName', 'Arial', 'YMinorTick', 'on', ...
    'YGrid', 'on', 'XGrid', 'off', 'XColor', [.3 .3 .3], ...
    'YColor', [.3 .3 .3],  'LineWidth', 1);
set(gcf,'position',fig_pos)
axis([-.1 1.1 ylim])
nsamples = min([size(y1,1) size(y2,1)  size(y3,1)]);
[p,tbl,stats] =anova2(vertcat(y1(1:nsamples,:),y2(1:nsamples,:), y3(1:nsamples,:)),nsamples);
multcompare(stats)




%    % Find center of gravity of cornea
%    st=regionprops(bw_cornea,'Centroid','ConvexImage');
%    cornea_rc = round([st.Centroid(2) st.Centroid(1)]);
%    
%    
%    % Make radial distance image from cornea center, display coordinates
%    % only with mask foreground pixels
%    bw_temp = false(size(bw_cornea));
%    bw_temp(cornea_rc(1),cornea_rc(2))=1;
%    ed_ctr = bwdist(bw_temp) .* bw_cornea;
%    max_dist = max(ed_ctr(:));
%    % Relative ED image from edge
%    rel_ed_edge = (1-ed_ctr./max_dist) .* bw_cornea;
%    
%    % EU of cornea from border
% %    ed_edge = bwdist(~bw_cornea);
%    
%    % 1-EU is distance from center
% %    ed_ctr_rel = ed_edge./max(ed_edge(:));
%    
%    % Sample distance from each 
%    lin_coords = round(sub2ind(size(bw_cornea),coords_tbl.X(:,1),coords_tbl.Y));
%    
%    % Get relative distance for each cell from edge to center
%    cell_ctr_reldist{n} = rel_ed_edge(lin_coords);
% %    cell_ctr_reldist{n} = rand(1,300);
%    
%    % Binn relative radial distances with histogram
%    [binned_cell_ctr_reldist(n,:),edges] = histcounts(cell_ctr_reldist{n},10);
%   
%    % Calculate Actual Area of each Bin
%    binned_area(n,:) = histcounts(rel_ed_edge, [1e-6 .1:.1:1]);
   
   
% end
% 
% keyboard


% all_cell_dist = vertcat(cell_ctr_reldist{:});
% plot(all_cell_dist,'.')

% binned_area = -diff(pi*(1-edges).^2);

anm_binned_cel_ctr_reldist = binned_cell_ctr_reldist ./ binned_area;

% anm_binned_cel_ctr_reldist = ...
%     bsxfun(@rdivide, double(binned_cell_ctr_reldist), binned_area);

% Plot cell count 
errorbar(edges(1:end-1)+diff(edges)/2,mean(binned_cell_ctr_reldist(tp==2,:),1),...
    std(binned_cell_ctr_reldist(tp==2,:),[],1),'r--');
hold on
errorbar(edges(1:end-1)+diff(edges)/2,mean(binned_cell_ctr_reldist(tp==21,:),1),...
    std(binned_cell_ctr_reldist(tp==21,:),[],1),'b--');
hold off
legend({'Day 2','Day 21'});
pos = get(gcf,'position');
set(gcf,'position', [pos(1:2) 400 275]);
[h,p] = ttest(mean(binned_cell_ctr_reldist(tp==2,:)),...
    mean(binned_cell_ctr_reldist(tp==21,:)));


% Plot cellcoutn normalized by area
errorbar(edges(1:end-1)+diff(edges)/2,mean(anm_binned_cel_ctr_reldist(tp==2,:),1),...
    std(anm_binned_cel_ctr_reldist(tp==2,:),[],1),'r--');
hold on
errorbar(edges(1:end-1)+diff(edges)/2,mean(anm_binned_cel_ctr_reldist(tp==21,:),1),...
    std(anm_binned_cel_ctr_reldist(tp==21,:),[],1),'b--');
hold off
legend({'Day 2','Day 21'});
pos = get(gcf,'position');
set(gcf,'position', [pos(1:2) 400 275]);

% Plot cell


xa=xlim; ya=ylim;
axis([xa(1) xa(2)+.05 ya]);


% end

