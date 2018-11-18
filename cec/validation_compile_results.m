function validation_compile_results()
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
csv_path = [out_path '/img_tbl.csv'];
tbl = readtable(csv_path);

% Lao dmanual coutn results and add to tabl
xlm_names = dir([out_path '/CellCounter_*.xml']);
for n=1:numel(xlm_names)
   counts = import_cellcounter_xmls([out_path '/' xlm_names(n).name]);

   % Remove cell counter string, find idex in tbl, add entries
   base_name = regexp(xlm_names(n).name,'CellCounter_(.*)_valid','tokens','once');
   
   ind = strmatch(base_name,tbl.img_name);
    tbl.manual_count_rfp_plus(ind) = counts(2);
    tbl.manual_count_rfp_minus(ind) = counts(4);
end

results_tbl=tbl(tbl.manual_count_rfp_plus~=0,:);


nsamples = numel(xlm_names);

% Plot paired unlabeled
plot(ones(nsamples,1),results_tbl.manual_count_rfp_minus, 'bo','Markersize',5); hold on
plot(2*ones(nsamples,1),results_tbl.auto_count_rfp_minus, 'bo','Markersize',5)
for n=1:nsamples
   plot([1 2], [results_tbl.manual_count_rfp_minus(n) ...
       results_tbl.auto_count_rfp_minus(n)],'Color', [.1 .1 .1])
end
xr = xlim; axis([0.5 2.5 ylim]); xticks([1 2])
ylabel('RFP- Cell Nuclei Count'); set(gcf,'position', [100 100 175 170])
xticklabels({'Manual','Auto'})

% Plot paired unlabeled
plot(ones(nsamples,1),results_tbl.manual_count_rfp_plus, 'ro','Markersize',5); hold on
plot(2*ones(nsamples,1),results_tbl.auto_count_rfp_plus, 'ro','Markersize',5)
for n=1:nsamples
   plot([1 2], [results_tbl.manual_count_rfp_plus(n) ...
       results_tbl.auto_count_rfp_plus(n)],'Color', [.1 .1 .1])
end
xr = xlim; axis([0.5 2.5 ylim]); xticks([1 2])
ylabel('RFP+ Cell Nuclei Count'); set(gcf,'position', [100 100 175 170])
xticklabels({'Manual','Auto'})



BlandAltman(results_tbl.manual_count_rfp_minus, results_tbl.auto_count_rfp_minus,...
    'RFP- Cellcount')
set(gcf,'position',[2.9951    7.1332   16.6053    5.2632])

BlandAltman(results_tbl.manual_count_rfp_plus, results_tbl.auto_count_rfp_plus,...
    'RFP+ Cellcount')
axis([xlim -5 5])
set(gcf,'position',[2.9951    7.1332   16.6053    5.2632])
hv = get(gca,'children');
set(hv(5:8),'MarkerEdgeColor',[1 0 0])
% set(
end

