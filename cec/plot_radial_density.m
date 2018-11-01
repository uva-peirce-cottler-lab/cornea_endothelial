function plot_radial_density(x1,y1, xtxt,ytxt,color_letter)
hf = figure('Units','Pixels');
% y1 = [meta_tbl.frac_mark_cells_q1(ix1) meta_tbl.frac_mark_cells_q2(ix1) ...
%     meta_tbl.frac_mark_cells_q3(ix1) meta_tbl.frac_mark_cells_q4(ix1) ...
%     meta_tbl.frac_mark_cells_q5(ix1)];
% y2 = [meta_tbl.frac_mark_cells_q1(ix2) meta_tbl.frac_mark_cells_q2(ix2) ...
%     meta_tbl.frac_mark_cells_q3(ix2) meta_tbl.frac_mark_cells_q4(ix2) ...
%     meta_tbl.frac_mark_cells_q5(ix2)];
errorbar(x1,mean(y1),std(y1),[ color_letter '+']) 
set(gca,'Box', 'off', 'TickDir' , 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'FontName', 'Arial', 'YMinorTick', 'on', ...
    'YGrid', 'on', 'XGrid', 'off', 'XColor', [.3 .3 .3], ...
    'YColor', [.3 .3 .3],  'LineWidth', 1);
fig_pos = [500 500 250 220];
set(gcf,'position',fig_pos)
ylabel(ytxt)
xlabel(xtxt)
axis([-.1 1.1 ylim])
hold on
% keyboard
mdl = fitlm((x1)', mean(y1)');
% keyboard
fprintf('f(x) = %.3fx + %.3f\n',mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1));
mdl_coeffs = mdl.coefCI;
fprintf('[%.3f %.3f]\n',mdl_coeffs(2,2),mdl_coeffs(1,2))

x_pred = linspace(-.1,1.1,200)';
y_pred = predict(mdl,x_pred);

plot(x_pred,y_pred,'k','Linewidth',2)


end

