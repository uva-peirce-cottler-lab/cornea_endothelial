function plot_cell_count(y1,y2, ytxt,x_cell)
hf = figure('Units','Pixels');
fig_pos = [500 500 250 220];

plot(1+(rand(size(y1,1),1)-.5)/10,y1,'ro', 'MarkerSize',4); hold on
add_errorbar(gcf, 1, y1, 0.4,'ABSOLUTE_WIDTH',1)
plot(2+(rand(size(y2,1),1)-.5)/10,y2,'bo', 'MarkerSize',4)
add_errorbar(gcf, 2, y2, 0.4,'ABSOLUTE_WIDTH',1)
hold off
ylabel(ytxt)
set(gca,'xtick',[1:2],'xticklabel',x_cell)
set(gca,'Box', 'off', 'TickDir' , 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'FontName', 'Arial', 'YMinorTick', 'on', ...
    'YGrid', 'on', 'XGrid', 'off', 'XColor', [.3 .3 .3], ...
    'YColor', [.3 .3 .3],  'LineWidth', 1);
set(gcf,'position',fig_pos)
axis([.5 2.5 ylim])
% plot(repmat(1:5,[5 1]),
[h, p] = ttest2(y1,y2);
fprintf('p: %.4f\n',p);


end

