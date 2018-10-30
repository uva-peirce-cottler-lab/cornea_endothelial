function  plot_2sample(y1,y2,y_text,x_groups)

plot(1+(rand(size(y1,1),1)-.5)/10,y1,'ko', 'MarkerSize',4); hold on
add_errorbar(gcf, 1, y1, 0.4,'ABSOLUTE_WIDTH',1)
plot(2+(rand(size(y2,1),1)-.5)/10,y2,'ko', 'MarkerSize',4)
add_errorbar(gcf, 2, y2, 0.4,'ABSOLUTE_WIDTH',1)
hold off
ylabel(y_text)
set(gca,'xtick',[1:2],'xticklabel',x_groups)
set(gca,'Box', 'off', 'TickDir' , 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'FontName', 'Arial', 'YMinorTick', 'on', ...
    'YGrid', 'on', 'XGrid', 'off', 'XColor', [.3 .3 .3], ...
    'YColor', [.3 .3 .3],  'LineWidth', 1);
fig_pos = get(gcf,'position');
set(gcf,'position',[fig_pos(1:2) 250 250])
axis([.5 2.5 ylim])
% plot(repmat(1:5,[5 1]),
[h, p] = ttest2(y1,y2);

end

