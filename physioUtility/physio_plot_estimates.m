function physio_plot_estimates(b,cols,names,colors,varargin)

npanels = length(cols) + 1;
tor_fig(1,npanels);

dorobust = 0;
if length(varargin) > 0, dorobust = varargin{1};,end

indx = 1;
for i = cols    % which columns

    tmp = squeeze(b(:,i,:))';

    subplot(1,npanels,indx);
    plot(tmp')
    title(names{i})

    subplot(1,npanels,npanels);
    h(indx) = tor_fill_steplot(tmp,colors{indx},dorobust);

    title('Group averages')
    drawnow

    indx = indx + 1;
end

legend(h,names(cols),7);

return