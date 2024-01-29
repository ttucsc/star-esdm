function fig_title(h, ttl)
% fig_title(h, ttl)         Places a title centered at the top of a figure
%                           useful for figures with multiple plots.
    ha = axes(h, 'Position',[0 0 1 1],'Visible','off');
    text(ha,.5,.975,ttl, "fontsize", 16,"HorizontalAlignment","center", "fontweight","bold");
end

