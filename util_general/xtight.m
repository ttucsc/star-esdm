function xtight(ax)
% sets x axis only to tight
if (~exist("ax","var")), ax = gca; end
  % Set axis tight only on y-axes
  yl=ylim(ax); % retrieve auto y-limits
  axis tight   % set tight range
  ylim(ax,yl)  % restore y limits 
end
