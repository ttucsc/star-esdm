function xyztight(xtight, ytight, ztight, ax)
% sets only some axes to tight.
    if (~exist("ax","var")), ax = gca; end
    
    if (~xtight)
        xl = xlim(ax);
    end
    if (~ytight)
        yl = ylim(ax);
    end
    if (exist("ztight","var") && ~isempty(ztight) && ~ztight)
        zl = zlim(ax);
    else
        ztight = true;  % so we don't try to reset z limits.
    end
    axis tight
    if (~xtight), xlim(ax, xl); end
    if (~ytight), ylim(ax, yl); end
    if (~ztight), zlim(ax, zl); end

end
