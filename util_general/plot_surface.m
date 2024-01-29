function [h, vpos] = plot_surface(fignum, valsurf, x, y, do_probs, subfig, ttl, do_clear, xlbl, ylbl, zlbl, legs, viewpos, do_lpf, xx, yy, anoms)
%
%   do_probs:  
%       1/true  : plot probs .001 to .999
%       2       : plot probs .0001 to .9999
%       3,>3    : plot probs .00001 to .99999
%       vector  : use vector values as line indexes.  draws fixed lines.

    if (~exist('do_lpf','var') || isempty(do_lpf)), do_lpf=false; end
    
    probs = [ .000001 .00001 .0001 .001  .01 .0228  .05  .10 .1587 .50 .8413  .90  .95  .9772 .99 .999   .9999, .99999, .999999];
    pclr    = [  'c',    'c', 'c',  'y', 'c', 'g',   'c', 'c', 'g',  'r', 'g', 'c', 'c', 'g', 'c', 'y',   'c',    'c',    'c'];
    pwid    = [   1,      1,   1,    2,  1.5,  2,   1.5,  1.5,  2,    2,   2,  1.5, 1.5, 1.5,  2,   2,     1,      1,      1 ];
    % .025, .975 are aprox 2 sigma in a normal distribution
    % .1587 and .8413 are +/- 1 sigma. in a normal distribution
    
    my_lines=false;
    if (length(do_probs)>1)
        nprobs = length(do_probs);
        if (mod(nprobs,2)==1)
            ctr = ceil(nprobs/2);
        else
            ctr = [floor((nprobs-1)/2),ceil((nprobs-1)/2)];
        end
        yvals = do_probs;
        do_probs = true;
        my_lines = true;
        pclr  = repmat('c',1,nprobs);
        pclr(ctr) = 'r';
        pwid = repmat(2,1,nprobs);
    elseif (do_probs==2)
        probs = probs(2:end-1);
        pclr  = pclr(2:end-1);
        pwid  = pwid(2:end-1);
    elseif (do_probs==1)
        probs = probs(3:end-2);
        pclr  = pclr(3:end-2);
        pwid  = pwid(3:end-2);
    end
    h = figure(fignum);
    if (exist('do_clear','var') && ~isempty(do_clear) && do_clear), clf; end
    if (exist('subfig', 'var') && ~isempty(subfig)),  h = subplot(subfig(1),subfig(2),subfig(3)); end

        if (exist('x','var') && ~isempty(x))
            if (size(valsurf,2)==1)
                plot(y,valsurf);
            else
                surf(x, y, valsurf,'edgecolor','none');
            end
        else
            if (size(valsurf,2)==1)
                plot(valsurf);
            else
                surf(valsurf,'edgecolor','none');
            end
        end
    lights_on; 
    grid on;
    
    if (exist('anoms','var') && ~isempty(anoms))
                % draw the actual anomaly values. 
        hold on;
        if (size(anoms,2)==1)
            scatter(yy, anoms, 10, 'markerfacecolor',[1,.5,.5], 'markeredgecolor',[1,.5,.5]);
        else
            scatter3(xx,yy, anoms, 10, 'markerfacecolor',[1,.5,.5], 'markeredgecolor',[1,.5,.5]);
        end
        hold off;
    end

    if (size(valsurf,2)~=1)
        if (do_probs)
            if (valsurf(end,1) > .5)
                is_cdf = true;
            else
                is_cdf = false;
            end
            hold on;
            if (my_lines)
                [problines, pdf_xlines] = find_lines(valsurf, yvals, y);
                nprobs = length(yvals);
            else
                if (do_lpf)
                    if (is_cdf)
                        [problines,  pdf_xlines] = find_problines([], valsurf,  probs, y, false, 1, 1e-14, false, 4.0, .5);
                    else
                        [problines,  pdf_xlines] = find_problines(valsurf,  [], probs, y, false, 1, 1e-14, false, 4.0, .5);
                    end
                else
                    if (is_cdf)
                        [problines,  pdf_xlines] = find_problines([], valsurf,  probs, y, false, 1, 1e-14, false, [], []);
                    else
                        [problines,  pdf_xlines] = find_problines(valsurf,  [], probs, y, false, 1, 1e-14, false, [], []);
                    end
                end
                nprobs = length(probs);
            end
            xlen = length(x);
            for i=1:nprobs
                plot3((1:xlen), problines(:,i), pdf_xlines(:,i), pclr(i), 'LineWidth',pwid(i));
            end
            hold off;
        end
    end
    
%    xyztight(true, true, false); 
    set(gca, 'XLimSpec', 'Tight');
    set(gca, 'ZLimSpec', 'Tight');
    
    if (exist('ttl',    'var') && ~isempty(ttl)),      title(ttl,'interpreter','none');    end
    if (exist('xlbl',   'var') && ~isempty(xlbl)),     xlabel(xlbl,'interpreter','none');  end
    if (exist('ylbl',   'var') && ~isempty(ylbl)),     ylabel(ylbl,'interpreter','none');  end
    if (exist('zlbl',   'var') && ~isempty(zlbl)),     zlabel(zlbl,'interpreter','none');  end
    if (exist('legs',   'var') && ~isempty(legs)),     legend(legs,'interpreter','none');  end
    if (exist('viewpos','var') && ~isempty(viewpos)),  view(viewpos); end

    vpos = view();

end

function [xvals, zvals] = find_lines(valsurf, yvals, y)

    [~,nc] = size(valsurf);
    nlocs = length(yvals);
    ylocs=nan(1,nlocs);
    for i=1:nlocs
        p = find(y==yvals(i),1); 
        if (~isempty(p)), ylocs(i) = p; end
    end
    ylocs(isnan(ylocs))=[];
    xvals = repmat(yvals,nc,1);
    zvals = valsurf(ylocs,:)';
end

