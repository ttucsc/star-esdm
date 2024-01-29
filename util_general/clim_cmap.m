function cmap = clim_cmap(len, typ, half, flag_minmax, fignum)
% cmap = clim_cmap(len, typ, half, flag_minmax, fignum)
%
%   Generates a cold-to-hot colormap  (blue -> cyan -> white -> yellowish -> red)
%       two versions:  default uses gaussian shapes to colors;  linear uses straight lines.
%       add fignum to display the colors in a figure, along with a colorbar.
%
%   Inputs:
%       len         (optional) length of colormap  [256]
%                       if len is even, then 2 midpoints are white;  if odd, only midpoint is white.
%       typ         (optional) "gaussian" or "linear"  ["gaussian"]
%       half        (optional) "high", "low" or [], for top half or bottom half (flipped) of colormap, or entire colormap []
%       flag_minmax (optional)  If true, marks the highest red value with magenta, and the lowest blue value with green.
%                                   Useful to flag out-of-range values.
%       fignum      figure to draw color map and color bar in, For seeing what the colormap looks like.
%                       if fignum is given, then colormap is set to generated map.
%   
%   This needs an option to specify where white should be in the range.  Currently white is always in the center.

    if (~exist("len","var") || isempty(len))
        len = 256;
    end
    if (~exist('typ','var') || isempty(typ))
        typ = "gaussian";
    end
    if (~exist("half",       "var")), half        = strings(0); end
    if (~exist("fignum"     ,"var")), fignum      = [];         end 
    if (~exist("flag_minmax","var")), flag_minmax = false;      end 

    if (len <= 20)
        scaling = 20;
    else
        scaling = 15;
    end
    is_odd = mod(len,2)==1;
    if (~isempty(half))
        len = len*2 - is_odd;
    end

    mid = ceil(len/2);
    r = zeros(len,1);
    
    if (strcmp(typ,"gaussian"))
    
        if (is_odd)
            g1 = gauss(len, len/6, mid)';
            g2 = gauss(len, len, mid)';
        else
            g1 = gauss(len-1, (len-1)/6, mid)';
            g1 = [g1(1:mid); g1(mid:end)];
            g2 = gauss(len-1, len-1, mid)';
            g2 = [g2(1:mid); g2(mid:end)];
        end
        g = g2 - min(g2);
        g = g / max(g);
        g1 = g1 / max(g1);
        
        r(1:mid)=g1(1:mid);
        r(mid+1:end) = g1(mid+1:end).^(1/scaling);

        b=flipud(r);
    else        
        p1 = ceil(len/4);
        l2 = mid-p1+1;
        p3 = floor(3*len/4);
    %   l3 = p3 - mid+1;
        l4 = len - p3;
        
        r(1:mid) = linspace(0,1,mid);
        r(mid+1:p3) = 1;
        if (strcmp(typ,"log"))
            error("not written yet!");
        else
            r(p3+1:end) = linspace(1,.8,l4);
            b=flipud(r);
        
            g = zeros(len,1);
            g(1:p1) = linspace(0,.875, p1);
            g(p1:mid) = linspace(.875, 1, l2);
            if (is_odd)
                g(mid:end) = flipud(g(1:mid));
            else
                g(mid+1:end) = flipud(g(1:mid));
            end
        end
%        g = [linspace(0,.75,p1)'; linspace(.75,1,l2)'; linspace(1,.75,l3)'; linspace(.75,0,l4)'];
%       g=[r(1:mid);b(mid+1:end)];
    end
    cmap = [r,g,b];
    
    if (flag_minmax)
        cmap(1,:) = [0,.75,0];
        cmap(end,:) = [1, 0, 1];
    end

    if (~isempty(half))
        if (strcmp(half,"high"))
            if (is_odd)
                cmap = cmap(mid:end,:);
            else
                cmap = cmap(mid+1:end,:);
            end
        else
            cmap = flipud(cmap(1:mid,:));
        end
    end


    if (~isempty(fignum))
        figure(fignum);
        len = size(cmap,1);     % in case we're returning only half the regular map...
        colormap(cmap);
        if (strcmp(typ,"gaussian"))
            plot(1:len, cmap(:,1), 'r-', 1:len, cmap(:,3),"b-", 1:len,cmap(:,2), 'g-', "linewidth", 2);
        else
            plot(1:len, cmap(:,1), 'r-', 1:len, cmap(:,3),"b-", 1:len,cmap(:,2), 'g--', "linewidth", 2);
        end
        colorbar;
        drawnow();
    end
end

