function cmap = clim_cmap(len, typ, half, flag_minmax, do_green, fignum)
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
    if (~exist("do_green",   "var")), do_green    = false;      end 

    if (len <= 25)
        scaling = 20;
    else
        scaling = 15;
    end
    if (do_green), scaling = scaling/3; end
    is_odd = mod(len,2)==1;
    if (~isempty(half))
        len = len*2 - is_odd;
    end

    mid = ceil(len/2);
    band1 = zeros(len,1);
    
    if (strcmp(typ,"gaussian"))
    
        if (is_odd)
            c1 = gauss(len, len/6, mid)';
            c2 = gauss(len, len, mid)';
        else
            c1 = gauss(len-1, (len-1)/6, mid)';
            c1 = [c1(1:mid); c1(mid:end)];
            c2 = gauss(len-1, len-1, mid)';
            c2 = [c2(1:mid); c2(mid:end)];
        end
        band2 = c2 - min(c2);
        band2 = band2 / max(band2);
        c1 = c1 / max(c1);
        
        band1(1:mid)=c1(1:mid);
        band1(mid+1:end) = c1(mid+1:end).^(1/scaling);

        band3=flipud(band1);
    else    
        p1 = ceil(len/4);
        l1 = p1;
        l2 = mid-p1+1;
        p3 = ceil(3*len/4);
        l3 = p3-mid;
        l4 = len - p3+1;
        
        if (len <= 25)
            endval = 0.8;
            midval1 = endval + (1-endval)*2/3;
            midval2 = .875;
        else
            endval = 0.6;
            midval1 = endval + (1-endval)*2/3;
            midval2 = .75;
        end
        if (do_green)
            endval = endval/2; 
            midval2 = midval2/1.1;
        end
        
        
        band1(1:mid) = linspace(0,1,mid);
        band1(mid+1:p3) = linspace(1, midval1, l3);
        band1(p3:end) = linspace(midval1,endval, l4);
        band3=flipud(band1);
        
        band2 = zeros(len,1);
        band2(1:p1) = linspace(0,midval2, l1);
        band2(p1:mid) = linspace(midval2, 1, l2);
        if (is_odd)
            band2(mid:end) = flipud(band2(1:mid));
        else
            band2(mid+1:end) = flipud(band2(1:mid));
        end
%        g = [linspace(0,.75,p1)'; linspace(.75,1,l2)'; linspace(1,.75,l3)'; linspace(.75,0,l4)'];
%       g=[r(1:mid);b(mid+1:end)];
    end
    cmap = [band1,band2,band3];
    
    if (flag_minmax)
        cmap(1,:) = [0,.75,0];
        cmap(end,:) = [1, 0, 1];
    end

    if (~isempty(half))
        if (any(strcmp(half,["high","top"])))
            if (is_odd)
                cmap = cmap(mid:end,:);
            else
                cmap = cmap(mid+1:end,:);
            end
        elseif (any(strcmp(half,["low","bottom"])))
            cmap = flipud(cmap(1:mid,:));
        else
            error("input 'half' ('%s') must be 'high' (or 'top') or 'low' (or 'bottom')", half);
        end
    end
    
    if (do_green)
        cmap(1:mid,:) = [cmap(1:mid,1), cmap(1:mid,3), cmap(1:mid,2)];
        cmap = flipud(cmap);
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

