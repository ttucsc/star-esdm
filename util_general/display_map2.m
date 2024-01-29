function [figh, subh, maph] = display_map2(img, ttl, R, fignum, subpos, ming, maxg, draw_coasts, cmap, lats_low, lons_low, figpos, fontsize, show_range, do_hide, figpause, nan_clr, sea_clr, keep_mask, nan_mask)
%        [figh, subh, maph] = display_map2(img, ttl, R, fignum, subpos, ming, maxg, draw_coasts, cmap, lats_low, lons_low, ...
%                                          figpos, fontsize, show_range, do_hide, figpause, nan_clr, sea_clr, keep_mask, nan_mask)
%
% draws a map in a figure, using the default map projection for the lat & lon range.
%   Uses matlab's worldmap(...) function to draw the map.
%
%   Inputs:
%       img         raster image of data to draw as map
%       ttl         title for map
%       R           Refmat object for map (R = make_refmat(lats, lons, nlats, nlons, dlat, dlon))
%                       or cell array with {lats, lons} marking center for each pixel in image.
%     optional:
%       fignum      figure to draw into, or handle of figure
%       subpos      if not empty, subplot position in figure
%       ming        min value for colormap
%       maxg        max value for colormap
%       draw_coasts if true, draw coastlines.  if == 2, draw N America coastlines and state boundaries.  if false or empty, 
%       cmap        desired colormap, or # of elements to use in standard CSC colormap, 
%                       which is based on jet(n), but has white at center (usually 0)
%       lats_low    for drawing a grid over the map.  Good to mark low-resolution grid spacing
%       lons_low        if empty, no low-res grid drawn on map.
%       figpos      figure position and size .  Will set fignum.Position = figpos.
%       fontsize    fontsize, in points [16]
%       show_range  t/f If true, adds data range to figure title [false]
%       do_hide     t/f If true, figure isn't shown.  Necessary when running w/o GUI, such as for batch jobs.    
%       figpause    pause time after drawing figure.  [0.1]; Set to longer for running interactively on HPCC
%       nan_clr     (background) color for portions of figure set to nan
%       mask_clr    (background) color for portions of the map masked by landmask.
%       keep_mask   binary mask; display all locations where keep_mask is true (optional)
%       nan_mask    binary mask; set all locations to nan where nan_mask is true.  (optional)
%       
%
%   Outputs
%       figh        figure handle
%       subh        handle to figure if drawn in a sub-position with subpos
%       maph        handle to map object within figh or subh.
%       
    if (~exist("keep_mask","var")), keep_mask = []; end
    if (~exist("nan_mask", "var")), nan_mask  = []; end
    img = squeeze(img);
    if (~isempty(keep_mask))
        myming = min(img(keep_mask>.1));
        mymaxg = max(img(keep_mask>.1));
    else
        myming = min(img(:));
        mymaxg = max(img(:));
    end

    if (~exist('ming','var') || isempty(ming)), ming = myming; end
    if (~exist('maxg','var') || isempty(maxg)), maxg = mymaxg; end
    if (~exist('fontsize','var') || isempty(fontsize)), fontsize = 16; end
    if (~exist('show_range','var') || isempty(show_range)), show_range = false; end
    if (~exist('do_hide','var') || isempty(do_hide)), do_hide = false; end
    if (~exist('figpause','var') || isempty(figpause)), figpause = 0.1; end
    if (~exist('sea_clr', 'var')), sea_clr = []; end
    
    [nlat,nlon] = size(img);
    cmap_len = size(cmap,1);
    IMG = nan(nlat,nlon,3);
    cmap_ix = (img-ming)/(maxg-ming)*(size(cmap,1)-1)+1;
    cmap_ix = max(1,min(cmap_len,round(cmap_ix)));

    for i=1:nlat
        for j=1:nlon
            if (~isempty(keep_mask) && keep_mask(i,j)<.5)
                IMG(i,j,:) = sea_clr;
            elseif (~isempty(nan_mask) && nan_mask(i,j)==true)
                IMG(i,j,:) = nan_clr;
            else
                try
                    IMG(i,j,:) = cmap(cmap_ix(i,j),:);
                catch
                    fprintf("oops:  %d %d  %.1f\n", i, j, cmap_ix(i,j))
                end
            end
        end
    end

        % fudge the corner pixels if we don't have the full range of pixel
        % values so we can label the colorbars properly.
    Imin = min(IMG(:));
    Imax = max(IMG(:));
    if (Imin > 0), IMG(1,1,:) = [0,0,0]; end
    if (Imax < 1), IMG(end,end,:) = [1,1,1]; end

    if (exist('fignum','var') && ~isempty(fignum))
        if (isa(fignum, "matlab.ui.Figure"))
            figh = fignum;
        else
            figh=figure(fignum);
        end
    else
        figh = gcf;
    end
    if (figh.Visible && do_hide)
        figh.Visible = false;
    end
    if (do_hide)
        pause(.05);
    end
    
    colorbar_below = true;
    if (exist('subpos','var') && ~isempty(subpos))
        if (iscell(subpos))
            subh = subplot(subpos{1}, subpos{2}, subpos{3}); 
            if (subpos{1} >= subpos{2})
                colorbar_below = false;
            end
        else
            subh = subplot(subpos(1), subpos(2), subpos(3));
            if (subpos(1) >= subpos(2))
                colorbar_below = false;
            end
        end
    else
        subh = figh;
    end
    
    maph = worldmap(zeros(size(img)), R);
    geoshow(maph, IMG, R);
    
%     maph = meshm(img,R);
%     set(maph, 'facealpha', 'texturemap', 'alphadata', double(~isnan(img)));
%     if (~isempty(bkgd_clr))
% %        set(maph, 'facealpha', 'texturemap', 'alphadata', ones(size(img)));
%         subh.Visible = true;
%         subh.Color = bkgd_clr;
%     end
% 
    
    if (exist('draw_coasts','var') && ~isempty(draw_coasts))
        if (draw_coasts)
            if (draw_coasts == 2)
                ncamerica_and_state_boundaries([],true);
            else
                load ('coastlines','coastlat','coastlon');
                plotm(coastlat, coastlon);
            end
        end
    end
            
    if (exist('lats_low','var') && ~isempty(lats_low))
        latbnd = R.LatitudeLimits;
        lonbnd = R.LongitudeLimits;
        ix1 = find(lats_low < latbnd(1),1,'last');
        if (isempty(ix1)), ix1=1; end
        ix2 = find(lats_low > latbnd(2),1);
        if (isempty(ix2)), ix2 = length(lats_low); end
        lats_low = lats_low(ix1:ix2);
%        lats_low(end+1) = nan;
        ix1 = find(lons_low < lonbnd(1),1,'last');
        if (isempty(ix1)), ix1 = 1; end
        ix2 = find(lons_low > lonbnd(2),1);
        if (isempty(ix2)), ix2 = length(lons_low); end
        lons_low = lons_low(ix1:ix2);
%        lons_low(end+1) = nan;

        
        for i=1:length(lats_low)
            plotm(lats_low(i)*ones(size(lons_low)), lons_low, '-','color',[.75,.75,.75]);
        end
        for i=1:length(lons_low)
            plotm(lats_low, lons_low(i)*ones(size(lats_low)), '-','color',[.75,.75,.75]);
        end
    end
    hold off;
%    if (exist('cmap','var') && ~isempty(cmap))
        colormap(subh, cmap);
%    end

%   imshow(flipud(img),[ming,maxg],'colormap',cmap);
    if (show_range) % (myming ~= ming || mymaxg ~= maxg)
        title(sprintf("%s range %.2f-%.2f", ttl, myming, mymaxg),'interpreter','none','fontsize',fontsize);
    else
        title(ttl,'interpreter','none','fontsize',fontsize);
    end

    minIMG = min(IMG(:));       % colorbar will pick up its range from the range of RGB values
    maxIMG = max(IMG(:));       % which, if the image is all white, will be very minimal.
    tick_marks = linspace(minIMG,maxIMG,11);
    tick_vals  = linspace(ming, maxg, 11);
    for i=1:length(tick_vals)
        if (maxg-ming < 1)
            tick_lbls(i) = sprintf("%.3f", tick_vals(i)); %#ok<AGROW> 
        elseif (maxg - ming <= 10)
            tick_lbls(i) = sprintf("%.1f", tick_vals(i)); %#ok<AGROW> 
        else
            tick_lbls(i) = sprintf("%d", round(tick_vals(i))); %#ok<AGROW> 
        end
            
    end
    if (colorbar_below)
        tbh = colorbar('location','southoutside');    % probably need tick marks and tick labels.
    else
        tbh = colorbar('location','eastoutside');
    end
    tbh.Ticks = tick_marks;
    tbh.TickLabels = tick_lbls;

    if (exist('figpos','var') && ~isempty(figpos))
        if (all(figpos<=3.0) && all(figpos >= -3))      % check fig position units.  +/-2 should be anywhere in the existing monitors of most systems for normalized units.
            figh.Units='normalized';
        else
            figh.Units='pixels';                          % otherwise, assume pixels.
        end
        figh.Position = figpos;
    end
%   drawnow();
    pause(figpause)

end



