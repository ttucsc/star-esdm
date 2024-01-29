function cmap = temp_gradient_1(len, nsteps, colorfirst, colorlast, color0, colorend)
%   cmap = temp_gradient_1(len, nsteps, colorfirst, colorlast, color0, colorend)
%
%   Generate a color-gradient for temperature data, going from white to yellow to red to purple,
%   with optionally replaceing the first and/or last elements, 
%   and optionally prepending and/or appending one or more colors
%
%   colorfirst replaces 1st  element(s) of main colormap
%   colorlast  replaces last element(s) of main colormap
%   color0     is prepended to the colormap
%   colorend   is appended  to the colormap
%
%   total length is len + size(color0,1) + size(colorend,1)

    if (~exist('colorfirst', 'var')), colorfirst = []; end
    if (~exist('colorlast',  'var')), colorlast  = []; end
    if (~exist('color0',     'var')), color0     = []; end
    if (~exist('colorend',   'var')), colorend   = []; end
    clenfirst = size(colorfirst,1);
    clenlast  = size(colorlast, 1)-1;

    cm1 = flipud(hot(400));
    cm2 = flipud(hsv(400));
    cm1_steps = round(nsteps/2.25);
    cm2_steps = nsteps - cm1_steps;
    cm1_ix = round(linspace(20,200,cm1_steps));
    cmap1 = cm1(cm1_ix,:);
    cm2_ix = round(linspace(11,81,cm2_steps));
    cmap2 = cm2(cm2_ix,:);
    cmap3 = [cmap1;cmap2];
        
    cmstep = len/nsteps;
    
    cmap = ones(len,3);

    for ii = 1:nsteps
        ix1 = floor((ii-1)*cmstep) + 1;
        ix2 = floor( ii   *cmstep);
        for i = ix1:ix2
            cmap(i,:) = cmap3(ii,:);
        end
    end
    
    cmap(           1:clenfirst,:) = colorfirst;
    cmap(end-clenlast:end,      :) = colorlast; 
    cmap = [color0; cmap; colorend];
end
