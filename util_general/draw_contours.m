function img = draw_contours(img, M, levels, contour_vals)

% draws contours into matrix mat using contour info in M. where
% 
%   M = contourc(mat,levels);

%   vals are the values to set the contour lines 2.    

    [ny,nx] = size(img);
    len = size(M,2);
    nvals = length(contour_vals);
    istart = 1;
    while istart < len
        lvl = M(1,istart);
        ilvl = find(levels <= lvl,1,'last');
%         if (istart == 1) 
%             prev = ilvl; 
%         elseif (ilvl ~= prev)
%             figure(1);  surf(img,'edgecolor','none');   view([0,90])         
%             prev = ilvl;    
%         end
        if (isempty(ilvl)), ilvl = 1; end
        val = contour_vals(min(nvals, ilvl));
        npts = M(2,istart);
        istart = istart+1;
        ix1 = istart;
        ix2 = istart+npts-1;
        coords = line_coords(M(:,ix1:ix2)',.5,true, [1,1;nx,ny], true);
        ix = sub2ind(size(img), coords(:,2), coords(:,1));
        try
            img(ix) = val;                
        catch
            oops();
        end
        istart = istart + npts;
    end
end