function coords = line_coords(pts, step_size, do_close, bounds, do_round)
% coords = line_coords(A, pts, step)
%   generates coordinates for points along lines specified by points given in n x 2 array of (x,y) pts
%   if (do_close is true, closes the polygon by drawing from last point to first point.
%
%   if do_round is true, coords are rounded to nearest integer. 
%   if bounds given [x1, y1; x2, y2] then all coordinates outside of bounds are removed.
%

    if (~exist('do_round', 'var')),   do_round = false;  end
    if (~exist('bounds',   'var')),   bounds   = [];     end
    if (~exist('do_close', 'var')),   do_close = false;  end
    if (~exist('step_size','var')),   step_size = 1;     end

    [nr,nc] = size(pts);
    if (isempty(pts))
        coords = [];
        return; 
    end
    if (nc~=2), error("error:  line_coords:  pts must be an n X 2 vector"); end
    if (nr==1)
        coords = pts;
        return;
    end
    npts = size(pts,1);
    msegs = npts-1;
    
    segments = cell(msegs+do_close,1);
    nsegs = 0;
    first_pt = [];
    lens = nan(msegs+do_close,1);
    for i=1:msegs
        if (any(isnan(pts(i:i+1,:)))), continue; end
        if (isempty(first_pt)), first_pt = i; end
        nsegs = nsegs + 1;
        segments{nsegs} = line_coords_sub(pts(i,:), pts(i+1,:), step_size, do_round);
        lastpt = i+1;
        lens(nsegs) = size(segments{nsegs},1);
    end
    if (isempty(first_pt))
        coords = [];
        return;
    end
    if (do_close)
        nsegs = nsegs + 1;
        segments{nsegs} = line_coords_sub(pts(lastpt,:), pts(first_pt,:), step_size, do_round);
        lens(nsegs) = size(segments{nsegs},1);        
    end

    nout = nansum(lens);
    
    coords = nan(nout,2);
    ix1=1;
    for i=1:nsegs
        ix2 = ix1 + lens(i)-1;
        coords(ix1:ix2,:) = segments{i};
        ix1 = ix1 + lens(i);
    end
    
        % remove any points outside of bounds
    if (~isempty(bounds))
        xkeepers = coords(:,1) >= bounds(1,1) & coords(:,1) <= bounds(2,1);
        ykeepers = coords(:,2) >= bounds(1,2) & coords(:,2) <= bounds(2,2);
        keepers = xkeepers & ykeepers;
        coords(~keepers,:) = nan;
    end
    
    coords = keep_unique(coords, do_round);
    
end

function coords = line_coords_sub(p1, p2, step_size, do_round)

    dxy = p2 - p1;
    len = sqrt(dxy(1)*dxy(1) + dxy(2)*dxy(2));
    npts = ceil(len / step_size) + 1;
    theta = atan2(dxy(2),dxy(1));
    l=(0:(npts-2))'*step_size;
    coords = [p1 + l*[cos(theta), sin(theta)]; ...
              p2];
    coords = keep_unique(coords, do_round);
end

function coords = keep_unique(coords, do_round)

        % remove any duplicates
    if (do_round)
        coords = round(coords);
    end
    xdifs = diff(coords(:,1));
    ydifs = diff(coords(:,2));
    dups = xdifs == 0 & ydifs == 0;
    coords(dups,:) = [];
    
        % remove any duplicate nan's
    xnans = isnan(coords(:,1));
    ynans = isnan(coords(:,2));
    mynans = (xnans | ynans);
    coords(mynans,:)=[];
%     dupnans = mynans(1:end-1) & mynans(2:end);
%     coords(dupnans) = [];
end
% 
%     if (isempty(a))
%         a = zeros(nr,nc);
%     else
%         [nr, nc] = size(a);
%     end
%     
%     m = -cot(theta);
%     b = rho/sin(theta);
%     
%     if (theta == 0)
%         c = round(rho);
%         if (c > 0 && c <= nc)
%             a(:,c)=val;
%         end
%     else
%         if (abs(tan(theta)) > 1)
%             for c=1:nc
%                 r = round(m*c + b);
%                 if (r>0 && r <= nr)
%                     a(r,c) = val;
%                 end
%             end
%         else
%             for r=1:nr
%                 c=round((r-b)/m);
%                 if (c>0 && c<=nc)
%                     a(r,c) = val;
%                 end
%             end
%         end
%     end
%     
% %     figure(98);
% %     aa = a;
% %     c = round(rho*cos(theta));
% %     r = round(rho*sin(theta));
% %     if (c>0 && c<=nc && r>0 && r <= nr)
% %         aa(r,c) = 2*val;
% %     end
% %     imshow(aa/max(max(aa)));
% %     fprintf('done\n');
% %         
% end
