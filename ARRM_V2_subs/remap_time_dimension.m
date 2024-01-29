function [hc_out, dayPos] = remap_time_dimension(hc, out_yrlen, distrib, dayPos)
% Inputs:
%       hc:  histogram counts, size nbins x yrlen.
%       dayPos:  nsets X yrlen, or 1 (see below).
%   pass in dayPos for mapping back to original timespace.  (with distrib set to empty)
%   pass in distribution (precip climatology, or cos(climatology) for temp) to histogram-equalize to new timespace
%   pass in distrib as 1 to remap with equal weights.
%
%   NOTE:  if using on PDFs or CDFs, be sure to renormalize after calling this function.

    nbins = size(hc,1);
    yrSize = size(hc, 2);
    nsteps = size(hc,3);
    
    if (~exist('dayPos','var') || isempty(dayPos))
        if (~exist('distrib','var') || isempty(distrib))
            distrib = ones(nsteps,yrSize);
        elseif (length(distrib)==1)
            distrib = repmat(distrib, nsteps,yrSize); 
        end
    end
    
    if (exist('dayPos','var') && ~isempty(dayPos))        % mapping back to original time frame.  Just interpolate    
                                                          % hc is now yrlen X nbins
        if (nsteps == 1)
            hc=hc';                                           % transpose histcounts matrix, because interp1 can only work on a matrix down the columns
            hc2 = [hc(end,:);hc;hc(1,:)];                                   % histograms are circular, and we need to be able to map around the year end.
                                                                            % dayPos can be will be between 0 & 1 or between yrSize and yrSize+1.
            hc_out = (interp1((0:yrSize+1)', hc2, dayPos))';                % transpose back so hc_out is size nbins X yrlen X nsets
        else                                                                % NOTE:  interp1 can only work on a matrix if you use linear interpolation.  Otherwise you have to do it with a loop.
            hc_out = zeros(nbins, out_yrlen, nsteps);
            for i=1:nsteps
%                hc2 = [hc(end,:,i);hc(:,:,i);hc(1,:,i)];
                hc2 = [hc(:,end,i), hc(:,:,i),hc(:,1,i)]';
                hc_out(:,:,i) = (interp1((0:yrSize+1)', hc2, dayPos(i,:)))';      % transpose back so hc_out is size nbins X yrlen X nsets      
            end
        end    
                                            % mapping to a smaller time frame.  interpolate each original day to its
                                            % mapped location and sum.
    else
        try
            [hc_out, dayPos] = histogram_equalize(hc, distrib, out_yrlen, yrSize);
        catch
            fprintf(2, "oops.  remat_time_dimension(...)\n");
        end
    end
end

function [hcEq, dayPos] = histogram_equalize(hc, distrib, out_yrlen, yrlen)
    nsteps = size(hc,3);
    if (nsteps == 1)
        if (out_yrlen == 1)
            hcEq = sum(hc,2);
            dayPos = ones(yrlen);
        else
            [hcEq,  dayPos] = hist_equalize_hcounts(hc,  distrib, out_yrlen, yrlen);   
        end
    else
        nbins = size(hc,1);
        hcEq = zeros(nbins,out_yrlen, nsteps); 
        dayPos = zeros(nsteps,yrlen);
        for i=1:nsteps
            [hcEq(:,:,i), dayPos(i,:)] = hist_equalize_hcounts(hc(:,:,i), distrib(i,:), out_yrlen, yrlen);
        end
    end
end

function [hcEq, dayPos] = hist_equalize_hcounts(hCounts, distrib, out_yrlen, yrlen)
%   hCounts         histogram, size nbins x yrlen
%   clim            climatology (# rainfall events/day)
%   yrSize          size to equalize hCounts to
%
%   hcEq            equalized hCounts
%   dayPos          location where each original day is located in hcEq

    nbins = size(hCounts, 1);
    distrib = max(.0001,distrib);     % in case Climatology goes negative.        % maybe s/b smaller, ian?
    hcEq   = zeros(nbins, out_yrlen);      % points interpolated.
    cumpdf = cumsum(distrib,'omitnan')/nansum(distrib);
    dayPos=zeros(1,yrlen);
    for i=1:yrlen
        x=out_yrlen*cumpdf(i);
        dayPos(i) = x;
        dx = mod(x,1.0);
        xi = floor(x);
        xj = xi+1;
        xi = mod(xi-1+out_yrlen,out_yrlen)+1;         % and adjust them so if they are between last day and 1st day,
        xj = mod(xj-1+out_yrlen,out_yrlen)+1;         % we wrap them at the end.
        try
            hcEq(:,xi)  = hcEq(:,xi) + hCounts(:,i) * (1-dx);
            hcEq(:,xj)  = hcEq(:,xj) + hCounts(:,i) * dx;
        catch
            fprintf('day-mapping-oops!\n');
        end
    end
end


