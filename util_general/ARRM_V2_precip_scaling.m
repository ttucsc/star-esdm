function [prcp_out, spr_min, dx, sigma, kurt, edges, bins] =  ARRM_V2_precip_scaling(prcp_in, scaling, pr_min, direction, dx, nbins, nsigmas)
%   [prcp_out, spr_min, dx, edges, bins] =  ARRM_V2_precip_scaling(prcp_in, scaling, pr_min, direction, dx)
%   
%   Input:
%       prcp_in         regular precip, in mm. (for forward scaling), and scaled precip for reverse scaling.
%                         Given the same values on input  e.g. scaling = 1/2, pr_min = 0.1, then
%                           Forward scaling calculates prcp_in .^ scaling + some offset
%                           Reverse scaling calculates (prcp_in - same offset).^ 1/scaling and recovers the original
%                           values (except trace precip stays at zero on reverse).
%                           
%       scaling         power to use (for forward scaling);  reverse scaling will raise to 1/scaling, so supply the same
%                           scaling value for forward or reverse
%       pr_min          min precip value.  Precip below this will be mapped to zero.  Typical values: .1 or .254 for .1 mm or 1/100 inch 
%                           Trace precip (below pr_min) is mapped to center bin  on forward, and stays at zero on reverse.
%                           Provide same value (e.g. .1 or .254) for reverse scaling.
%       direction       "forward" or "reverse"
%       dx              for forward scaling:  scaled bin width to be used for histogramming.  If empty will be
%                           calculated based on nbins.  If present, will override nbins.
%                           (best practice is to leave blank, an use edges & bins provided on output.) 
%                       for reverse scaling;  must be the dx value used to generate the output to recover original values!
%       nbins           (optional)  # of evenly spaced edges for binning.  Default is 1000 bins.
%       nsigmas         (optional)  max # of std devs to bin over.  Should be at least 7 for ARRM_V2 to work.
%
%   Output
%       prcp_out        scaled (for forward) or unscaled (for reverse) precip values
%       spr_min         scaled value of pr_min
%       dx              bin & edge step size (needed to recover original values on reverse)
%               output useful only when doing forward scaling (are empty on return for reverse scaling.
%       sigma           std dev. of resulting rescaled precip when mirrored
%       kurt            kurtosis of resulting rescaled precip when mirrored.
%       edges           reasonable set of edges for ARRM_V2 to histogram properly.  1000 steps, out to +/- 10 sigma.
%       bins            mid-point bin values of edges.
%                           edges will straddle 0 at midpoint of edges array (index 501 of 1000)
%                           bins(501) will be the bin for dry days (precip less than pr_min)
%
%   Calculation done:
%           prcp_out = prcp_in.^(scaling) - pr_min^(scaling) + dx/2         forward  (note: spr_min = pr_min^(scaling) )
%           prcp_out = (prcp_in - dx/2 + pr_min^(scaling)).^(1/scaling)     reverse
%       
%   Forward scaling raises prcp_in to the scaling power.  So use 1/2 to take square root, 1/3 for cube root, etc.
%   Note that it shifts the data nearer to the origin, so a mirrored histogram will be centered @ 0, with dry days
%   counted in bin(501)  (==bin value 0, for precip < pr_min.
%
%   NOTE:  to safely histogram such that precip of pr_min is not binned in zero bin, edges must have a maximum delta of
%   dx.  Returned edges satisfy this requirement.
%
%   NOTE: to reverse the mapping, you MUST provide the original output's scaling, dx and pr_min.  And when mapping model
%   output onto observation's distribution, rescale the mapped output using the observation's scaling, dx and pr_min.
%
%   NOTE:  trace precip which maps to < dx/2 is NOT mapped to 0;  it WILL be binned to the zero-bin if the output is
%   histogrammed using the edges output.  precip which maps to below zero will be truncated at 0, and will stay at zero
%   on reverse scaling, but trace amounts that map to between 0 and dx/2 will map back to their original values, so set
%   these to zero explicitly if that's what you want.

    if (~exist("dx",     "var")),                     dx=[];        end
    if (~exist("nbins",  "var") || isempty(nbins)),   nbins = 1000; end
    if (~exist("nsigmas","var") || isempty(nsigmas)), nsigmas = 10; end
        

                % for ARRM_V2:  check we're not scaling the wrong way!
        if (scaling > 1), error("error:  forward scaling, but scaling term > 1:  %.4f", scaling); end
        
    if (strcmp(direction, "forward"))
    
        nanmap = isnan(prcp_in);
                % make sure there are no negative values, and scale in the input.
        prcp_in = max(0, prcp_in);      % this will set nan's to 0.  We'll set them back later.
        spr_min = pr_min ^ scaling;
        sprcp = prcp_in .^ scaling;
        
                % find range of scaled values, and create some temporary bins covering the data range.
        spr_max = nanmax(sprcp);        % max scaled precip value
        edges1 = linspace(spr_min, spr_max, 500);
        dx1 = edges1(2)-edges1(1);
            % get std dev of distrib, using temporrary set of bins just wide enough.
        hcnts = histcounts(sprcp, edges1);
        hcnts = [fliplr(hcnts),hcnts];
        bb = (edges1(1:end-1)+edges1(2:end))/2 - spr_min + dx1/2;
        bins1 = [fliplr(-bb), bb];
        [~,sigma,~,kurt] = pdf_stats(hcnts, bins1, [], true);
        
            % binning range for ARRM_V2
        if (~isempty(dx) && ~isempty(nbins))
            binmax = dx * nbins;
        else
            if (isempty(dx))
                dx = nsigmas*sigma/nbins;
            end
            binmax = max(3*spr_max, nsigmas*sigma);
                % adjust binmax and # of bins
            nbins = ceil(binmax/dx);
            binmax = nbins*dx;
        end
        
            % now get edges for full binning range
        edges = linspace(0, binmax, nbins+1)-dx/2;        % edges straddle 0, so original +/- dx go into bin(1).
        bins  = edges(1:(end-1))+dx/2;                         % and bin(2) is centered on pr_min's scaled value. 
        
            % shift the scaled data down so it is just above 0, so pr_min would be binned at the bottom (center?) of the 1st bin
            % above 0.  prcp values in the half-open range [0,prmin) will be binned in to the zero bin.
            % (maybe at center...check this, Ian!)
            
        prcp_out = sprcp - spr_min + dx/2;
        prcp_out = max(0, prcp_out);
        prcp_out(nanmap) = nan;          % reset NA's to nan
        
    elseif (strcmp(direction, "reverse"))
        
        zero_map = prcp_in <=0;
        
        spr_min = pr_min ^ scaling;
        
            % shift the scaled data back to the right place, the undo the scaling.
        sprcp = prcp_in + spr_min - dx/2;
        
        prcp_out = sprcp .^ (1/scaling);
        prcp_out(zero_map) = 0;
        
        sigma = [];
        kurt  = [];
        bins  = [];
        edges = [];
            
    else
        error("error:  bad direction:  %s\n", direction);        
    end
        
end