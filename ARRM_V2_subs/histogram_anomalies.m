function [hist_tbl, lohi ] = histogram_anomalies( anoms, na_map, yr_range, rolling_steps, pdf_yrs, edges, do_fourier, do_gauss, yrlen, rolling_range )
%   returns one or more histograms of model anomalies
%
%   If yrstep == 0 or yrstep >= # of years, then a single (yrlen X nedges) histogram is created and returned.
%   Otherwise, yr_step and pdf_yrs define the number of histograms to create and the number of years' data to include in
%   each histogram.
%       nyrs  (# of years' data) = length(anoms) / yrlen;
%       nsets (# of histograms)  = ceil((nyrs-pdf_yrs) / yr_step).  
%           first histogram uses year 1 to year pdf_yrs
%           last  histogram uses year (nyrs-pdf_yrs+1) to nyrs
%
%   Inputs:
%       anoms           anomalies to be histogrammed
%       na_map          map of location of NAs in the data.  These data points will be excluded in calculating the
%                           histograms.  Care should be taken that the presence of NAs will not bias the results.  
%                           NOTE:  No correction for NA's presence is done here.
%       yr_range        range of years to histogram over (relative to start of data.)  1-based.  1 = start at beginning (1st year) of anoms.
%       rolling_steps   years (relative to yr_range(1)) where rolling histograms should be calculated.  usually
%                           calculated via a call to Dataparams.calc_rolling_steps(RunParams.pdf_yrstep).
%                           if empty, a single histogram is calculated, using all data in yr_range.
%       pdf_yrs         # of years to include in each histogram
%       edges           bin edges to use for histogramming
%       do_fourier      boolean.  if true, uses (circular) fourier transform, after extended the data to avoid circular
%                           overlap problems.
%       do_gauss        boolean.  if true, uses gaussian kernel with equivalent smoothing as a rectangular filter of
%                           width pdf_yrs.  
%       yrlen           # of values per year.  Usually 365.
%       rolling_range   start, end indexes of data to use for each rolling step.  (needed only when doing rolling histograms
%
%   Outputs
%       hist_tbl        histograms.  1 or nsets.  If 1 histogram,         size is (nedges-1) X yrlen
%                                                 If multiple histograms, size is (nedges-1) X yrlen X nsets
%                       units are counts of occurrences.  Not normalized to 1.
%       lohi            2-element count of the number of data points outside the bounds of edges.
%
%   ToDo:
%
%       add a histogram-warping function to do histogram equalization or some other remapping of day-or-year
%       to adjust for climatology.
%
%________________________
%
%       MATH Notes:
%       The following gives an explanation of how to get a sigma for a gaussian with the equivalent standard deviation
%       of a rectangular filter of a given length.  There are 2 sigmas of interest:  the sigma for a time-domain
%       gaussian to convolve with the time-domain signal, or the sigma in the frequency domain, to multiply
%       (point-by-point) with the fourier transform of the time-domain signal.  There is a useful value, sigma-0, which
%       is the sigma at which the time domain and frequency domain sigmas are the same.
% notes re sigmas, gaussians & fft's of gaussians:
%   1.  std dev (sigma) of uniform distribution of length n: 
%           for continuous:  n/2/sqrt(3) = n/sqrt(12)
%           for discrete:  sqrt(((n^2-1))/12) , 
%           which, as m -> large, -> n/sqrt(12).
%   1.a. n, corresponding length of uniform (rectangular) distribution for given sigma  (solving discrete for n):  
%           n = sqrt(12*(sigma^2) + 1) 
%   2.  sigma-0 of gaussian s.t. sigma(fft(gaussian)), in freq = sigma in time (for n points): sqrt(n/(2*pi))
%   3.  sigma (freq. domain) for a sigma=1 (time domain) is .1592 = 1/(2*pi)
%   4.  sigma-fft for arbitrary time domain sigma:  sigma-fft = sigma0^2 / sigma(time)
%   4.a.    and same for given sigma-fft:  sigma(time domain) = sigma0^2 / sigma(freq domain)
%
%_______________________

    if (~exist('do_fourier','var') || isempty(do_fourier)), do_fourier = false; end
    if (~exist('do_gauss',  'var') || isempty(do_gauss  )), do_gauss   = false; end
    if (~exist('yrlen','var'     ) || isempty(yrlen)     ), yrlen      = 365;   end

            % get the locations of NAs, if not provided
    if (isempty(na_map))
        na_map = isnan(anoms);
        nnas = sum(n_map);
    else
        nnas = sum(na_map);
        if (nnas > 0), anoms(na_map) = nan; end
    end
    
        % check for non-linear binning.
    ddx = diff(edges);  
    dx = edges(2) - edges(1);
    if (range(ddx) > 1e-6*dx)
        nonlinear_binning = true;
    else
        nonlinear_binning = false;
    end
    
    npts = length(anoms);
    totyrs = npts/yrlen;
    nyrs = yr_range(2)-yr_range(1)+1;
    do_single=isempty(rolling_steps);
    if (mod(nyrs, 1) ~= 0),  error('histogram_anomalies:  error:  invalid yr_range'); end
    if (mod(totyrs,1) ~= 0), error('histogram_anomalies:  error:  data is not even multiple of yrlen'); end
    if (do_single && do_gauss)
        rolling_steps = nyrs/2;
        if (isempty(pdf_yrs) || pdf_yrs==0 || pdf_yrs >= nyrs), pdf_yrs = nyrs; end
    end
    
    lohi = [sum(anoms(:)<edges(1)),sum(anoms(:)>edges(end))];       % # of data points out of range.
    nbins = length(edges)-1;                                        % # of histogram bins
    
    if (do_single && (nonlinear_binning || ~do_gauss))      % just one histogram requested, and gaussian smoothing not requested.  histogram each day and we're done.
        hist_tbl = zeros(nbins,yrlen);
        anoms = reshape(anoms, yrlen,totyrs);
        for i=1:yrlen
            try
                hist_tbl(:,i) = histcounts(anoms(i,yr_range(1):yr_range(2)), edges);
            catch
                oops();
            end
                
        end
    elseif (nonlinear_binning)
%            [start_ix, end_ix, nyrs] = obj.DP.using_range(yr_type, rolling_set, obj.RP.pdf_yrstep);
        nrolling_steps = length(rolling_steps);
        hist_tbl = zeros(nbins, yrlen, nrolling_steps);
        for j=1:length(rolling_steps)
            bstart = rolling_range(j,1);
            bend   = rolling_range(j,2);
            data = reshape(anoms(bstart:bend),yrlen, []);
            for i=1:yrlen
                hist_tbl(:,i,j) = histcounts(data(i,:), edges);
            end
        end
                
            
    else

        pdf_yrs = floor(pdf_yrs);   % make sure it's an integer!
        if (mod(pdf_yrs,2) ~= 1), pdf_yrs = pdf_yrs + 1; end   % make sure it's odd.

%       nbins = length(edges)-1;                       % # of histogram bins

        fudge=1e-15*(edges(end)-edges(1));  % computer-arithmetic truncation error for double-precision divides to get bin indexes.
                                            % Avoids error when value is exactly on bin edge.

        ix=(0:(npts-1))';
        iyrs = floor(ix/yrlen) + 1;
        doy = mod(ix,yrlen)+1;
        anom_ix = max(0, min(nbins-1,floor((anoms-edges(1))/dx + fudge)))+1;      % integer index for each anomaly.
                                                                                  % NOTE:  NAs will be set to max bin.
            % We'll make the table as small as possible so histogramming is
            % faster.  Then append zeros at the end.
        min_ix = min(anom_ix);
        max_ix = max(anom_ix);
        startz = min_ix-1;                  % # of zeros to add at low end of hist_tbl when we're done
        endz   = nbins  - max_ix;           % # of zeros to add at top end of hist_tbl when we're done
        nix    = max_ix - startz;           % # of anom_ix's

        anom_ix = anom_ix - startz;

        try
            lcl_tbl = zeros(nix,yrlen,totyrs);        % nans will go in the extra row, and we'll zero them later.
        catch
            oops();
        end
%         for i=1:npts
%             lcl_tbl(anom_ix(i), doy(i), iyrs(i)) = 1;
%         end
%         fprintf('in histogram anomalies, min,max are: %d %d   min,max anomaly are: %.3f %.3f\n', min_ix, max_ix, min(anoms(:)), max(anoms(:))); 
        
        if (nnas > 0)               % remove NA [points if we have any.
            anom_ix(na_map) = [];
            doy(na_map) = [];
            iyrs(na_map) = [];
        end
            
        lcl_tbl(sub2ind([nix, yrlen, totyrs], anom_ix, doy, iyrs))=1;           % this should be faster than iterating through the array.
        
        
        half_pdf_yrs = floor(pdf_yrs/2);

        if (do_fourier)
                % fft is fastest w/ exact power of 2, so 
                % find 1st power of 2 >= # of years we need.
            minyrs = totyrs + 2*pdf_yrs;
            myrs=2;
            while (myrs<=minyrs), myrs = myrs*2; end
            mstart = floor((myrs-totyrs)/2);                  % # of years to pad with at start
            mend = myrs-totyrs-mstart;                        % # of years to pad with at end

                % pad either end with mean of 1st or last 1/2*pdf_yrs years.
            mean_start = mean(lcl_tbl(:,:,1:half_pdf_yrs),3);
            mean_start = repmat(mean_start, 1,1,mstart);
            mean_end   = mean(lcl_tbl(:,:,(end-half_pdf_yrs+1):(end)),3);
            mean_end   = repmat(mean_end,1,1,mend);

            lcl_tbl = cat(3, mean_start, lcl_tbl, mean_end);

                % smooth via FFT
            nz = size(lcl_tbl, 3);
            if (~do_gauss)            
                KRNL = fft(circshift([ones(pdf_yrs,1);zeros(nz-pdf_yrs,1)],-half_pdf_yrs,1));
            else
                    % get gaussian w/ equivalent smoothing as pdf_yrs-year filter
                    % see math notes in climatology_3d_gauss_V2.m
                sigma = sqrt((pdf_yrs^2-1)/12);
                KRNL = GAUSS_fft(nz,sigma, true);
            end
            KRNL = reshape(KRNL,1,1,nz);

                    % note:  for versions of Matlab earlier than R2016b, 
                    % would have to repmat the KRNL to the full size for the
                    % .* :
    %         KRNL = repmat(KRNL,nc,nr,1);

                    % do convolutions in fourier domain.
            tbl = real(ifft(fft(lcl_tbl,nz,3) .* KRNL,nz,3));   % NOTE: tbl has one plane for each day...
                                                                % so effectively each day is a pdf, normalized to 1.

            % now keep just 1 set every yr_step years.

            keepers = false(nz,1);
            keepers(mstart + yr_range(1) + rolling_steps) = true;
            hist_tbl = tbl(:,:,keepers);
            outlen = size(hist_tbl,3);
            
            hist_tbl = hist_tbl * pdf_yrs;      % scale hist_tbl to actual number of events in each histogram.
            hist_tbl = max(0,hist_tbl);
        else
                % fft is fastest w/ exact power of 2, so 
                % find 1st power of 2 >= # of years we need.
            mstart = half_pdf_yrs;                          % # of years to pad with at start
            mend   = pdf_yrs - mstart;                      % # of years to pad with at end

                % pad either end with mean of 1st or last 1/2*pdf_yrs years.
            mean_start = mean(lcl_tbl(:,:,1:half_pdf_yrs),3);
            mean_start = repmat(mean_start, 1,1,mstart);
            mean_end   = mean(lcl_tbl(:,:,(end-half_pdf_yrs+1):end),3);
            mean_end   = repmat(mean_end,1,1,mend);

            lcl_tbl = cat(3, mean_start, lcl_tbl, mean_end);        % pad the local histogram table.
            tbl_len = size(lcl_tbl, 3);

                % smooth by averaging w/ pdf_yrs_year rectangular filter

            outlen = length(rolling_steps);       
            hist_tbl = zeros(nix, yrlen, outlen);            
            for i=1:outlen
                istart = mstart + yr_range(1) + rolling_steps(i) - half_pdf_yrs;     % each set starts            
                iend   = min(tbl_len, istart+pdf_yrs-1);
                istart = min(istart, iend-pdf_yrs+1);
                hist_tbl(:,:,i) = sum(lcl_tbl(:,:,istart:iend),3);

            end
        end

    %         % add zeros to fill out the table to size nyvals;
        if (startz > 0)
            zs = zeros(startz,yrlen,outlen);
            hist_tbl = cat(1, zs,hist_tbl);
        end
        if (endz > 0)
            zs = zeros(endz,yrlen,outlen);
            hist_tbl = cat(1, hist_tbl, zs);
        end
        

    end    
    
%     bins = edges(1:end-1);
%     if (ndims(hist_tbl)~=3)
%         figure(17);
%         surf(1:365,bins,hist_tbl,'edgecolor','none');
%         view([-89,1]);
%     else
%         for i=1:size(hist_tbl,3)
%             figure(i);
%             surf(1:365,bins,hist_tbl(:,:,i),'edgecolor','none');
%             view([-89,1]);
%         end
%         figure(16);
%         tothist = sum(hist_tbl,3);
%         surf(1:365,bins,tothist,'edgecolor','none');
%         view([-89,1]);
%         for i=16:-1:1; figure(i); end; figure(16); 
%         
%     end
% 
%     fprintf('\n');
%     elapsed=toc(t1);
%     fprintf('histogramming:                %8.4f\n', elapsed);  
end

