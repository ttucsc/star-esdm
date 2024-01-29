function pdfs_adj = adjust_kde_pdfs(pdfs, sigs_orig, bins)
% function pdfs_adj = adjust_kde_pdfs(pdfs, slvrsigs, bins)
% Adjusts the pdfs to remove the spreading effect of convolving the initial histograms with a gaussian in the KDE
% process.

    [nbins, yrlen] = size(pdfs);
    x = bins;
    if (isrow(x)), x=x'; end

%     mus1 = nansum(x .* pdfs);   % mean daily value, s/b around midpoint of bins.
% 
%     sigs1 = sqrt(nansum(x.^2 .* pdfs) - mus1.^2); % 
    
    [sigs, mus] = pdf_sigmas(pdfs, bins, [], false);        % sigs is widened from original sigmas by KDE smoothing.
    
%     sigs_adj_old = sqrt(sigs1.^2 + slvrsigs.^2);
    sigs_adj = sigs;
                % now adjust for spreading from convolving w/ silversigs.
    pdfs_adj = zeros(nbins, yrlen);

    for i=1:yrlen
            % calculate an adjusted x axis representing the 'true' x axis for the 
            % kde-widened distribution.  We then interpolate this curve to find
            % interpolate this curve at the standard bin locations.
%         x_adj_old  = mus1(i) + (x-mus1(i)) * sigs(i)/sigs_adj_old(i);
        x_adj  = mus(i) + (x-mus(i)) * sigs_orig(i)/sigs_adj(i);
%               means   +  slightly narrowed x values... gives us a corrected x-axis for the KDE-widened pdf on ith day.
        flags = ~isnan(pdfs(:,i)); 
        pdf = pdfs(:,i);            
%       try
            pdfs_adj(:,i) = interp1(x_adj(flags), pdf(flags), x,'pchip',nan);
%         catch me
%             fprintf("oops! adj_kde_pdfs:  problem adjusting for convolving gaussian with gaussian\n");
%             msgtext=getReport(me);
%             fprintf('%s\n', msgtext);
%             fprintf('------\n');        
%         end
    end  
         % normalize adjusted pdfs to make them pdfs.   
    negs = pdfs_adj(:) < 0;
    if (any(negs))
%       if (min(pdfs_adj(:)) < -1e-18)
        if (min(pdfs_adj(:)) < -1e-14)
            fprintf(2, 'warning:  adjust_kde_pdfs: %d negatives (most negative:  %.6g\n', sum(negs), min(pdfs_adj(negs)));
        else
            fprintf(      'note:  adjust_kde_pdfs: %d negatives (most negative:  %.6g\n', sum(negs), min(pdfs_adj(negs)));            
        end
        pdfs_adj(negs)=0;
    end
    ss = nansum(pdfs_adj);
    pdfs_adj = pdfs_adj ./ ss;
end
