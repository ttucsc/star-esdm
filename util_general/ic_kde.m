function [pdf_out] = ic_kde(my_pdf, npts, do_correction)
%does KDE on my_pdf (which can be either a pdf or a histogram)
%
%   NOTE:  my_pdf must be sampled at even intervals. 
%   if do_correction is true, then it corrects for the spreading that occurs from convolving two gaussians.
%   Set to true IF:
%       1.  data is relatively close to gaussian
%       2.  silversig is more than 1/10 std dev of your data.
%               gaussian spreading:  sqrt(sig1^2 + sig2^2)/sig1^2  

    if (~exist('do_correction','var') || isempty(do_correction))
        do_correction = false;
    end
    if (isempty(npts))
        npts = sum(~isnan(my_pdf));
        if (npts < 2), error("error:  my_pdf (sum=%.2f) must be a histogram if npts is empty", npts); end
    end
    if (iscolumn(my_pdf))
        my_pdf = my_pdf';
    end

    nbins = length(my_pdf);
    bins = 1:nbins;
    s = sum(my_pdf);
    [sig,mu] = pdf_sigmas(my_pdf, bins, [],true);
    silversig =  1.06 * sig .* (npts.^(-.2));
    GS = GAUSS_fft(nbins, silversig, true);
    if (1-sum(GS) > 1e-8), error("error:  silvermans' sigma is too narrow to convolve properly with the data"); end
    pdf_mid = abs(ifft(GS .* fft(my_pdf)));    
    
    if (do_correction)
        try
            bins_adj  = mu + (bins-mu) * sig/sqrt(sig^2+silversig^2);
            pdf_out = interp1(bins_adj, pdf_mid, bins,'makima',0);
            pdf_out = pdf_out * s / sum(pdf_out);   % normalize
        catch
            oops();
        end
    else
        pdf_out = pdf_mid * s / sum(pdf_mid);   % normalize
    end
        
end

