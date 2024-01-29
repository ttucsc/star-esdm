function [FILT, filt] = calc_filter(nterms, sig_term, myrs, yrlen)
%FILT = calc_filter(nterms, sig_term, myrs)
%
%   Returns fourier domain filter for an ideal rectangular filter with nterms
%   convolved with gaussian of sigterm, for myrs*yrlen days.
%   nterms is the number of cycles per yrlen of the starting ideal filter.
%   Filter is based on a 365-day year, so an nterms of 5 with sig_term of 3 means start from
%   an ideal filter passing frequencies of up to 5 cycles per 365-day year,
%   then convolve it with a guassian with a sigma of 3.
%
%   Inputs:
%       nterms          # of terms per year to keep 
%                           (positive terms only, not including DC or negative terms.)
%                           nterms=5 starts with a rect. filter of length 11
%                           nterms is relative to 1 year.  Will be scaled up by mterms.
%       sig_term        sigma for gaussian to convolve with rect. filter
%                           s/b specified relative to 1 year.  Will be scaled up by mterms
%       myrs            # of years to create filter for. 
%                           (1, 41, 141, etc., for 1 year, or 41 yrs of obs data, or 141 yrs of model data...)
%
%   Returns:
%       FILT            Fourier domain filter of length yrlen*myrs.
%       filt            circular spatial domain filter of length yrlen*myrs.
%                           note:  filt is centered at 1, intended for use
%                           with cconv(...), circular convolution.
%                           Centering at 1 means filter does not shift
%                           signal.
%                           filt is not calculated if only 1 return
%                           variable is detected.

    if (nargin < 4), yrlen=365; end

        % the rectangular (ideal) low-pass filter
    len = yrlen*myrs;
    nterms = floor(nterms*myrs);
    sig = sig_term*myrs; 
    Frect = zeros(1,len);

    Frect([1:(nterms+1), (end-nterms+1):end]) = 1;
    
                % the gaussian
    midpt = ceil((len+1)/2);        % 11: -> 6.  12:  ->7.
    GF=gauss(len,sig, midpt);
    GF=circshift(GF,1-midpt,2);     % shift peak back to GF(1).
                % now convolve Frect with GF to make filter
    FILT = cconv(Frect,GF,len);
    FILT = FILT/FILT(1);          % make sure that DC is passed at full power.
    
    if (nargout > 1)
        filt = real(ifft(FILT));    % *should* be purely real, but may have tiny fractional imag. part due to computational limitations...~ 10^-15.
    end
end