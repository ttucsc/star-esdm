function [ shift ] = ic_corr_1d( base, changed, nterms, nterp, no_DC, fignum )
%   Does fft-based correlation to get subpixel registration on 1-D data
%
%   Correlation done by standard fourier-domain correlation, equivalent of circular correlation is spatial domain.
%   This version works well if signal is cyclical and changed signal is reasonably similar to base signal.
%   Does not produce a spike at correlation point...instead provides a smooth continuum, and simply returns the point
%   with the highest correlation value.
%
%   does:  ifft( fft(changed) .* conj(fft(base)) ), then finds max.
%
%               using only the first nterms terms of the fft.
%
%   If nterp > 1, expands signal by factor of nterp, so effectively interpolates to sub-cell precision.
%
%   Assumes changed and base are same size, and of similar signal levels.
%
%   To interpolate to sub-step ('sub-pixel') accuracy, set nterp to expand signal as needed.  (nterp = 10 will return result to
%   nearest 1/10 cell.
%   
%       (As an alternative, see e.g. phase_corr_1d(...), based on Foroosh, "Extension of Phase Correlation to Subpixel
%       Registration", IEEE 2002.)
%
% intputs:
%   base        original signal
%   changed     shifted signal.  response is amount that this has shifted relative to base
%
%   nterms      # of terms to keep (1 = DC, 2 = DC + 1 cycle, etc.);  0 or empty:  keep all terms
%   nterp       # of times to expand signal to find max.  
%   no_DC       exclude DC term (sets average of each signal to zero)
%   fignum      figure to plot in.  if absent, no plot is generated.
%
% 
    if (~exist('nterms','var') || isempty(nterms) || nterms < 1), nterms= 0;     end
    if (~exist('nterp', 'var') || isempty(nterp)  || nterp  < 1), nterp = 1;     end
    if (~exist('no_DC', 'var')),                                  no_DC = false; end
    
    if (any(isnan(base)) || any(isnan(changed)))
        shift = nan;
        return;
    end
    
    if (~isrow(base))
        base = base';
    end
    if (~isrow(changed))
        changed = changed';
    end
    
    len = length(base);
    
    B = fft(base);
    C = fft(changed);
    if (no_DC)
        B(1)=0;
        C(1)=0;
    end
    
    BC = C .* conj(B);
    
    if ( nterms > 0)
        bc = ifft_resample(BC, nterp*len, nterp*(nterms-1));
    else
        bc = ifft_resample(BC, nterp*len);
    end
    
    ix = find(bc == max(bc),1);
    
    if (exist('fignum','var'))
        bf = fft_filt(base,nterms-1);       % fft_filt's nterms doesn't count DC term.
        cf = fft_filt(changed, nterms-1);
        figure(fignum);
        subplot(2,1,1);
        plot(1:length(base), base, 'b-', 1:length(bf), bf, 'b:', 1:length(changed), changed, 'r-', 1:length(cf), cf, "r:");
        subplot(2,1,2);
        plot(1:length(bc), bc, [ix,ix],[min(bc),max(bc)]);
    end

    shift = (ix-1)/nterp;
    
    if (shift > len/2)
        shift = shift - len;
    end
end
