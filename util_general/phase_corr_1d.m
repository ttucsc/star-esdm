function [ shift ] = phase_corr_1d( base, changed, gravwidth, nterms, maxshift, fignum )
%   Does Phase-based correlation to get subpixel registration on 1-D data
%
%       Based on Foroosh, "Extension of Phase Correlation to Subpixel Registration", IEEE 2002.
%
%   Correlates via cross-power spectrum:
%       exp(-i(ux)) = Y2(u) .* conj*Y1(u) ./ magnitude(Y1(u) .* conj(Y1(u))
%
%   THIS IS A WORK IN PROGRESS.  USE gravwidth=2, ian!
%
%   Works well if there are sharp edges in the signal, or if changed is very similar to base, except shifted.  
%   Is rather sensitive to noise otherwise.
%   Basically divides out the magnitude of each fourier term, so that higher frequency terms have as
%   much influence on the result as lower frequency terms, even though the low-frequency terms have much 
%   more power in the straight fourier transform.
%   
% intputs:
%   base        original signal
%   changed     shifted signal.  response is amount that this has shifted relative to base
%   gravwidth   gravity width:  # of elements over which to calculate the ctr of gravity
%                           probably should be odd number...
%                   1   gives nearest cell (distance of peak from ctr)
%                   2   returns ctr of gravity of peak and next highest pt
%                   3   ctr of gravity & 1 either sie of peak
%
%   nterms      # of terms to keep (1 = DC, 2 = DC + 1 cycle, etc.)
%
% 
    if (~exist('nterms','var') || isempty(nterms) || nterms==0); nterms=[]; end;
    if (~isrow(base))
        base = base';
    end
    if (~isrow(changed))
        changed = changed';
    end
    
    len = length(base);
    
    B = fft(base);
    C = fft(changed);
    
    BC = (C .* conj(B)) ./ abs(B .* conj(B));
    
    if (~isempty(nterms))
        keepers = BC(1:nterms);
        BC = [keepers, fliplr(conj(keepers(2:end)))];
        bc = ifft(BC);
    else
        bc = ifft(BC);
    end
    bclen = length(bc);
            % mask off any points outside acceptable range
    if (exist('maxshift','var') && ~isempty(maxshift))
        maxshift = ceil(maxshift*bclen/len);
        range = mod((-maxshift:maxshift),bclen)+1;
        flt = zeros(1,bclen);
        flt(range)=1;
        bc = bc .* flt;
    end
    ix = find(bc == max(bc),1);
    
    if (exist('fignum','var'))
        figure(fignum);
        subplot(2,1,1);
        plot(1:length(base), base, 1:length(changed), changed);
        subplot(2,1,2);
        plot(1:length(bc), bc, [ix,ix],[min(bc),max(bc)]);
    end

    shift = ix-1;

    if (gravwidth > 1)
        bcabs = abs(bc);
        
        wid1 = floor(gravwidth/2);
        ix1 = ix-wid1;
        ix2 = ix+wid1;
        wts = -wid1:wid1;
        rng = mod((ix1:ix2)-1,bclen)+1;
%        fprintf('range:  ');  fprintf('%d ',rng); fprintf('\n');
        if (mod(gravwidth,2)==0)
            prv = rng(wid1);
            nxt = rng(wid1+2);
            if (bc(prv) > bc(nxt))
                wts=wts(1:end-1);
                rng=rng(1:end-1);
            else
                wts=wts(2:end);
                rng=rng(2:end);
            end
        end
                % calculate ctr of gravity for pts of interest
        frac = sum( wts .* bcabs(rng)) / sum(bcabs(rng));
        
        shift = shift + frac;
    end
    
    if (shift > length(bc)/2)
        shift = shift - length(bc);
    end
    shift = shift  * len / bclen;
end
