function [yout, YOUT] = lpf_FILT(y, FILT)

%   returns y circularly filtered with fourier-domain filter FILT
%   also returns Y_OUT, the fourier transform of y_out.
%

        % make sure y & filt are same shape...
    if ((iscolumn(y) && isrow(FILT)) || (isrow(y) && iscolumn(FILT)))
        FILT = FILT';
    end
    
    Y=fft(y);
    YOUT = Y .* FILT;
    yout = real(ifft(YOUT));
end