function d_out = lpf_2d(d, nx, ny, keep_mean, sigx, sigy)
%   Does 2-D low-pass filter of FFT of input rectangle.
%       need to add rounding off in x & y by sigx and sigy.
%

    if (~exist('keep_mean','var') || isempty(keep_mean)), keep_mean = true; end

    D=fft2(d);
    [sy,sx] = size(D);
    
    if (exist('sigx','var') && ~isempty(sigx))
        DF = make_filt(nx, ny, sigx, sigy, sx, sy);
    else
        DF=zeros(size(D));
        DF(      1:ny+1,       1:nx+1)=1;
        DF(sy-ny+1:sy,         1:nx+1)=1;
        DF(      1:ny+1, sx-nx+1:sx)=1;
        DF(sy-ny+1:sy,   sx-nx+1:sx)=1;
    end
    
    if (~keep_mean), DF(1,1)=0; end
    D = D .* DF;
    d_out = real(ifft2(D));
    
%     figure(99);  
%     subplot(3,1,1);
%     surf(d,'edgecolor','none');
%     subplot(3,1,2);  
%     surf(DF,'edgecolor','none');
%     subplot(3,1,3);
%     surf(d_out,'edgecolor','none');
%     lights_on;
end

function F = make_filt(nx, ny, sigx, sigy, sx, sy)

    F = zeros(sy,sx);

    xFILT = calc_filter(nx, sigx, 1, sx);
    yFILT = calc_filter(ny, sigy, 1, sy);
    
    for i=1:sy
        F(i,:) = xFILT * yFILT(i);
    end
end
    

