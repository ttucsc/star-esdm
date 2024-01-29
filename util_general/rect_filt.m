function dout = rect_filt(d, filtlen, dim, doGauss)
% circular low-pass filter
%   inputs:
%       d       data to filter  (can be 1D, 2D or 3D.  
%       filtlen length of rectangular filter, or rect. equiv. length for gaussian filter
%       dim     (optional) dimension to filter along.  if empty or missing, filter along 1st non-singleton dimension
% %       doGauss (optional) boolean.  If true, filter with gaussian pseudo-equivalent of rect. filter.

%       MATH Notes for gaussian:
% notes re sigmas, gaussians & fft's of gaussians:
%   1.  sigma of uniform distribution of length n: 
%           for continuous:  (n/2)/sqrt(3) = n/sqrt(12)
%           for discrete:  sqrt((n^2-1)/12) , 
%           which, as m -> large, -> n/sqrt(12).
%   1.a. n, length of uniform dist for given sigma:  
%           n = sqrt(12*(sigma^2) + 1) 
%   2.  sigma-0 of gaussian s.t. sigma(fft(gaussian)), in freq = sigma in time (for n points): sqrt(n/(2*pi))
%   3.  sigma (freq. domain) for a sigma=1 (time domain) is .1592 = 1/(2*pi)
%   4.  sigma-fft for arbitrary time domain sigma:  sigma-fft = sigma0^2 / sigma(time)
%   4.a.    and same for given sigma-fft:  sigma(time domain) = sigma0^2 / sigma(freq domain)
%   

    [ny,nx,nz] = size(d);
    sz = size(d);
    len = sz(dim);
    
    if (~exist('dim','var') || isempty(dim))
        if (ny>1),          dim=1;
        elseif (nx > 1),    dim=2; 
        else,                dim=3; 
        end
    end
    
    if (~exist('doGauss','var') || isempty(doGauss))
        doGauss = false;
    end

    if (~doGauss)
            % rectangular filter
        filt=zeros(len,1, class(d));
        filtlen=round(filtlen);
        if (mod(filtlen,2)==0), filtlen=filtlen+1; end
        filtlen=min(filtlen, ny);
        filt(1:filtlen)=1;
        filt=circshift(filt,-floor(filtlen/2),1);
        filt=filt/sum(filt);
        FILT = fft(filt);
    else
        % gaussian filter, with equivalent sigma of rectangular filter with length filtlen 
        filtlen=min(filtlen,len);
        sig=sqrt((filtlen^2-1)/12);         
        FILT=GAUSS_fft(len, sig, true)';
    end
    
    if (dim==1)
        dout = squeeze(real(ifft(fft(d,ny,1) .* repmat(FILT,1,nx,nz),ny,1)));
    elseif (dim==2)
        dout = squeeze(real(ifft(fft(d,nx,2) .* repmat(FILT',ny,1,nz),nx,2)));
    else
        FILT = reshape(FILT,1,1,[]);
        dout = squeeze(real(ifft(fft(d,nz,3) .* repmat(FILT,ny,nx,1),nz,3)));
    end
end


% old version:  worked for 2d only
% % circular low-pass filter
% %   inputs:
% %       d       data to filter
% %       filtlen length of rectangular filter, or rect. equiv. length for gaussian filter
% %       dim     (optional) dimension to filter along.  if empty or missing, filter along 1st non-singleton dimension
% % %       doGauss (optional) boolean.  If true, filter with gaussian pseudo-equivalent of rect. filter.
% 
% %       MATH Notes for gaussian:
% % notes re sigmas, gaussians & fft's of gaussians:
% %   1.  sigma of uniform distribution of length n: 
% %           for continuous:  n/2/sqrt(3) = n/sqrt(12)
% %           for discrete:  sqrt((n^2-1)/12) , 
% %           which, as m -> large, -> n/sqrt(12).
% %   1.a. n, length of uniform dist for given sigma:  
% %           n = sqrt(12*(sigma^2 + 1)) 
% %   2.  sigma-0 of gaussian s.t. sigma(fft(gaussian)), in freq = sigma in time (for n points): sqrt(n/(2*pi))
% %   3.  sigma (freq. domain) for a sigma=1 (time domain) is .1592 = 1/(2*pi)
% %   4.  sigma-fft for arbitrary time domain sigma:  sigma-fft = sigma0^2 / sigma(time)
% %   4.a.    and same for given sigma-fft:  sigma(time domain) = sigma0^2 / sigma(freq domain)
% %   
% 
%     [ny,nx] = size(d);
%     if (~exist('dim','var') || isempty(dim))
%         if (ny>1),  dim=1;
%         else,       dim=2; end
%     end
%     
%     if (~exist('doGauss','var') || isempty(doGauss))
%         doGauss = false;
%     end
%     
%     if (~doGauss)
%         if (dim==1)
%             filt=zeros(ny,1);
%             filtlen=round(filtlen);
%             if (mod(filtlen,2)==0), filtlen=filtlen+1; end
%             filtlen=min(filtlen, ny);
%             filt(1:filtlen)=1;
%             filt=circshift(filt,-floor(filtlen/2),1);
%             filt=filt/sum(filt);
%             FILT=repmat(fft(filt),1,nx);
%             D=fft(d,ny,1);
%             dout=real(ifft(D .* FILT));
%         else
%             filt=zeros(1,nx);
%             filtlen=round(filtlen);
%             if (mod(filtlen,2)==0), filtlen=filtlen+1; end
%             filtlen=min(filtlen, nx);
%             filt(1:filtlen)=1;
%             filt=circshift(filt,-floor(filtlen/2),2);
%             filt=filt/sum(filt);
%             FILT=repmat(fft(filt),ny,1);
%             D=fft(d,nx,2);
%             dout=real(ifft(D .* FILT));
%         end
%     else        
%         if (dim==1)
%             filtlen=min(filtlen,ny);
%             sig=sqrt((filtlen^2-1)/12);
%          
%             G=GAUSS_fft(ny, sig, true)';
%             G=repmat(G,1,nx);
%             D=fft(d,ny,1);
%             dout=real(ifft(D .* G));
%         else
%             filtlen=min(filtlen,nx);
%             sig=sqrt((filtlen^2-1)/12);
%             G=GAUSS_fft(nx, sig, true);
%             G=repmat(G,ny,1);
%             D=fft(d,nx,2);
%             dout=real(ifft(D .* G, nx, 2));
%         end        
%     end
% end
% 
