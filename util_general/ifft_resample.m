function [yout] = ifft_resample(Y, outlen, terms, dim)
% [yout] = ifft_resample(Y, outlen, terms, dim)

%   Resamples y to length outlen using fft, and optionally filter to keep only some terms of FFT.
%   operates along 1st non-singleton dimension of y.
%   returns:  
%       yout        resampled data
%
%   inputs:
%       Y           FFT of data to be resampled or filtered
%       outlen      length of desired output
%       terms       either an array of 1's identifying frequency terms to keep, or
%                   a scalar representing the highest frequency term to keep.  
%                   -1 :  keep all terms
%                       note:  term 0 is DC, 1 is 1 cycle, 2 is 2 cycles, etc.
%                       note:  because of ambiguity, terms == 1 means 'DC + first term'.
%                       note:  if terms is non-scalar, it can either be a vector (applied to all fft sets)
%                              or matrix with flags unique to each set.
%                       note:  terms only needs to be as long as needed to specify all terms desired.  Will be filled
%                                   with zeros as needed.
%

    if (length(size(Y))>2)
        throw(MException('IFFT_RESAMPLE','input dimension too high.  Can only operate on vectors or 2-D matrices'));
    end
      
    if (~exist('terms','var'))
        terms=[];
    end
    if (~exist('dim','var') || isempty(dim))
        dim = find(size(Y)>1,1);
    end
    n=size(Y,dim);

    if (isempty(terms))
        terms=-1;
    end
    if (isscalar(terms))
        if (terms == -1 && outlen==n)      % nothing to do.  copy input & return.
            yout=ifft(Y,n,dim);
            return;
        end

        if (terms < 0)      % keep all terms
            nn = min(n, outlen);
            terms=ones(1,ceil(nn/2));
        else        % if scalar, then keep all terms up to single value given in terms
            terms=ones(1,terms+1);
        end
    end

    [nrows,ncols] = size(Y);
    
    if (outlen == n)
        Yout = Y;
    else
        nzeros = max(0,outlen-n);
        nkeep = min(outlen, n);             % # frequencies to copy.
        m1 = ceil(nkeep/2);                 % pos. freq's to keep
        m2 = n - floor(nkeep/2)+1;      % start of neg freq's to keep
        if (outlen > n)                 
            if (dim == 1)
                myzeros=zeros(nzeros,ncols);    % array of zeros to insert if needed.
                Yout = [Y(1:m1,:);myzeros;Y(m2:end,:)];
            else
                myzeros=zeros(nrows, nzeros);    % array of zeros to insert if needed.
                Yout = [Y(:,1:m1),myzeros,Y(:,m2:end)];
            end
        else
            if (dim == 1)
                Yout = [Y(1:m1,:);Y(m2:end,:)];
            else
                Yout = [Y(:,1:m1),Y(:,m2:end)];
            end
        end
    end

    
            % if terms-to-keep were supplied, turn them into a mask to zero out terms not wanted.
    if (~isempty(terms))  
        if (sum(terms ~= 0) ~= sum(terms == 1))
            terms(terms ~= 0) = 1;      % make sure terms are all zeros or ones.
        end
        if (isvector(terms))            % if given a vector, replicate it into a matrix.
            if (dim==1)
                if (isrow(terms))       % could be wrong orientation.  fix it!
                    terms=terms.';
                end
                terms=repmat(terms, 1,ncols);
            else
                if (iscolumn(terms))       % could be wrong orientation.  fix it!
                    terms=terms.';
                end
                terms=repmat(terms, nrows,1);
            end
        end
        
        nt = size(terms, dim);                  % terms is now matrix of positive-frequency terms.
        nzeros = max(0, outlen - (2*nt-1));     % now clone it to the negative frequencies.
        if (dim==1)
            terms=[terms;zeros(nzeros,ncols);flipud(terms(2:end,:))];
        else
            terms=[terms,zeros(nrows,nzeros),fliplr(terms(:,2:end))];
        end
        Yout = Yout .* terms;
    end
        
    yout = ifft(Yout);
    
end
    
% function [yout] = fft_resample(y, outlen, dim)
% %[yout] = fft_resample(y, outlen)
% %   Resamples y to length outlen using fft.
% %   operates along 1st non-singleton dimension of y.
% %function [yout, yout2, y2] = fft_resample(y, outlen)
% %   returns:  
% %       yout        resampled data
% %       yout2       resampled, subsampled back to original size.  
% %       y2          copy of orig data (fft'd, reverse fft'd.)  S/b identical to orig data...
% %                       yout2 & y2 *SHOULD* be identical to input y, within rounding errors.
% %       you can safely ignore yout2 and y2, unless you're checking to see how much you've lost in rounding.
% %
% 
%     if (~exist('dim','var') || isempty(dim))
%         dim = find(size(y)>1,1);
%     end
%     n=size(y,dim);
% 
%     [nrows, ncols] = size(y);
% 
%     if (outlen==nrows)      % nothing to do.  copy input & return.
%         if (flip)
%             yout=y';
%         else
%             yout=y;
%         end        
%         return;
%     end
% 
%     sizeout =size(y);
%     sizeout(dim) = outlen;
%     
%     yout=nan(sizeout);
% 
%             % might be able to do this as a single step.  try it, ian!
%      
%     Y = fft(y, n, dim);
%     Yout = zeros(sizeout);
%     
% %      YY=fft(y,outlen);
% %      yy_out=ifft(YY) * outlen/nrows;
% 
%     if (dim == 1)
%         Yout = Y(
%     else
%     end
%     
%     for i=1:ncols
%         Y = fft(y(:,i));
%         y2(:,i) = ifft(Y);
%         if (outlen < nrows)           % reducing.  keeping only low-freq terms
%             m1 = fix((outlen+1)/2);
%             m2 = m1 + nrows - outlen;
%             Y = [Y(1:m1);Y(m2:end)];
%         else                        % expanding.  zero-pad FFT        
%             if (mod(nrows,2)==1)   
%                 m1 = fix((nrows+1)/2);
%                 m2=m1+1;
%                 nz=outlen-nrows;
%             else
%                 m1=fix((nrows+1)/2)+1;
%                 m2=m1;           % if input length is even, we need to duplicate midpoint to keep fft conjugate.
%                 nz  = outlen-nrows-1;
%             end
%             Y = [Y(1:m1);zeros(nz,1);Y(m2:end)];
%         end 
%                 % reverse FFT & rescale
%         yout(:,i) = ifft(Y) * outlen/nrows;
%         ix=fix((0:(nrows-1))*outlen/nrows)+1;
%         yout2(:,i) = yout(ix,i);  % subsampled back to original size
%     end
%             
%     if (flip)
%         yout=yout'; % transpose back to col vector
%         yout2=yout2';
%         y2=y2';
%     end
% end
