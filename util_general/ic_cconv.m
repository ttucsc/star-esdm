function a = ic_cconv(A, b, dim)
%   circularly convolves every row/column/whatever (specified by dim) of A with vector b.
%   if dim is not specified, works along first non-singleton dimension of A.
%
%   Note:  if A is a row or column vector, you might as well use cconv(...).

    sz = size(A);
    if (~exist('dim','var'))
        dim = find(sz>1,1);
    end
    
    n = sz(dim);
    if (length(b) < n)
        b(end+1:n)=0;
    elseif (length(b) > sz(dim))
        b=b(1:n);
    end
   
    a = ifft(fft(A,n,dim) .* fft(b),n,dim);
    
    if (isreal(A) && isreal(b) && ~isreal(a)), a = real(a); end       % in case fft stuff left a with a small residual complex part.
end
    

