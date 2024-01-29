function s = jc_sigmoid(x, k, p)
%   Returns sigmoid on x, approaching -1 for negative values of x, and 1 for positive values.
%   For k==1, sigmoid(1) = .5.
%   higher values of k create steeper sigmoid.
%   if p is 'abs', then output is absolute value ('the bird')  ----v----
%   if p is numeric, then it should be a positive integer, and will raise the sigmoid ^p. ____/----
%   if k and p are 1 (or not provided), sigmoid converges to zero faster than x.
    if (nargin==1), k=1; end
    s = 2/pi * atan(k*x);
    if (nargin == 3)
        if (strcmp(p,'abs'))
            s=abs(s);
        else
            s=s.^p; 
        end
    end
end