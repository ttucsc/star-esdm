function [ X ] = gauss( n, sigma, mu )
% [ X ] = gauss( n, sigma, mu )  returns gaussian of size n,
%   n is size of gaussian to generate
%   sigma is sigma of gaussian.  sigma=0 produces an impulse function
%   mu is mean.  if mu not provided, mu is midpoint.
%   results are normalized so sum(X)=1;
%
%   note:    if n is even and sigma=0, then impulse isn't centered.
%   note 2:  if mu is not near the middle and sigma is large relative to length n, 
%            resulting then one tail will be truncated, and though normalized to 1,
%            probably is not what you want...
%
%   ian scott-fleming, 9/2010
%
%   This could be calculated with:  X = normpdf(x,mu, sigma);
%   except normpdf has a problem if sigma is 0.
%   Also, for normpdf, if sigma < stepsize in x, then area does not sum to 1.
 
    n=fix(double(n));
    if (sigma==0.0)
        X = zeros(1,n);
        X(1, fix(n/2)+1)=1;        % note:  even-sized kernel will shift image!
        return;
    end

    x = double(1:n);
    if (nargin < 3)
        mu = double((x(1)+x(n))/2.0);
    else
        mu = double(mu);
    end
            % full equation:
    X = 1.0/(sqrt(2*pi)*sigma) * exp(-((x-mu).^2)/(2.0*double(sigma)^2));
            % The above *should* sum to 1.
            % But with finite precision rounding, the above is generally off slightly, so the
            % area is not exactly 1 (typically off by around 1e-17).  
            % Thus, we normalize it:
    X = X / sum(X);         % normalize to total area of 1.  

end

