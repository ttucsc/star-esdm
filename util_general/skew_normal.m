function [y,sn_stats, sn_cdf] = skew_normal(alpha,x,mu, sig, do_ctr)
%returns skew normal(alpha) for a gaussian with mean mu, std dev sig, calculated over range x.
% See https://en.wikipedia.org/wiki/Skew_normal_distribution
%
%   Inputs:
%     Required:
%       alpha                       shape factor.  alpha < 0 -> left-skewed  alpha > 0, right-skewed
%                                       alpha = 0 -> std. normal curve.
%     Optional: [defaults]
%       x       [-10:.1:10]         independent variable range
%       mu      [mid-point of x]    mean for initial gaussian
%       sig     [1.0]               sigma for initial gaussian
%       do_ctr  [false]             if true, repositions output so mean is set @ mu.
%                                       if false, mean will be offset to left (negative alpha) or right (positive alpha)
%   Return values
%       y           skewed gaussian, adjusted to sum to 1.0
%                       note:  if range is inappropriate, results will still be adjusted to sum to 1, but shape of
%                       distribution will not be as expected.
%       sn_stats    skewed-normal stats struct, with fields mu, sigma, skew and kurt.
%                               note:  xkurt is excess kurtosis (kurtosis-3)  xkurt==0 for pure gaussian (alpha=0)
%       sn_cdf      CDF, integral of y, ( cumsum(y) ) normalized to sk_cdf(end) == 1
%                       As with y, the pdf, the shape will not be as expected if range of x is inappropriate.  
    
    if (~exist('do_ctr','var') || isempty(do_ctr)), do_ctr = false;     end
    if (~exist('sig',   'var') || isempty(sig)),    sig = 1.0;          end
    if (~exist('mu',    'var') || isempty(mu)),     mu  = 0.0;          end
    if (~exist('x',     'var') || isempty(x)),      x   = [-10,.1,10];  end
        
    xrange  = range(x);    
    len     = length(x);
 %   xstep   = xrange/(len-1);
    sigsteps= (len-1)/(xrange/sig);        % steps per std. deviation.
    dx = xrange/(len-1);

    xxmin = min(x(1), -10*sig+mu);
    xxmax = max(x(end),10*sig+mu);
    myrange = range([xxmin,xxmax]);
    
    mylen = floor(myrange*sigsteps);
    if (mod(mylen,2)==0)
        mylen=mylen+1; 
    end
    xxmax = xxmin + (mylen-1)*dx;
    myx = linspace(xxmin, xxmax, mylen);

        %gaussian w/ area of 1, with std. deviation of sigsteps, length pf mylen, and centered.
    g = gauss(mylen, sigsteps) / dx;
    phi = .5*(1+erf(alpha*((myx-mu)/sig)/sqrt(2)));    
    yy = 2*g .* phi;
    if (do_ctr)
        gmu = pdf_mus(yy, myx,true) - mu;
        myx = myx - gmu;
    end
    
    y = interp1(myx,yy,x,'makima');
    y = y / sum(y) / dx;     % normalize to be a true pdf

    if (nargout > 1)
        sn_stats = struct;
        [sn_stats.mu, sn_stats.sigma, sn_stats.skew, sn_stats.xkurt] = pdf_stats(y,x,[],true);
        if (nargout > 2)
            sn_cdf = cumsum(y);
            sn_cdf = sn_cdf / sn_cdf(end);
        end
    end

end

