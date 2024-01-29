function gval = gaussval(x, sigma, mu)

    if (nargin<2), sigma=1; end
    if (nargin<3), mu=0;
    gval = exp(-((x-mu).^2/(2.0*double(sigma)^2)));
end

