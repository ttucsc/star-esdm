function weights = calc_weights(lat_pt, lon_pt, lats, lons, interp_method)
%       quad_lats  = [latout(1); latout(1); latout(2); latout(2)];
%       quad_lons  = [lonout(1); lonout(2); lonout(1); lonout(2)];

        % if we match lat or lon exactly, then lats or lons will only have one element.  We need two here.
    if (length(lats)==1)
        lats = [lats-.0001, lats+.0001];
    end
    if (length(lons)==1)
        lons = [lons-.0001, lons+.0001];
    end
    quad_lats  = [lats(1); lats(2); lats(1); lats(2)];
    quad_lons  = [lons(1); lons(1); lons(2); lons(2)];
    lat_pts = repmat(lat_pt,4,1);
    lon_pts = repmat(lon_pt,4,1);
    dists = distance(quad_lats, quad_lons, lat_pts, lon_pts);
    dinv = 1./(dists.^2);   % weights based on 1/square(distance)

    if (strncmpi(interp_method,'bil',3))              % 'bilinear_sampling':  use standard ARRM_V2 sampling with bilinear weights.
        weights = bilinear_weights_ic(lat_pt, lon_pt, lats, lons);            

            %  special case:  latpt,lonpt is exactly on one of the gridcell centers.
            %  Set weight for that gridcell to almost 1, and set tiny weights for the others based on inverse
            %  distance. This allows us to use alternate sampling when main gridpoint is NA.
            %   This is same as zero-distance case in inverse_distance weighting below.
        if (any(weights==1.0))
            jx=find(dists==0,1);
            dinv(jx) = nan;
            weights = 1e-8*dinv/sum(dinv,"omitnan");
            weights(jx) = 1.0 - sum(weights,"omitnan");
        end

    elseif (strncmpi(interp_method,'inv',3))            % 'inverse_distance':  use standard ARRM_V2 sampling with inverse_distance weights.
            %  special case:  latpt,lonpt is exactly on one of the gridcell centers.
            %  Set weight for that gridcell to almost 1, and set tiny weights for the others
            %   This allows us to use alternate sampling when main gridpoint is NA.

        if (any(char(interp_method)=='1'))
            pwr=1;
        else
            pwr=2;
        end
        if (any(dists==0))  
            jx=find(dists==0,1);
            dinv = 1./(dists.^pwr);   % weights based on 1/square(distance)
            dinv(jx) = nan;
            weights = 1e-8*dinv/sum(dinv,"omitnan");
            weights(jx) = 1 - sum(weights,"omitnan");
        else
            dinv = 1./(dists.^pwr);   % weights based on 1/square(distance)
            weights = dinv/sum(dinv);
        end
        weights = reshape(weights,2,2);
    else
        error('error:  %s:  unknown sampling method: %s', mfilename, interp_method);
    end

end
