function  [vdata,nc]  = ncdf_read_gridded_data(fnames, varname, lat_pt, lon_pt, daterange, outcalendar, lbl)
% [vdata,nc]  = ncdf_read_gridded_data(fnames, varname, lat_pt, lon_pt, daterange, outcalendar, lbl)
% function to read closest four gridpoints of data from a netcdf file and interpolate the results to the specified lat/lon.
%   This should duplicate the input used by ARRM_V2 for the same lat/lon
%
%   data returned brackets the latrange & lonrange specified, so if a single point is specified, will generally return
%   the four surrounding gridpoints from the input data.
%       Use ncdf_read_nearest_single_location(...) to get the data for the single closest data point.
% Returns data in netcdf order (lon, lat, time)
%   which is the matlab orientation of indata(time, lat, lon).
%
%   Note:   this will be REAAAALLLLY slow if file is not organized with data contiguous in time 
%           i.e., like CSC's *.llt.nc files, or for netcdf4 files, if chunked sensibly 
%
%   Assumes:  lats & lons are identical between files.  
%
%   Required Inputs:
%           fnames          filenames to read (or ncdf objects with (unmodified) meta-data of a previously opened file
%           varName         variable to read from file
%           lat_pt,lon_pt   location to find data for.  S/B single lat & lon.
%                               if you need range of lats & lons, use ncdf_read_files(...)
%   Optional Inputs:
%           daterange       either [yr1, yr2] or [yr1, mo1, day1; yr2, mo2, day2]  specifying start and end date range
%                               if empty, then reads all dates
%           outcalendar     optional.  can change calendar if desired.  if empty, calendar is not changed.
%           lbl             string,  If lbl is not empty, it is printed, along with the lats & lons of the gridpoints and their weights.
%                               useful for getting data for specific city or station.
%
%   Returns:
%           vdata           The interpolated data for the lat/lon pt and dates specified.
%           nc              ncdf object containing appropriate metadata, including lat, lon, time, & vardata. 
%
%       NOTE:  vdata returned first.  ncdf_read_files returns [nc, vdata], in reversed order to this function.  
%
%   10/21/2021  icsf  modified to return vdata and nc.  Also updated header comments.
%
%-------------------------------------------------
 

    if (~exist("daterange","var"))
        daterange = [];
    end
    if (~exist("outcalendar","var"))
        outcalendar = strings(0);
    end
    if (~exist("lbl","var") || isempty(lbl))
        lbl = strings(0); 
    else
        lbl = string(lbl);
    end


    if (ischar(fnames))
        fnames = string(fnames);
    end
    if (numel(fnames)==1)
        nc = ncdf_read_file(fnames, varname, lat_pt, lon_pt, daterange, outcalendar);
    else
        nc = ncdf_read_files(fnames, varname, lat_pt, lon_pt, daterange, outcalendar);
    end

    [lats, lons] = ncdf_get_latlons(nc);

    all_mdldata = nc.getvardata(varname);

    if (numel(lats)==1 && numel(lons)==1)   % exact match to lat/lon location.
        weights = 1;
        ix = 1;
    else                                    % need to interpolate.      
        if (length(lats)==1)                % make sure we have 4 points.  If exact match on lat or lon, duplicate that data so we have 4 points to calculate the weights..                           
            lats = [lats-.0001, lats+.0001];
            all_mdldata = repmat(all_mdldata,1,2,1);
        elseif (length(lons)==1)
            lons = [lons-.0001, lons+.0001];
            all_mdldata = repmat(all_mdldata,1,1,2);
        end    
        all_mdldata = reshape(all_mdldata,[],4,1);
            
        weights = calc_weights(lat_pt, lon_pt, lats, lons, "bilinear");        
        weights = reshape(weights,4,1);
        vdata = weighted_sum(all_mdldata, weights);
        
        ix = [1, 1; 2,1; 1,2; 2,2];
    end

    if (strlength(lbl)>0)
        if (size(ix,1)==1)
            fprintf("%s:  (exact match to location):", lbl)
            fprintf("   lat       lon    weight\n");
            fprintf("%9.4f %9.4f %6.3f\n", lats, mod(lons+180,360)-180, weights);
        else
            fprintf("%s Grid Points:\n", lbl);
            fprintf("   lat       lon    weight\n");
            for i=1:4
                fprintf("%9.4f %9.4f %6.3f\n", lats(ix(i,1)), mod(lons(ix(i,2))+180,360)-180, weights(i));
            end
        end
    end
    
end

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

