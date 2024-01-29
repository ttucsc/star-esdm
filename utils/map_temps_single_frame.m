function map_temps_single_frame(fname, varname, lbl, varargin)
% function map_temps_single_frame(fname, varname, lbl, varargin)
%                       fname,      netcdf filename to read, or matname to read
%                       varname,    variable to animate from file
%                       lbl,        text description for figure title
%   Optional keyword/value pairs:
%       "is_us",        is_us,      if true, draws state boundaries.  if false, draws only continent and island boundaries 
%       "dates",        yr_range,   years to display.
%       "latrange",     latrange,   region to animate
%       "lonrange",     lonrange
%       "mode",         mode,       0,1,2,3  animation type (see below)
%       "filterwidth",  kernelwidth, # of days to average into each frame.
%       "minmax",       minmax      min & max vals to use for range instead of calculating.
%       "clim_yrs",     [yr1,yr2]   years to use for climatology;  Usually only 1950-2010 if the data goes to 2100.
%                                       if empty, uses all years to calculate the climatology
%       "dir",          outdir      output folder name (ignored if fname is fully qualified)  [.]
%       "fignum",       fignum      figure number or handle to draw into [create new figure]
%       "figname",      figname     name of file to save figure to [do not save fig]
%       "lats_low",     lats_low    low-res lats for drawing lat/lon grid.  Usually from original model gridding 
%       "lons_low",     lons_low    low-res lons for drawing lat/lon grid.  
%
%       mode        0       display temps or precip values as is
%                   1       difference from overall average.            Will show winter & summer, but not latitude-dependent differences
%                   2       difference from previous day avg or sum     high-pass filter in time
%                   3       * difference from mean climatology          (averaged) anomaly.  Will show weather patterns moving
% %                 4       * ** difference from dynamic climatology    (not yet implemented.  Will required doing an
%                                                                           FFT filtering on each gridpoint.
%                               * these take a while to compute.
%                               ** not implemented yet.

%   note on lon specification:  lonrange must be in the same range as lons in the file:  i.e. -180 to 180 or 0-360.

    [fname, varname, lbl, is_us, mode, year_range, dates, kernelwidth, isPrecip, isTemps, usemat,...
            latrange, lonrange, lats, lons, ...
            datevecs, cal, verbose, minmax, fignum, clim_yrs, fullrange, figname, lats_low, lons_low] = init_params(fname, varname, lbl, varargin{:});
    if (is_us)
        border_flag=2;
    else
        border_flag=1;
    end
    

    if (usemat)
          
        fprintf("reading data from matfile %s\n", fname);
        [v, vmax, vmin, datevecs, lats, lons, filename] = my_matload(fname, varname, lats, lons, datevecs, "day-365", latrange, lonrange, year_range,fullrange);         
        fprintf("data originally from %s\n", filename)
        
        if (~isempty(minmax))
            vmin = minmax(1);
            vmax = minmax(2);
        end
        
    else
        
        if (calendar_length(cal) ~= 365)
            fprintf("changing calendar to 365-day from %s\n", cal);
            cal = "365-day";
        end
        
        [v, ~, vmin, datevecs, lats, lons] = my_ncload(fname, varname, latrange, lonrange, year_range, fullrange, verbose);

        datevecs = datevecs(:,1:3);
        
        if (isTemps)
            if (vmin > 150)
                v = v - 273.15;
            end
%             if (vmin > 150)
%                 vmax = vmax - 273.15;
%                 vmin = vmin - 273.15;
%             end                
        end
        
%       mode        0       display temps or precip values as is
%                   1       difference from overall average
%                   2       difference from previous day avg or sum
%                   3       * difference from mean climatology
%                   4       * ** difference from dynamic climatology
%                               * these take a while to compute.
%                               ** not implemented yet.
        if (mode == 1)
            v = v - nanmean(v);
        elseif (mode == 3 || mode == 4)
            fprintf("removing climatology...");
            [v, ~, ~]  = remove_climatology(v, clim_yrs, mode);
            fprintf("done\n");
        end
        
        if (~isempty(kernelwidth) && kernelwidth > 1)        %!!!   
            fprintf("low-pass-filtering data...");
            v = rect_filt(v, kernelwidth,1, false);
            fprintf("done");
        end
        
        if (mode == 2)      % difference from one step to the next
            v = diff(v);
        end
            
        if (~isempty(minmax))
            vmin = minmax(1);
            vmax = minmax(2);
        else
            vmax = nanmax(v(:));
            vmin = nanmin(v(:));
        end        
        
    end
    
    nmaps = size(dates,1);
        
    cmap = get_cmap(isPrecip, kernelwidth);
%   colormap(cmap);
    R =  setup_map(lats, lons);     

    dnums = datenum_cal(datevecs, cal);
    
%   h = figure(fignum);

    if (nmaps > 1)
        nrows = ceil(sqrt(nmaps));
        ncols = ceil(nmaps/nrows);
    else
        subpos = [];
    end
    
    for i=1:nmaps
        if (nmaps > 1)
            subpos = [nrows, ncols, i];
        end
        
        dnum = datenum_cal(dates(i,:),cal);
        ix = dnum - dnums(1)+1;
        if (isPrecip)
            [img, imvmin, imvmax] = precip_scale(squeeze(v(ix,:,:)), vmin, vmax, true);
        else
            img = squeeze(v(ix,:,:));
            imvmin = vmin;
            imvmax = vmax;
        end
                
        ttl = sprintf("%s %s",lbl,datestr_cal(dnum,cal, "yyyy-mm-dd"));
        h = display_map(img,ttl,R,fignum,subpos,imvmin, imvmax, border_flag, cmap, lats_low, lons_low);
        figure(h);  % bring it to the front...
    end

    if (~isempty(figname))
        fprintf("saving figure as %s...", figname);
        savefig(h, figname);
        fprintf("done\n");
    end
    figure(h);
end

function  R = setup_map(lats, lons)
    
        % now display movie of v
    nlats = length(lats);
    nlons = length(lons);
    if(lats(1)==-90)
        dlat=0;
    else
        dlat=lats(5)-lats(4);
    end
    if (lats(end)-lats(1) > 175)
        lats(1)=-90;
        lats(end)=90;
    end
    if (lons(end)-lons(1) > 355)
        lons(1)=0;
        lons(end)=360;
    end
    
    if (lons(1)==0 || lons(1)==-180)
        dlon=0;
    else
        dlon = lons(5)-lons(4);
    end
%   R = make_refmat(lats, lons, nlats,nlons,lats(10) - lats(9),lons(10) - lons(9))  ; 
    R = make_refmat(lats, lons, nlats,nlons,dlat,dlon)  ; 
end

function [img, imvmin, imvmax] = precip_scale(img, vmin, vmax, dolog)

    nanmap = isnan(img);
    img(nanmap)=0;
    imvmin = log10(1+vmin);
    imvmax = log10(1+vmax);
    if (dolog)
        img = log10(1+img);     % vmin and vmax are never log scale yet, but img is in 2nd half of program.
    end    
    img = min(img, imvmax);

    
end
    
function     [fname, varname, lbl, is_us, mode, year_range, dates, kernelwidth, isPrecip, isTemps, usemat,...
            latrange, lonrange, lats, lons, ...
            datevecs, cal, verbose, minmax, fignum, clim_yrs, fullrange, figname, lats_low, lons_low] = init_params(fname, varname, lbl, varargin)


    varnames = ["tasmax","tmax","tasmin","tmin","tas","tavg","pr","prcp","prec","precip","rh","relhum"];
    tempnames =["tasmax","tmax","tasmin","tmin","tas","tavg"];
    prcpnames =["pr","prcp","prec","precip"];
                    % parse input
    p = inputParser;

    
    addRequired(p,"fname",                  @(s) ischar_s(s) || isa(fname, "ncdf"));
    addRequired(p,"varname",                @(s) ischar_s(s) && any(strcmpi(s,varnames)));
    addRequired(p,"lbl",                    @(s) ischar_s(s));   
    addParameter(p,"latrange",  [],         @(s) isempty(s)   || (isnumeric(s) && length(s)==2 && s(1)>= -90 && s(2) >= s(1) && s(2) <= 90));
    addParameter(p,"lonrange",  [],         @(s) isempty(s)   || (isnumeric(s) && length(s)==2 && s(1)>= -180 && s(2) >= s(1) && s(2) <=360 && s(2)-s(1) <= 360));
    addParameter(p,"is_us",     [],         @(s) isempty(s)   || islogical(s));
    addParameter(p,"dates",     [],         @(s) isnumeric(s) && size(s,2)==3 && s(1) > 1850 && s(1) < 2100);
    addParameter(p,"mode",      0,          @(s) isnumeric(s) && s>= 0 && s <= 5);
    addParameter(p,"filterwidth",[],        @(s) isempty(s)   || (isnumeric(s) && s>0 && s <= 365));
    addParameter(p,"dir",       [],         @(s) isempty_s(s) || isfolder(s));
    addParameter(p,"verbose",   true,       @(s) islogical(s) || (s==0 || s==1));
    addParameter(p,"fignum",    1,          @(s) isnumeric(s) && s > 0);
    addParameter(p,"minmax",    [],         @(s) isnumeric(s) && (isempty(s) || length(s) == 2));
    addParameter(p,"clim_yrs",  [],         @(s) isnumeric(s) && (isempty(s) || length(s) == 2));
    addParameter(p,"figname",   [],         @(s) isempty_s(s) || ischar_s(s));
    addParameter(p,"lats_low",  [],         @(s) isempty(s)   || isnumeric(s));
    addParameter(p,"lons_low",  [],         @(s) isempty(s)   || isnumeric(s));
            
    parse(p, fname, varname, lbl, varargin{:});

    latrange = p.Results.latrange;
    lonrange = p.Results.lonrange;
    fname       = p.Results.fname;
    varname     = p.Results.varname;
    lbl         = p.Results.lbl;
    mode        = p.Results.mode;
    dates       = p.Results.dates;
    is_us       = p.Results.is_us;
    kernelwidth = p.Results.filterwidth;
    basedir     = p.Results.dir;
    verbose     = p.Results.verbose;
    minmax      = p.Results.minmax;
    fignum      = p.Results.fignum;
    clim_yrs    = p.Results.clim_yrs;
    figname     = p.Results.figname;
    lats_low    = p.Results.lats_low;
    lons_low    = p.Results.lons_low;
    
        % add base folder if fname or figname aren't fully qualified.
    if (~isempty_s(basedir))
        if (ischar_s(fname) && ~strcmp(extractBefore(fname,2),"/"))
            fname=fullfile(basedir,fname); 
        end
        if (~isempty_s(figname) && strcmp(extractBefore(figname,2),"/"))
            figname=fullfile(basedir,figname); 
        end
    end
    if (ischar_s(fname)) 
        if (~isfile(fname)), error("error:  file %s does not exist or is not readable", fname); end
        [~,~,ext] = fileparts(fname);
        usemat = strcmpi(ext,".mat");
    else
        usemat = false;
    end
    if (usemat)
        T = load(fname,   "stepsize", "kernelwidth","lats","lons","cal","datevecs");
        cal = T.cal;
        datevecs = T.datevecs;
        tstamps = datenum_cal(datevecs,cal);
        if (~isempty(kernelwidth))
            if (T.kernelwidth ~= kernelwidth), error("specified kernelwidth (%f) doesn't match matfile's kernelwidth (%f).  Please use .nc file instead of matfile", kernelwidth, T.kernelwidtth); end
        else
            kernelwidth = T.kernelwidth;
        end
        lats = T.lats;
        lons = T.lons;
    else

        if (isempty(kernelwidth))
            kernelwidth = 1;
        end
        [lats,lons] = ncdf_get_latlons(fname);
        [~,cal, ~,~,~,tstamps] = ncdf_get_time_info(fname,varname);
        datevecs = datevec_cal(tstamps, cal);
    end

    tstamps = floor(tstamps);
    file_tstamprange = [tstamps(1),tstamps(end)];
    file_latrange    = [lats(1), lats(end)];
    file_lonrange    = [lons(1), lons(end)];
    
    isPrecip = any(strcmpi(varname,prcpnames));
    isTemps  = any(strcmpi(varname,tempnames));
   
        % get latrange, lonrange
    
    if (isempty(latrange))
        latrange = file_latrange;
    else
        latrange = [max(latrange(1),file_latrange(1)),min(latrange(2),file_latrange(2))];
    end
    if (isempty(lonrange))
        lonrange = file_lonrange;
    else
        lonrange = [max(lonrange(1), file_lonrange(1)), min(lonrange(2),file_lonrange(2))];
    end
    if (isempty(is_us))
        is_us = ((lonrange(1) >= -180 && lonrange(2) < -30) || (lonrange(1) >= 180 && lonrange(2) <= 330)) && ...
                 (latrange(1) >=    0 && latrange(2) <= 90);
    end
    
    if (~isempty(lats_low))
        keepers = lats_low >= latrange(1) & lats_low <= latrange(end);
        lats_low = lats_low(keepers);
    end
    if (~isempty(lons_low))
        keepers = lons_low >= lonrange(1) & lons_low <= lonrange(end);
        lons_low = lons_low(keepers);
    end
    
        % get dates to display and range of dates to pull.
    
    if (isempty(dates))
        dnums = tstamps(1) + floor(kernel_width/2);
    else
        dnums = datenum_cal(dates, cal);
    end
    
    dates = datevec_cal(dnums, cal);
    if (isempty(clim_yrs) || mode < 3)
        year_range = [min(dates(:,1)), max(dates(:,1))];
    else
        year_range = [min(clim_yrs(1), dates(1,1)), max(clim_yrs(2), dates(end,1))];
    end   
    
    date_range = [year_range(1),1,1;  year_range(end),12,31-1*(calendar_length(cal)==360)];

    tstamp_range = [datenum_cal(date_range(1,:),cal), datenum_cal(date_range(2,:),cal)];
    
    tstamp_range = [max(tstamp_range(1),file_tstamprange(1)), min(tstamp_range(2),file_tstamprange(end))];
    
    if (any(dnums < tstamp_range(1)) || any(dnums > tstamp_range(end)))
        error("error:  date(s) are outside range of dates in file");
    end
    
    if (all(tstamp_range == file_tstamprange) && all(latrange == file_latrange) && all(lonrange == file_lonrange))
        fullrange = true;
    else
        fullrange = false;
    end
    
end

function [v, vmax, vmin, datevecs, lats, lons] = my_ncload(fname, varname, latrange, lonrange, year_range, fullrange, vbose)

    my_vbose = ncObj.verbose(vbose);
    
    if (ischar_s(fname))
        [~,fn,fext] = fileparts(fname);
        myfname = strcat(fn,fext);
        fprintf("reading netcdf file %s\n", myfname);
    end

    cal = '365-day';
    if (fullrange)
        nc = ncdf_read_file(fname, varname, [], [], [], cal);
    else
        nc = ncdf_read_file(fname, varname, latrange, lonrange, year_range, cal,false, varname, false, -1);
    end
    fprintf("done reading\n");
    v = nc.getvardata(varname);
    vmax = nanmax(v,[], 'all');
    vmin = nanmin(v,[], 'all');
    
    lats = nc.getvardata("lat");
    lons = nc.getvardata("lon");
    
    [~,mycal, ~, ~, ~, tstamps] = ncdf_get_time_info(nc);
    datevecs = datevec_cal(tstamps, mycal);
    
    ncObj.verbose(my_vbose);

end

function [v, vmax, vmin, datevecs, lats, lons, filename] = my_matload(fname, varname, lats, lons, datevecs, cal, latrange,lonrange,year_range,fullrange)

    matobj = matfile(fname);
    
    vname = matobj.vname;
    if (~strcmp(vname, varname)), error("error:  varname %s doesn't match varname %s in matfile %s", varname, vname, fname); end
    vmin = matobj.vmin;
    vmax = matobj.vmax;
    isPrecip = matobj.isPrecip;
    filename = matobj.filename;
    fprintf("reading %s\n", varname);
    if (fullrange)
        v = matobj.(vname);        
        daterange = [matobj.datevecs(1,1:3); matobj.datevecs(end,1:3)];
    else
        daterange = [year_range(1),1,1; year_range(2),12,31-1*(calendar_length(cal)<365)];
        [timix1, timix2, latix1, latix2, lonix1, lonix2] = range_indexes(lats, lons, datevecs, cal, latrange, lonrange, daterange);    
        v = matobj.(vname)(timix1:timix2, latix1:latix2, lonix1:lonix2);
    end
    fprintf("done reading\n");
    
    datevecs = datevec_cal(datenum_cal(daterange(1,:),cal):datenum_cal(daterange(2,:),cal),cal);

    mynans = v==65535;
    v = single(v) / single(65534) * single((vmax - vmin)) + single(vmin);
    v(mynans) = nan;
    if (isPrecip)
        vmax = 10^vmax-1;
        vmin = 10^vmin-1;
    end
          
end
    
function [timix1, timix2, latix1, latix2, lonix1, lonix2] = range_indexes(lats, lons, datevecs, cal, latrange, lonrange, daterange)

    tstamps = datenum_cal(datevecs, cal);
    trange  = datenum_cal(daterange, cal);
    got_times=false;
    got_lats = false;
%   got_lons = false;
    try
        timix1 = find(tstamps >= trange(1),1);
        timix2 = find(tstamps <= trange(2),1,'last');
        got_times = true;
        
        latix1 = find(lats >= latrange(1),1);
        latix2 = find(lats <= latrange(2),1,'last');
        got_lats = true;

        lonix1 = find(lons >= lonrange(1), 1);
        lonix2 = find(lons <= lonrange(2), 1, 'last');
    catch
        if (~got_times)
            error("no data in specified date range %s - %s.  File:  %s - %s", datestr_cal(datevecs), datestr_cal(tstamps(1)), datestr_cal(tstamps(end)));
        elseif (~got_lats)
            error("no data in specified lat range %.4f - %.4f.  File:  %.4f - %.4f", latrange, lats(1), lats(end));
        else
            error("no data in specified lon range %.4f - %.4f.  File:  %.4f = %.4f", lonrange, lons(1), lons(end));
        end
    end
end
    
function cmap = get_cmap(isPrecip, kernelwidth)

    cmap = jet(256);
    if (isPrecip)
        cmap(120:129,:) = 1;
        cmap(256,:) = [.25 .75 1];
            % set max to 10 inches per week average.
        vmax =  log10(kernelwidth * 2540);
        daymax = log10(kernelwidth * 10/7*25.4);
        ix = ceil((daymax+vmax)/(2*vmax)*256);
        n=256-ix+1;
        cmap(ix:end,:)=repmat([.25,.75,1],n,1);
    end    
end

function [vout, vmin, vmax] = remove_climatology(v, clim_yrs, mode)

    vclim = calc_climatology(v, clim_yrs, mode);
    
    ndays = size(v,1);
    nyrs = ndays/365;
    if (mode == 3)
        vout = v - repmat(vclim,nyrs,1,1);
    else
        vout = v - vclim;
    end

    vmin = nanmin(vout(:));
    vmax = nanmax(vout(:));
end
function vclim = calc_climatology(v, clim_yrs, mode)
    % does fft-based low pass filter of 365-day daily average of data, using data from yr1 to yr1.  
    %   v s/b of size ndays x nlats x nlons.
    %       ndays must be even multiple of 365-days long.
    %   yr1 & yr2 are ordinal, not calendar years:  1 is 1st year of data, so v,1,10  uses average of first 10 365-day years.
    %       leave yr1 empty to start at beginning.  
    %       leave yr2 empty to end at the end
    
    if (mode ~= 3), error("error:  mode %d not implemented yet\n", mode); end

    [ndays,nlats,nlons] = size(v);
    nyrs = ndays/365;
    if (mod(nyrs,1) ~= 0), error("error: data is not even multiple of 365 days: %d days (s/b %d days)", ndays, floor(nyrs)*365);  end
    if (isempty(clim_yrs) || (clim_yrs(1)==1 && clim_yrs(2)==nyrs))
        v =  reshape(v,365,[],nlats,nlons);
    else
        v = reshape(v(365*(clim_yrs(1)-1)+1:365*clim_yrs(2),:,:),365,[],nlats,nlons);
    end 
%   clim1 = nanmean(v,2);
%   C = fft(clim1);
    
    FILT = cast(calc_filter(6, 2, 1, 365)','like',v);

        % just keep the real part of the fft-filtered data.  ideally it should be purely real;  limited-precision math
        % usually leaves a very small imaginary part.

    vclim = real(ifft(fft(squeeze(nanmean(v, 2))) .* repmat(FILT,1,nlats, nlons)));
    
end

