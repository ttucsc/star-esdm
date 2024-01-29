function animate_temps(fname, varname, lbl, varargin)
% function animate_temps(
%                       fname,      netcdf filename to read, or matname to read
%                       varname,    variable to animate from file
%                       lbl,        text description for figure title
%   Optional keyword/value pairs:
%       "is_us",        is_us,      if true, draws state boundaries.  if false, draws only continent and island boundaries 
%       "year_range",   yr_range,   years to display.
%       "pauselen",     pauselen,   in seconds, between frames
%       "latrange",     latrange,   region to animate
%       "lonrange",     lonrange
%       "mode",         mode,       0,1,2,3  animation type (see below)
%       "stepsize",     stepsize,   time period between frames, in days
%       "filterwidth",  kernelwidth, # of days to average into each frame.
%       "minmaxrange",  minmax      min & max vals to use for range instead of calculating.
%       "clim_yrs",     [yr1,yr2]   years to use for climatology;  Usually only 1950-2015 if the datagoes to 2100.
%                                       if empty, uses all years to calculate the climatology
%       "dir",          basedir     folder name where data is read from, and where mat file is written to
%                                       (ignored if fname or matname is fully qualified)  [.]
%       "figname",      figname     name of file to write output to.
%                                       should have extension ".mpeg", ".mp2" or ".avi".
%
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
%
%   updated 10/21;  replaced nanmean,nanmin & nanmax with updated calls to mean, min & max.
%                   replaced drawing to the screen with writing to a video file.
%
%   STILL NEEDS:  spatial downsampling?
%                 currently reads the entire file at once.
%                 should read only the timeslices needed, work on them, then get the next
%                 Ability to use livneh or nclimgrid (original gridded obs) for baseline instead of self.

    [fname, varname, lbl, is_us, mode, year_range, frameRate, stepsize, kernelwidth, isPrecip, isTemps, ...
                                        latrange, lonrange, lats, lons, datevecs, cal, fullrange, verbose, minmax, fignum, clim_yrs] = init_params(fname, varname, lbl, varargin{:}); %#ok<ASGLU>
    if (is_us)
        border_flag=2;
    else
        border_flag=1;
    end
    

%     if (usemat)
%           
%         fprintf("reading data from matfile %s\n", fname);
%         [v, vmax, vmin, datevecs, lats, lons, filename] = my_matload(fname, varname, lats, lons, datevecs, "day-365", latrange, lonrange, year_range,fullrange);         
%         fprintf("data originally from %s\n", filename)
%         
%         if (~isempty(minmax))
%             vmin = minmax(1);
%             vmax = minmax(2);
%         end
%         
%     else
        
        if (calendar_length(cal) ~= 365)
            fprintf("changing calendar to 365-day from %s\n", cal);
            cal = "365-day";
        end
        
        [~,~,fext] = fileparts(figname);
        if (strcmpi(fext,".mj2"))
            vid = VideoWriter(figname, "Motion JPEG 2000");
            vid.CompressionRatio = 10;
        elseif (strcmpi(fext,".avi"))
            vid = VideoWriter(figname, "Motion JPEG AVI");
            vid.Quality = 33;
        elseif (strcmpi(fext,".mp4"))
            vid = VideoWriter(figname, "MPEG-4");
            vid.Quality = 33;
        end
%       vid.Colormap=cmap;
        vid.FrameRate=frameRate;
        open(vid);
        
        
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
            v = v - mean(v,"omitnan");
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
        
        if (stepsize ~= 1)
            ndays = size(v,1);
            keepers = ceil(stepsize/2):stepsize:ndays;
            v = v(keepers,:,:);
            datevecs = datevecs(keepers,:);
        end
        
        if (mode == 2)      % difference from one step to the next
            v = diff(v);
        end
            
        if (~isempty(minmax))
            vmin = minmax(1);
            vmax = minmax(2);
        else
            vmax = max(v(:),[],"omitnan");
            vmin = min(v(:),[],"omitnan");
        end        
        
%         if (~isempty_s(matname))
%             dir=fileparts(fname);
%             savemat(fullfile(dir,matname), varname, v, lats, lons, datevecs, cal, stepsize, kernelwidth, fname, isPrecip);
%         end
%   end
    
    fprintf("data range: %6.2f - %6.2f\n", min(v(:),[],"omitnan"), max(v(:),[],"omitnan"));
    
    nsteps = size(v,1);
        
    cmap = get_cmap(isPrecip, kernelwidth);
    colormap(cmap);
    [coastlat, coastlon, R] = setup_map(lats, lons, is_us);     
    if (isPrecip)
        [img, imvmin, imvmax] = precip_scale(squeeze(v(1,:,:)), vmin, vmax, false);
    else
        img = squeeze(v(1,:,:));
        imvmin = vmin;
        imvmax = vmax;
    end
    
    h = first_map(datevecs, cal, lbl, img, R, imvmin,imvmax, border_flag,coastlat,coastlon, cmap, fignum);
    f = getframe(h);
    writeVideo(vid,f);           
   
    for i=2:nsteps     
            
        if (isPrecip)
            [img, imvmin, imvmax] = precip_scale(squeeze(v(i,:,:)),vmin, vmax, false);
        else
            img = squeeze(v(i,:,:));            
        end
            
        update_map(datevecs, cal, i, img, R, lbl, coastlat, coastlon, imvmin, imvmax);
        f = getframe(h);
        writeVideo(vid,f);           
        
        pause(pauselen);

        show_progress(i,nsteps);
    end

    close(vid);
    
end

function [coastlat, coastlon, R] = setup_map(lats, lons, is_us)
    
    latrange = [lats(1),lats(end)];
    lonrange = [lons(1),lons(end)];
    
    if (is_us)
        [coastlat,coastlon, slat, slon] = ncamerica_and_state_boundaries([],false, [], latrange, lonrange);         % bug when lat,lon provided...needs some extra nan's! 
%       
%       [coastlat,coastlon, slat, slon] = ncamerica_and_state_boundaries([],false, []);
        coastlat = [coastlat;nan;slat];
        coastlon = [coastlon;nan;slon];
    else
        load("coastlines","coastlat","coastlon");
    end
    
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

function h = first_map(datevecs, cal, lbl, img,R, imvmin,imvmax, border_flag,coastlat,coastlon, cmap, fignum)

    ismin = img <= imvmin;
    ismax = img >= imvmax;
    if (any(ismin(:)))
        img(ismin) = imvmin; 
    else
        img(end,end) = imvmin;
    end
    if (any(ismax(:)))
        img(ismax) = imvmax;
    else
        img(1,1) = imvmax;
    end
    yr = datevecs(1,1);
    month = datevecs(1, 2);
    doy = floor(datenum_cal(datevecs(1,:), cal) - datenum_cal([datevecs(1,1),1,0],cal));
    ttl = sprintf("%s %4d %s",lbl, yr, get_monlbl(month, doy));
    h = display_map(img,ttl,R,fignum,[],imvmin, imvmax,border_flag,cmap);
    hold on;
    plotm(coastlat,coastlon);
    hold off;
    drawnow();

end

function update_map(datevecs, cal, i, img, R, lbl, coastlat, coastlon, imvmin, imvmax) 

    ismin = img <= imvmin;
    ismax = img >= imvmax;
    if (any(ismin(:)))
        img(ismin) = imvmin; 
    else
        img(end,end) = imvmin;
    end
    if (any(ismax(:)))
        img(ismax) = imvmax;
    else
        img(1,1) = imvmax;
    end
    yr = datevecs(i,1);
    month = datevecs(i, 2);
    doy = floor(datenum_cal(datevecs(i,:), cal) - datenum_cal([datevecs(i,1),1,0],cal));
    meshm(img,R);
    title(sprintf("%s %4d %s",lbl, yr, get_monlbl(month, doy)),'interpreter','none','FontName','FixedWidth',"FontSize",15);
    hold on;
    plotm(coastlat,coastlon);
    hold off;
    drawnow();
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
%   img(nanmap) = imvmin;
    
end
    


function mon_lbl = get_monlbl(month, doy)

    mons=['JAN';'FEB';'MAR';'APR';'MAY';'JUN';'JUL';'AUG';'SEP';'OCT';'NOV';'DEC'];
    

    lll=' +...+...+...+...+...+...+...+...+...+...+...+...+ ';
    len = length(lll)-2;
    ix1 = 1+max(1,min(len-1,floor(len*doy/366-1)));
    ix2 = ix1+2;
    mon_lbl = lll;
    mon_lbl(ix1:ix2) = mons(month,:);
end

function [fname, varname, lbl, is_us, mode, year_range, frameRate, stepsize, kernelwidth, isPrecip, isTemps, ...
                                        latrange, lonrange, lats, lons, datevecs, cal, fullrange, verbose, minmax, fignum, clim_yrs] = init_params(fname, varname, lbl, varargin)


    varnames = ["tasmax","tmax","tasmin","tmin","tas","tavg","pr","prcp","prec","precip","rh","relhum"];
    tempnames =["tasmax","tmax","tasmin","tmin","tas","tavg"];
    prcpnames =["pr","prcp","prec","precip"];
                    % parse input
    p = inputParser;

    
    addRequired(p,"fname",                  @(s) ischar_s(s));
    addRequired(p,"varname",                @(s) ischar_s(s) && any(strcmpi(s,varnames)));
    addRequired(p,"lbl",                    @(s) ischar_s(s));   
    addParameter(p,"latrange",  [],         @(s) isempty(s)   || (isnumeric(s) && length(s)==2 && s(1)>= -90 && s(2) >= s(1) && s(2) <= 90));
    addParameter(p,"lonrange",  [],         @(s) isempty(s)   || (isnumeric(s) && length(s)==2 && s(1)>= -180 && s(2) >= s(1) && s(2) <=360 && s(2)-s(1) <= 360));
    addParameter(p,"is_us",     [],         @(s) isempty(s)   || islogical(s));
    addParameter(p,"year_range",[],         @(s) isnumeric(s) && s(1) >= 1850 && s(2) >=s(1) && s(2) <=2150);
    addParameter(p,"mode",      0,          @(s) isnumeric(s) && s>= 0 && s <= 5);
    addParameter(p,"pauselen",  [],         @(s) isempty(s)   || (isnumeric(s) && s <= 2));
%   addParameter(p,"matname",   "",         @(s) isempty_s(s) || ischar_s(s));
    addParameter(p,"stepsize",  28,         @(s) isempty(s)   || (isnumeric(s) && s>0 && s <= 365));    
    addParameter(p,"filterwidth",[],        @(s) isempty(s)   || (isnumeric(s) && s>0 && s <= 365));
    addParameter(p,"dir",       [],         @(s) isempty_s(s) || isfolder(s));
    addParameter(p,"verbose",   true,       @(s) islogical(s) || (s==0 || s==1));
    addParameter(p,"fignum",    1,          @(s) isnumeric(s) && s > 0);
    addParameter(p,"minmax",    [],         @(s) isnumeric(s) && (isempty(s) || length(s) == 2));
    addParameter(p,"clim_yrs",  [],         @(s) isnumeric(s) && (isempty(s) || length(s) == 2));
            
    parse(p, fname, varname, lbl, varargin{:});

    latrange = p.Results.latrange;
    lonrange = p.Results.lonrange;
    fname       = p.Results.fname;
    varname     = p.Results.varname;
    lbl         = p.Results.lbl;
    mode        = p.Results.mode;
    year_range  = p.Results.year_range;
    is_us       = p.Results.is_us;
%   matname     = p.Results.matname;
    stepsize    = p.Results.stepsize;
    kernelwidth = p.Results.filterwidth;
    pauselen    = p.Results.pauselen;
    basedir     = p.Results.dir;
    verbose     = p.Results.verbose;
    minmax      = p.Results.minmax;
    fignum      = p.Results.fignum;
    clim_yrs    = p.Results.clim_yrs;
    
    if (isempty(kernelwidth))
        kernelwidth = stepsize;
    end

    if (~exist('pauselen','var'))
        pauselen=[]; 
    elseif (pauselen > 2)
        error("pauslen (%.2f) > 2 second.  To use a pause of longer than 1 second, use minus the pause desired", pauselen);
    else
        pauselen = abs(pauselen);
    end
    c=extractBefore(fname,2);    
    if (~isempty_s(basedir) && c ~= filesep), fname=fullfile(basedir,fname); end
    if (~isfile(fname)), error("error:  file %s does not exist or is not readable", fname); end
    [~,~,ext] = fileparts(fname);
    usemat = strcmpi(ext,".mat");
    if (usemat)
        T = load(fname,   "stepsize", "kernelwidth","lats","lons","cal","datevecs");
        cal = T.cal;
        datevecs = T.datevecs;
        tstamps = datenum_cal(datevecs,cal);
        if (~isempty(stepsize))
            if (T.stepsize ~= stepsize), error("specified stepsize (%f) doesn't match matfile's stepsize (%f).  Please use .nc file instead of matfile", stepsize, T.stepsize); end  
        else
            stepsize = T.stepsize;
        end
        if (~isempty(kernelwidth))
            if (T.kernelwidth ~= kernelwidth), error("specified kernelwidth (%f) doesn't match matfile's kernelwidth (%f).  Please use .nc file instead of matfile", kernelwidth, T.kernelwidtth); end
        else
            kernelwidth = T.kernelwidth;
        end
        lats = T.lats;
        lons = T.lons;
    else
        if (isempty(kernelwidth)), kernelwidth = stepsize; end
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
    
    if (isempty(year_range)), year_range = [datevecs(1,1),datevecs(end,1)]; end
    if (calendar_length(cal) >= 365)
        date_range = [year_range(1),1,1; year_range(2),12,31];
    else
        date_range = [year_range(1),1,1; year_range(2),12,30];
    end
    if (~isempty(clim_yrs))
        clim_yrs = clim_yrs - year_range(1)+1;
    end
    

    tstamp_range = [datenum_cal(date_range(1,:),cal), datenum_cal(date_range(2,:),cal)];
    
    tstamp_range = [max(tstamp_range(1),file_tstamprange(1)), min(tstamp_range(2),file_tstamprange(end))];
    
    if (all(tstamp_range == file_tstamprange) && all(latrange == file_latrange) && all(lonrange == file_lonrange))
        fullrange = true;
    else
        fullrange = false;
    end
    
    if (isempty(pauselen))
        if (stepsize > 30)
            pauselen=.2;
        elseif (stepsize > 15)
            pauselen=.1;
        elseif (stepsize > 7)
            pauselen=.025;
        else
            pauselen=0;
        end    
    end
    frameRate = ceil(30/(1+10*pauselen));
end

function [v, vmax, vmin, datevecs, lats, lons] = my_ncload(fname, varname, latrange, lonrange, year_range, fullrange, vbose)

    my_vbose = ncObj.verbose(vbose);
    
    [~,fn,fext] = fileparts(fname);
    myfname = strcat(fn,fext);
    fprintf("reading netcdf file %s\n", myfname);

    cal = '365-day';
    if (fullrange)
        nc = ncdf_read_file(fname, varname, [], [], [], cal);
    else
        nc = ncdf_read_file(fname, varname, latrange, lonrange, year_range, cal,false, varname, false, -1);
    end
    fprintf("done reading\n");
    v = nc.getvardata(varname);
    vmax = max(v,[], 'all', "omitnan");
    vmin = min(v,[], 'all', "omitnan");
    
    lats = nc.getvardata("lat");
    lons = nc.getvardata("lon");
    
    [~,mycal, ~, ~, ~, tstamps] = ncdf_get_time_info(nc);
    datevecs = datevec_cal(tstamps, mycal);
    
    ncObj.verbose(my_vbose);

end

% function [v, vmax, vmin, datevecs, lats, lons, filename] = my_matload(fname, varname, lats, lons, datevecs, cal, latrange,lonrange,year_range,fullrange)
% 
%     matobj = matfile(fname);
%     
%     vname = matobj.vname;
%     if (~strcmp(vname, varname)), error("error:  varname %s doesn't match varname %s in matfile %s", varname, vname, fname); end
%     vmin = matobj.vmin;
%     vmax = matobj.vmax;
%     isPrecip = matobj.isPrecip;
%     filename = matobj.filename;
%     fprintf("reading %s\n", varname);
%     if (fullrange)
%         v = matobj.(vname);        
%         daterange = [matobj.datevecs(1,1:3); matobj.datevecs(end,1:3)];
%     else
%         daterange = [year_range(1),1,1; year_range(2),12,31-1*(calendar_length(cal)<365)];
%         [timix1, timix2, latix1, latix2, lonix1, lonix2] = range_indexes(lats, lons, datevecs, cal, latrange, lonrange, daterange);    
%         v = matobj.(vname)(timix1:timix2, latix1:latix2, lonix1:lonix2);
%     end
%     fprintf("done reading\n");
%     
%     datevecs = datevec_cal(datenum_cal(daterange(1,:),cal):datenum_cal(daterange(2,:),cal),cal);
% 
%     mynans = v==65535;
%     v = single(v) / single(65534) * single((vmax - vmin)) + single(vmin);
%     v(mynans) = nan;
%     if (isPrecip)
%         vmax = 10^vmax-1;
%         vmin = 10^vmin-1;
%     end
%           
% end
    
% function [timix1, timix2, latix1, latix2, lonix1, lonix2] = range_indexes(lats, lons, datevecs, cal, latrange, lonrange, daterange)
% 
%     tstamps = datenum_cal(datevecs, cal);
%     trange  = datenum_cal(daterange, cal);
%     got_times=false;
%     got_lats = false;
% %   got_lons = false;
%     try
%         timix1 = find(tstamps >= trange(1),1);
%         timix2 = find(tstamps <= trange(2),1,'last');
%         got_times = true;
%         
%         latix1 = find(lats >= latrange(1),1);
%         latix2 = find(lats <= latrange(2),1,'last');
%         got_lats = true;
% 
%         lonix1 = find(lons >= lonrange(1), 1);
%         lonix2 = find(lons <= lonrange(2), 1, 'last');
%     catch
%         if (~got_times)
%             error("no data in specified date range %s - %s.  File:  %s - %s", datestr_cal(datevecs), datestr_cal(tstamps(1)), datestr_cal(tstamps(end)));
%         elseif (~got_lats)
%             error("no data in specified lat range %.4f - %.4f.  File:  %.4f - %.4f", latrange, lats(1), lats(end));
%         else
%             error("no data in specified lon range %.4f - %.4f.  File:  %.4f = %.4f", lonrange, lons(1), lons(end));
%         end
%     end
% end
    
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

% function savemat(matname, varname, v, lats, lons, datevecs, cal, stepsize, kernelwidth, filename, isPrecip)
% 
% %     if (isPrecip)
% %         v = log10(1+v);       % already done in main...
% %     end
%     vmin = min(v(:),[],"omitnan");
%     vmax = max(v(:),[],"omitnan");            
% % 
%     vname = varname; 
%     fprintf("Saving averages as %s.  See README in matfile to recover original values\n", matname)
%     ui_v = uint16(round((v-vmin)/(vmax-vmin)*65534));
%     ui_v(isnan(v)) = 65535; %#ok<NASGU>
%     eval (sprintf("%s = ui_v;", varname));
%     if (isPrecip)
%         README = sprintf("data is stored as unsigned shorts (with NAN's stored as 65535);  to recover (approx) original values: vals = double(%s)/65534.0*(vmax-vmin)+vmin); prcp = 10.^(vals) - 1; prcp(%s==65535)=nan;  ", varname, varname); %#ok<NASGU>
%     else
%         README = sprintf("data is stored as unsigned shorts (with NAN's stored as 65535);  to recover (approx) original values: vals = double(%s)/65534.0*(vmax-vmin)+vmin); vals(%s==65535)=nan;", varname, varname); %#ok<NASGU>
%     end
% 
%     save(matname, "vmin", "vmax", "cal", "stepsize", "kernelwidth", "filename", "isPrecip", "lats", "lons", "datevecs", "vname", varname);
% 
% end

function [vout, vmin, vmax] = remove_climatology(v, clim_yrs, mode)

    vclim = calc_climatology(v, clim_yrs, mode);
    
    ndays = size(v,1);
    nyrs = ndays/365;
    if (mode == 3)
        vout = v - repmat(vclim,nyrs,1,1);
    else
        vout = v - vclim;
    end

    vmin = min(vout(:),[],"omitnan");
    vmax = max(vout(:),[],"omitnan");
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

    vclim = real(ifft(fft(squeeze(mean(v, 2,"omitnan"))) .* repmat(FILT,1,nlats, nlons)));
    
end

