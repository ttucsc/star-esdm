function tmax_stats_maps(ncname, varargin)
% function tmax_stats(ncname, figname, downsample, lats, lons, yr_range, yr_step, thresh, consec_cnt, varname, fignum)
%
%   program to produce maps of temperatures above a given threshold (default 100 deg F)
%       Saves output to png files.
%
%   written as a quickie to get some maps for John Zak.
%   Defaults:
%       ncname          (no default).   name of nc file to read.
%       figname         fig_%03d.png    template for figure names: 
%       downsample      4               downsampling to smaller size.     Will average over eac 4x4 region.  For data downscaled to
%                                                                         livneh gridded obs, this takes it from 1/16 to 1/4 degree.
%       lats, lons      all lats, lons in file
%       yr_range        all years in file
%       yr_step         [10, 10];       (1 or 2 values) 
%                                       if one value, # of years to average and # of years between steps.
%                                       if two values, 1st is # of years to average, 2nd is # of years between steps.
%       thresh          100 F (37.78 C) temperature threshold:  
%       consec_cnt      1               map # of times consecutive days are above threshold  (Other values than 1: not tested!)
%       varname         (from file)     [will look for 1st standard variable in file]
%       fignum          1:15            if figname ends in .mj2 default is 1
%
%   Uses 9 worker threads to create the maps.  Doesn't write any maps until all mapping is done.
%
%       Ian:    for making a movie, see:  https://stackoverflow.com/questions/17409752/matlab-create-movie
%               also, Answer 1 shows how to export movie frame by frame, rather than create from stored movie.
%       Todo:       1.  Add histogram of counts on right
%                   2.  Perhaps low-pass filter, and matlab contour(...) instead of colormap trick. 


    [ncname, figname, downsample, lats, lons, yr_range, yr_step, thresh, consec_cnt, varname, units, fignum, do_movie, is_conus, lbl] = init_params(ncname, varargin{:});
    
    nlats = length(lats);
    nlons = length(lons);
    end_yr = yr_range(2)-yr_step(2)+1;
    nyrs  = end_yr - yr_range(1)+1;
    nsteps = ceil(nyrs/yr_step(1));
    st_yrs = yr_range(1):yr_step(1):end_yr;
    
    mlats = nlats/downsample;
    mlons = nlons/downsample;
    
    latrange = [lats(1),lats(end)];
%   lonrange = [lons(1),lons(end)];
    maps = zeros(mlats,mlons, nsteps, "uint16");
    
%   ARRM_V2_start_workers(14);
    
    logname = "tmax_stats_map.log";
    fprintf("logname: %s\n", logname)
    fidlog = fopen(logname,"w");
    fprintf("%s logfile, %d loops\n", mfilename, mlons);
    fclose(fidlog);
    ncdf.verbose(false);
    
    parfor lonix = 1:mlons
        if (downsample == 1)
            lon = lons(lonix);
            nc = ncdf_read_file(ncname, varname, latrange, lon, yr_range, '365-day');
            d = squeeze(nc.getvardata(varname));
        else            
            ix1 = (lonix-1)*downsample+1;
            ix2 = lonix*downsample;
            lon = lons([ix1,ix2]); %#ok<PFBNS>
            nc = ncdf_read_file(ncname, varname, latrange, lon, yr_range, '365-day');
            d = nc.getvardata(varname);
%           d = squeeze(nanmax(d,[],3));
            d = squeeze(max(d,[],3,"omitnan"));
            ndays = size(d,1);
            if (downsample > 1)
                d = reshape(d,ndays,downsample,mlats);
%               d = squeeze(nanmax(d,[],2));
                d = squeeze(max(d,[],2));
            end   
        end
        if (strcmp(units,"K"))
            d = d - 273.15;
        elseif (any(strcmp(units,["F","degF","deg F"])))
            d = (d-32)*5/9;
        end
        d = fillmissing(d,"nearest");   % replace nan's w/ nearest non-nan value.
        
        counts = calc_counts(d, st_yrs, yr_step, thresh, consec_cnt, do_movie); 
        maps(:,lonix,:) = counts'; 
        fidlog = fopen(logname,"a");
        fprintf(fidlog, "%6d of %6d\n",lonix, mlons);
        fclose(fidlog);
        [~,cnt] = system(sprintf("cat %s | wc -l  ", logname));
        cnt = str2double(cnt) - 1;
        pct = cnt/mlons * 100;
        if (mod(cnt, ceil(mlons/20))==0), fprintf("%4.f %% completed\n", pct); end
%         show_progress(cnt, nlons); 
    end
    fprintf("\n");
            
    save_figs( maps, st_yrs,  yr_step, fignum,  figname, lats, lons, thresh, do_movie, is_conus, lbl);

end

function cnts = calc_counts(data, st_yrs, yr_step, thresh, consec_cnt, do_movie)
    
    st_yrs = st_yrs - st_yrs(1)+1;      % make st_yrs relative, so it starts at 1 instead of 1951.
    [ndays,nlats] = size(data);
    nyrs = ndays/365;
    data = single(data >= thresh);
    nsteps = length(st_yrs);
    if (consec_cnt > 1)
        n = ndays + consec_cnt;
        krnl = make_krnl(consec_cnt,n);
        data = ifft(fft(data,n,dim) .* fft(krnl),n,dim);
        data = data(1:ndays,:);
        data = mod(data+1,consec_cnt)==0;

    end
    data = reshape(data,365,nyrs,nlats);
    data = squeeze(sum(data,1));
    
    cnts = zeros(nsteps,nlats,"uint16");
    if (do_movie)
        data = [repmat(data(1,:),2*yr_step(2),1); ...   % extend data so we don't have overlap with fft
              data; ...
              repmat(data(end,:),2*yr_step(2),1)];
        myrs = nyrs + 4*yr_step(2);
        
        sig = sqrt(yr_step(2)^2-1)/sqrt(12);
        G = GAUSS_fft(myrs, sig, true)';
        data = abs(ifft(fft(data) .* G));
        istart = 2*yr_step(2)+1;
        iend   = istart + nyrs-1;
        data = data(istart:iend,:);
        for i=1:nsteps
            ix = st_yrs(i) + floor(yr_step(2)/2);
            cnts(i,:) = data(ix,:);
        end
    else

        for i=1:nsteps
            ix1 = st_yrs(i);
            ix2 = ix1 + yr_step(2)-1;
            try
                cnts(i,:) =  round(sum(data(ix1:ix2,:),1)/yr_step(2));
            catch
                oops();
            end
        end
    end
%     data=reshape(data,yr_step,nsteps,nlats);
%     cnts = squeeze(sum(data));
%     cnts = cnts / myrs;
end

function krnl = make_krnl(len, nyrs)

    krnl = zeros(nyrs,1);
    krnl(1:len)=1;
    krnl = circshift(krnl,-floor(len/2),1);
end


function save_figs(maps, st_yrs,  yr_step, fignum,  figname, lats, lons, thresh, do_movie, is_conus, lbl)

    [nlats,nlons,nmaps] = size(maps);
    
%   save("test.mat");
    
    
    coast_flag = 2;
%     cm1 = flipud(hot(100));
%     cm2 = flipud(hsv(300));
%     cmap = [cm1(1:50,:);cm2(1:150,:);zeros(1,3)];
    cmap = temp_gradient_1(200,20,[.94,.97,1],[],[], []);
    ming = 0;
    maxg = 200; %max(maps(:))
%   maxg2 = maxg; %max(maps(:))

    if (do_movie)
        [~,~,fext] = fileparts(figname);
        if (strcmpi(fext,".mj2"))
            v = VideoWriter(figname, "Motion JPEG 2000");
            v.CompressionRatio = 10;
        elseif (strcmpi(fext,".avi"))
            v = VideoWriter(figname, "Motion JPEG AVI");
            v.Quality = 33;
        elseif (strcmpi(fext,".mp4"))
            v = VideoWriter(figname, "MPEG-4");
            v.Quality = 33;
        end
%       v.Colormap=cmap;
        v.FrameRate=8;
        open(v);
    end
    
    
    
    dlat = lats(10)-lats(9);
    dlon = lons(10)-lons(9);
    R = make_refmat(lats, lons, nlats, nlons, dlat, dlon);
    substeps = do_movie;
    for i=1:nmaps
        for jj=0:1*substeps
            if (i == nmaps || jj==0) 
                j=i; 
            else
                j=i+1; 
            end
            figix = min(length(fignum), i);
            fign = fignum(figix);
            figure(fign);
            fthresh = round(thresh*9/5 + 32);
            ctr_yr = st_yrs(i)+floor(yr_step(2)/2);
%           decade = sprintf("%4d's",ctr_yr-mod(ctr_yr,10));
            if (do_movie)
                ttl = sprintf("Days per year above %d F, 10-yr avg, %d, %s", fthresh, ctr_yr, lbl);
            else            
                ttl = sprintf("Days per year above %d F, %d - %d, %s", fthresh, st_yrs(i), st_yrs(i)+yr_step(2)-1, lbl);
            end
            if (jj==0)
                mymap = squeeze(maps(:,:,i));
            else
                w1 = (1-jj/2);
                w2 = 1-w1;
                mymap = squeeze(w1*maps(:,:,i) + w2*maps(:,:,j));   % add an in-between frame so movie is smoother
            end
            mymap(mymap>maxg) = maxg;
    %       lvls = [25,100];
    %         M = contourc(double(mymap),lvls);
    %         mymap = draw_contours(mymap, M, lvls, [200,201]);

            h = display_map(  mymap, ttl,  R, fign, [],       ming,  maxg, coast_flag, cmap, [], [], [100, 600, 1500,1000], 24, false);
            hold on;
    %       [~,hh(1)] = contourm(double(mymap), R, [ 25, 25], '-',"linewidth",2,'color',[.25,.65,.25]);
            [~,hh   ] = contourm(double(mymap), R, [ 0,0],    '-',"linewidth",1,'color',[.98,.98,.98]);
            [~,hh(2)] = contourm(double(mymap), R, [ 30, 30], '-',"linewidth",1,'color',[0.25,.75,1]);
            hh(1).DisplayName = "avg <1 day/yr";
            hh(2).DisplayName = ">  30 days/yr";
            if (sum(mymap>100,'all') > 0)
                [~,hh(3)] = contourm(double(mymap), R, [90,90], '-',"linewidth",2,'color',[.80,.80,.80]);
                hh(3).DisplayName = ">  90 days/yr";
%                 if (sum(mymap>180,'all')> 0)
%                     [~,hh(4)] = contourm(double(mymap), R, [180,180], '-',"linewidth",1,'color',[.15, .5, .15]);
%                     hh(4).DisplayName = "> 180 days/yr";
%                 end
            end
            legend(hh,'FontSize',16);
                % write decade in Gulf of Mexico, and move it over slowly year by year.
            textstep = 3.5/nmaps;
            textlon = -96 + textstep*i;
            fsize = 64 - 12*is_conus;
            textm(27.5,textlon, sprintf("%d", ctr_yr),"fontsize", fsize, "color",[.75,.75,.75])
            hold off;
            drawnow();
            pause(.1);
            if (do_movie)
                f = getframe(h);
                writeVideo(v,f);           
                if (jj==0), show_progress(i, nmaps); end
            else
                outname = sprintf(figname,st_yrs(i));
                saveas(h, outname);
            end
        end        
    end
    
    if (do_movie)
        close(v);
    end
end
    
function     [ncname, figname, downsample, lats, lons, yr_range, yr_step, thresh, consec_cnt, varname, units, fignum, do_movie, is_conus, lbl] = init_params(ncname, varargin)
    
    p = inputParser;
    addOptional(p, "figname",      "fig_%03.png", @(s) ischar_s(s));
    addOptional(p, "downsample",    4, @(s) mod(s,1)==0 && s > 0 && s < 32);
    addOptional(p, "latrange",     [], @(s) isnumeric(s) && any(length(s)==[0,2]));
    addOptional(p, "lonrange",     [], @(s) isnumeric(s) && any(length(s)==[0,2]));
    addOptional(p, "yr_range",     [], @(s) isnumeric(s) && any(length(s)==[0,2]));
    addOptional(p, "yr_step", [10,10], @(s) all(s >= 1 & s <= 50))
    addOptional(p, "thresh",  37.7778, @(s) s > -30 && s < 150) ;  % 37.7778 C = 100 F;
    addOptional(p, "consec_cnt",    1, @(s) s >= 1 && s <= 50);
    addOptional(p, "varname",      "", @(s) ischar_s(s));
    addOptional(p, "fignum",        []);
    addOptional(p, "lbl",           "");

    parse(p, varargin{:});
    
    figname     = p.Results.figname;
    downsample  = p.Results.downsample;
    latrange    = p.Results.latrange;
    lonrange    = p.Results.lonrange;
    yr_range    = p.Results.yr_range;
    yr_step     = p.Results.yr_step;
    thresh      = p.Results.thresh;
    consec_cnt  = p.Results.consec_cnt;
    varname     = p.Results.varname;
    fignum      = p.Results.fignum;
    lbl         = p.Results.lbl;
    
    [ok, nc] = isnetcdf(ncname,true);  
    if (~ok), error("error:  %s is not a netcdf file", ncname); end
    
    if (isempty_s(varname)) 
        varname = ncfind_varinfo_ic(ncname); 
    end
    units = string(nc.get(sprintf("/Variables/%s/Attributes/units/Value",varname)));
    
    % get label from filename if not provided
    
    if (isempty_s(lbl))
        if (~contains(ncname,"downscaled"))
            lbl = "Livneh Gridded Obs";
        else
            run_info = ARRM_V2_parse_netcdf_filenames(ncname); 
            ss = extractAfter(run_info.scenario, "rcp");
            scen = str2double(ss)/10;
            lbl = sprintf("%s RCP %.1f" , run_info.model, scen);
        end
        
    end
    % get thresh in deg C
    
    if (thresh > 170)
        thresh = thresh - 273.15;
    elseif (thresh > 65)
        thresh = (thresh-32)*5/9;
    end
    
    % get lats,lons and dates in range
    
    [lats, lons] = ncdf_get_latlons(nc);
    
    if (~isempty(latrange)), lats = lats(lats >= latrange(1) & lats <= latrange(2)); end
    nlats = length(lats);
    mlats = floor(nlats/downsample);
    lats = lats(1:(mlats*downsample));
    
    if (any(lonrange < 0) && ~any(lons < 0)), lonrange = lonrange + 360; end
    is_conus = true;
    if (~isempty(lonrange))
        lons = lons(lons >= lonrange(1) & lons <= lonrange(2));
        is_conus = false;
    end
    nlons = length(lons);
    mlons = floor(nlons/downsample);
    lons = lons(1:(mlons*downsample));
    
    [~,cal,~,~,~,tstamps] = ncdf_get_time_info(nc);
    
    dvecs = datevec_cal(tstamps, cal);
    yrs = unique(dvecs(:,1));
    
    if (~isempty(yr_range))
        yr_range = [max(yr_range(1), yrs(1)), min(yr_range(2), yrs(end))];
    else
        yr_range = [yrs(1),yrs(end)];        
    end
    
        % truncate data to exact multiple of yr_steps
    if (length(yr_step) == 1), yr_step = [yr_step, yr_step]; end    
    nyrs = yr_range(2) -yr_step(2) - yr_range(1)+1;
    nsets = floor(nyrs/yr_step(1));
    yr_range(2) = yr_range(1) + nsets*yr_step(1) + yr_step(2)-1;
    
    [~,~,fext] = fileparts(figname);
    do_movie = any(strcmpi(fext,[".mj2",".avi",".mp4"]));
    if (isempty(fignum))
        if (do_movie)
            fignum = 1;
        else
            fignum = 1:nsets; 
        end
    end
    
end
    

    


    
    
    
        