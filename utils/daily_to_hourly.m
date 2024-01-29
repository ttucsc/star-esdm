% function hourly = daily_to_hourly(hourly, minvals, maxvals, isRH, finalMinutes, inputMinutes, allResults)
function [hourly, dailyClim, obsmean, mmPts, mmVals, mmPtsobs, mmValsobs, mmPtsobsraw, mmValsobsraw, obsHourly, nptsPerDay] = daily_to_hourly(hourly, minvals, maxvals, isRH, finalMinutes, inputMinutes)
% function to expand tmin/tmax or rhmin/rhmax data to hourly, by creating daily templates from obsHourly.
%

%   hourly          hourly observations.  Full years, 365-day calendar, starting 1/1, ending 12/31
%   minvals,maxvals daily min/max to be expanded to hourly, full years, 365 day calendar, 1/1 - 12/13
%   isRH            boolean.  true for Relative Humidity.  False for temperature data.
%   finalMinutes    (optional)  # of minutes between samples in output.  [60].
%   inputMinutes    (optional)  # of minutes between samlpes in input.   [60].
%
%   Assumption is that tmin and tmax are sampled on a local 24-hour clock.
%   For temperature data, tmin should occur in the early morning and tmax occurs mid-afternoon.  
%   For relative humidity, rhmax should occur in the early morning (usually around time of tmin) 
%   and rhmax in the afternoon (usually around time of tmax).
%
%   Program calculates the mean hourly values over the obsHourly data, smooths to determine the daily signal over the
%   365-day year, and applies this daily shape to recreate an hourly signal from minvals and maxvals.
%
%   finalMinutes can be provided to get results on a finer or courser grain than hourly.
%   If finalMinutes is smaller than the sampling period of the input data (assumed hourly)

%   Usually you just want to use
%       hourly = daily_to_hourly(...)
%   The additional return parameters are for producing plots of the output, which is done by the program obs_hourly_surface2
%   It will run a little faster if you don't grab all the return variables.

    if (nargout < 7)
        mmValsobs = [];
        mmPtsobsraw = [];
        mmValsobsraw = [];
    end

        % run params for creating hourly table
    nterms_day = 6;
    sigterm_day = 1;
    nterms_yr = 6;
    sigterm_yr = 3;
    minutes = 15;       % resample obsHourly to 15-minute steps.  This better locates the position of daily min & max.
                        % also allows min & max to exceed observed values (via fourier resampling).    
    if (~exist('finalMinutes','var') || isempty(finalMinutes))
        finalMinutes = 60;      % expand to hourly.
    end
    
    if (~exist('inputMinutes','var') || isempty(inputMinutes))
        inputMinutes = 60;      % expand to hourly.
    end
    minutes = min(inputMinutes, minutes);
    inptsPerDay = 24*inputMinutes/60;    
    inptsHourly = numel(hourly);
    inyrsHourly = inptsHourly/inptsPerDay/365;
    inptsDaily  = numel(minvals);
%   inyrsDaily  = inptsDaily/365;
    if (mod(inptsHourly,365) ~= 0 || mod(inptsDaily,365) ~= 0 || numel(maxvals) ~= inptsDaily)
        error('hourly_to_daily:  error:  input data sizes not compatible.  Must be 365-day calendar, full year.');
    end

     hourly  = fixnans(hourly);
     minvals = fixnans(minvals);
     maxvals = fixnans(maxvals);
     
    resampleRate = inputMinutes/minutes;
    nptsPerDay = resampleRate*24;
    if (resampleRate ~= 1)
        outlen = inptsHourly*resampleRate;
        obsHourly = fft_resample(hourly(:),  outlen);
        if (isRH)
            obsHourly = min(100,max(0,obsHourly));
        end
    else
        obsHourly = hourly;
    end

    if (nargout > 7)
        [mmPtsobsraw, mmValsobsraw] = find_minmax_points(obsHourly, nptsPerDay);      % for displaying results
    end

    obsHourly = reshape(obsHourly, nptsPerDay, 365, inyrsHourly);

    obsmean = mean(obsHourly, 3);

    if (nargout > 6)
        [mmPtsobs, mmValsobs] = find_minmax_points(obsmean, nptsPerDay);              % for displaying results
    end
    
    dailyClim = smooth_2d(obsmean, nterms_day, sigterm_day, nterms_yr, sigterm_yr);
    [mmPts, mmVals, mmWts] = calc_weights_and_points(dailyClim, isRH, nptsPerDay);  % for displaying results
%   [mmPts, ~,      mmWts] = calc_weights_and_points(dailyClim, isRH, nptsPerDay);

    reconstructed = reconstruct_daily(minvals, maxvals, mmWts, mmPts, isRH);

    if (minutes ~= finalMinutes)
        outlen = numel(reconstructed)*minutes/finalMinutes;
        hourly = fft_resample(reconstructed, outlen);
    else
        hourly = reconstructed;
    end

end

function  signal  = reconstruct_daily(minvals, maxvals, mmWts,  mmPts, isRH)
% let n = nptsPerDay = size(mmWts,1);
% tmin, tmax    array of daily tmins & tmaxs
% mmWts         tmin & tmax weights for each point in year.             Size: n x 365 x 2  1=tminWt, 2=tmaxWt
% mmPts         index of point for each day when tmin and tmax occur.   Size:     365 x 2  1=tminPt, 2=tmaxPt
    ndays = length(minvals);
    nptsPerDay = size(mmWts,1);
    
    totpts = ndays * nptsPerDay;
    signal = nan(totpts,1);
    
    minPts = mmPts(:,1);
    maxPts = mmPts(:,2);
    minWts = mmWts(:,:,1);
    maxWts = mmWts(:,:,2);
    
    minvals(end+1) = minvals(end);
    maxvals(end+1) = maxvals(end);
    minNow = minvals(1);               % mmWts scale current minval and maxval to get current signal level
    maxNow = maxvals(1); 

    for iday=1:ndays

        idoYr = mod(iday-1,365)+1;
        inext = iday+1;

        if (~isRH)
            for itime=1:nptsPerDay
                if (itime == minPts(idoYr))

                    maxNow = maxvals(iday);
                end
                if (itime == maxPts(idoYr))

                    minNow = minvals(inext);
                end
                ix = (iday-1)*nptsPerDay + itime;

                signal(ix) = minWts(itime, idoYr)*minNow + maxWts(itime,idoYr)*maxNow;

            end
        else
            for itime=1:nptsPerDay
                if (itime == minPts(idoYr))

                    maxNow = maxvals(inext);
                end
                if (itime == maxPts(idoYr))

                    minNow = minvals(iday);
                end
                ix = (iday-1)*nptsPerDay + itime;
    
                signal(ix) = minWts(itime, idoYr)*minNow + maxWts(itime,idoYr)*maxNow;

            end
        end
    end
end

function signal = smooth_2d(data, nterms_y, sigterm_y, nterms_x, sigterm_x)
%   This isn't a true 2-D smooth...it treats the data as one long time series, and first smooths along the entire
%   series (assumes continuity along the y-dimension.  Then it reshapes it back to the origal shape, and smooths each
%   row individually along the x direction.

    [ny,nx] = size(data);
            % filter over the entire time period.  
            % nterms_y and sigterm_y are for 1 day (ny pts), so calc_filter creates a filter for the entire period
            % and we filter it as one sequence.  This is valid for our time series data, since each day leads into the
            % next.
    FILT = calc_filter(nterms_y, sigterm_y, nx, ny);
    
    signal  = lpf_FILT(data(:),  FILT);
    
    signal  = reshape(signal, ny, nx);
    
        % now smooth along the day-or-year axis by
        % calculating climatology for each hour of the day
    signal  = horizontal_smooth(signal,  nterms_x, sigterm_x); 
    
end   

function clim = horizontal_smooth(data, nterms, sig_term)
% data is n x 365 

    [ny, nx] = size(data);
    
    clim=zeros(size(data));
    
    FILT = calc_filter(nterms, sig_term, 1, nx);
    for i=1:ny
        clim(i,:) = lpf_FILT(data(i,:), FILT);
    end    
end

function [mmPts, mmVals, mmWts] = calc_weights_and_points(signal, isRH, nptsPerDay)
        % interpolate to 5-minute steps so we can identify Tmax and Tmin points
        
    [mmPts, mmVals]  = find_minmax_points(signal, nptsPerDay);
    mmWts  = calc_mmWts(signal,  mmPts(:,1),  mmPts(:,2), mmVals(:,1), mmVals(:,2), isRH);
    
end

function mmWts = calc_mmWts(signal, minPts, maxPts, minVals, maxVals, isRH)

    [nptsPerDay,ndays] = size(signal);
    
    mmWts = nan(nptsPerDay, ndays, 2);
    
    minVals(end+1) = minVals(end); % assume day after end's tmin is same as last day's tmin.     
    maxVals(end+1) = maxVals(end); % assume day after end's tmin is same as last day's tmin.     
    minNow = minVals(1);      
    maxNow = maxVals(1);      % assume yesterday's max is same as today's.
    for iday=1:ndays
        inext = iday+1;
        if (iday == 400 || iday == 406)
            fprintf('here %d\n', iday);
        end
        if (~isRH)
            for itime=1:nptsPerDay
                if (itime == minPts(iday)), maxNow = maxVals(iday); end
                if (itime == maxPts(iday)), minNow = minVals(inext); end
                mmWts(itime, iday, 1) = (maxNow - signal(itime, iday))/(maxNow - minNow);
                mmWts(itime, iday, 2) = 1- mmWts(itime, iday, 1);
            end
        else
            for itime=1:nptsPerDay
                if (itime == minPts(iday)), maxNow = maxVals(inext); end
                if (itime == maxPts(iday)), minNow = minVals(iday); end
                mmWts(itime, iday, 1) = (maxNow - signal(itime, iday))/(maxNow - minNow);
                mmWts(itime, iday, 2) = 1- mmWts(itime, iday, 1);
            end
        end
    end


end

function [minVals, maxVals] = find_minmax(data, nptsPerDay)

    ndays = numel(data)/nptsPerDay;
    data = reshape(data, nptsPerDay, ndays);
    
    minVals = min(data);
    maxVals = max(data);
    
end


function [mmPts, mmVals]  = find_minmax_points(signal, nptsPerDay)

    ndays = numel(signal)/nptsPerDay;
    signal = reshape(signal, nptsPerDay, ndays);
    
    mmPts  = nan(ndays, 2);
    mmVals = nan(ndays, 2);
    
    [mins, maxs] = find_minmax(signal, nptsPerDay);
    
    for i=1:ndays
        mmPts(i,1) = find(signal(:,i) == mins(i), 1);
        mmPts(i,2) = find(signal(:,i) == maxs(i), 1);
    end
    
    mmVals(:,1) = mins;
    mmVals(:,2) = maxs;
end


function [obs, missing] = fixnans(obs, maxlen, use_these)

    if (~exist('maxlen','var') || isempty(maxlen)), maxlen = 24; end
    if (~exist('use_these','var')), use_these = []; end
    
    orig_shape = size(obs);
    if (isempty(use_these))
        if (isnan(obs(1)))   % special cases:  nans at start or end.
            ix = find(~isnan(obs),1);
            if (ix > maxlen), error("error:  too many missing values at start."); end
        end
        if (isnan(obs(end)))
            ix = find(~isnan(obs),1,'last');
            if (length(obs)-ix > maxlen), error("error:  too many missing values at end."); end
        end
        [obs, missing] = fillmissing(obs, 'linear','EndValues','nearest');
    else
        
%       [missing, longest, seq] = find_nans(obs);
        [missing, longest, ~  ] = find_nans(obs);

        if (longest > maxlen)
            obs(missing) = use_these(missing);
        else
            [obs, missing] = fillmissing(obs, 'linear', 'EndValues','closest');
        end
    end
    obs = reshape(obs, orig_shape);
end           
    
function [missing, longest, sequences] = find_nans(obs)

    missing = isnan(obs);
    
    sequences=1.0*missing(:);
    
    m1 = find(missing);
    
    if (sequences(1) == 1)
        istart=2;
    else
        istart=1;
    end
    for i=istart:length(m1)
        ix = m1(i);
        sequences(ix) = sequences(ix-1)+1;
    end
    longest = max(sequences);
end
    

