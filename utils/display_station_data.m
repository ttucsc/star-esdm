function display_station_data(stns, ncName, yrRange, stns_per_plot, doSkipFig, savePDFs, saveFigs, figdir, logname)
%   Program to display station data from a station netcdf file
%
%   Inputs:
%       stns            stations to be displayed.  cell array or array of strings
%       ncName          name of station netcdf file
%                           note:  netcdf file has Variable Name (Tmax, Tmin, Prec) in UserData.
%       yr_range        year-range to display [start_year, end_year]  (optional)
%       stns_per_plot   number of stations per figure.  (optional)
%       doSkipFig       boolean.  If true, figure is not displayed.
%                           if false or empty or missing, figure is displayed.
%       savePDFs        boolean.  If true, saves PDFs as well.  
%                           NOTE:  this takes several seconds per figure.
%       saveFigs        boolean.  If true, saves matlab figs
%       figdir          directory where to save figure & pdf.  
%                           if empty or missing, save to '.'

    if (~exist('yrRange','var') || isempty(yrRange))
        yrRange = [1950,2017];
        noDates = true;
    else
        noDates = false;
    end
    if (~exist('doSkipFig','var') || isempty(doSkipFig))
        doSkipFig = true;
    end
    
    if (~exist('stns_per_plot','var') || isempty(stns_per_plot))
        stns_per_plot=9;
    end
    if (~exist('savePDFs','var') || isempty(savePDFs))
        savePDFs = false;
    end
    if (~exist('saveFigs','var') || isempty(saveFigs))
        saveFigs = false;
    end
    if (~exist('figdir','var') || isempty(figdir))
        figdir = '.';
    end
    if (~exist('logname','var'))
        logname=[];
    end
    
    if (ischar(stns)), stns = {stns}; end      % in case there's only 1 station and user didn't give it to us as a cell array or a string.
    
    
    stnTbl = QC_get_site_table(ncName);
    vname = stnTbl.Properties.UserData.varName;

    h = [];
    firstfig = true;
    figcount = 0;
    
    nsites=length(stns);
    if (nsites==1)          % see if user wants a report on the failing sites only
        if (strcmp(stns{1},"fails"))
            stns = make_stn_fail_list(stnTbl, yrRange(1),yrRange(2));
            nsites=length(stns);
        else
            [stns,nsites]=QC_get_site_list_from_file(stns);     % either reads file, or returns stns if not a filename.
        end
    end
    
    if (~isempty(logname))
        if (exist(logname,'file')==2), delete(logname); end
        diary(logname); 
    end
    
    fprintf('stnID       st_yr endyr 50%%yr    %%   75%%yr     %%    Jan  Feb  Mar  Apr  May  Jun  Jul  Aug  Sep  Oct  Nov  Dec ord valid  Site\n');
        
    stns_per_plot   = max(1,min(stns_per_plot, length(stns)));
    fignames        = cell(ceil(nsites/stns_per_plot),1);
    pdfnames        = cell(ceil(nsites/stns_per_plot),1);
    hh              = cell(ceil(nsites/stns_per_plot),1);
    
    for i=1:nsites
        try
            if (iscell(stns))
                stn = stns{i}; 
            else
                stn = stns(i); 
            end
            stnInfo = QC_get_site_table(stnTbl, "stnID", stn, "removeLeaps", true, "searchType","any");
            stnInfo = stnInfo(1,:);
%             dvecs = datevec(stnInfo.Properties.UserData.dates);
            site_id=stnInfo.stnID(1);
            site=stnInfo.stnName(1);
            sdate=datevec(stnInfo.startDate(1));
            edate=datevec(stnInfo.endDate(1));
            if (noDates)
                yr1 = sdate(1);
                yr2 = edate(1);
            else
                yr1 = yrRange(1);
                yr2 = yrRange(2);
            end
            syr = max(sdate(1), yrRange(1));
            eyr = min(edate(1), yrRange(2));
%            stnInfo = QC_get_site_table(stnTbl, yr1, yr2, "stnID", site_id, "removeLeaps",true,"searchType","stnID");
            if (calendar_length(stnInfo.Properties.UserData.calendar) ~= 365), remove_leaps = true; else, remove_leaps = false; end
            stnInfo = QC_get_data(stnInfo, [], yr1, yr2, remove_leaps);     % might be able to replace [] with nc if we read it somehwere earlier, Ian.
            yr_50=[];
            yr_75=[];
            temps=stnInfo.data(1,:);
            npts=length(temps);
            if (mod(npts,365)~=0)
                mpts = ceil(npts/365)*365;
                temps(end+1:mpts)=nan;
            end
            npts=length(temps);
            dv=reshape(temps,365,npts/365);
            y_end=yr2;
            if (y_end == edate(1) && edate(2) < 12)
                y_end=y_end-1;          % end on last 'complete' year.
            end
%             nyrs=yr2-yr1+1;
            y_start=yr1;
            iyr=1;
            
            [~, valids, ord_txt, validFlag] = check_valid(temps, yr1, yr2);
            
                % find 1st year where we have at least 50% or 75% valid points from there to end of data.
            while(isempty(yr_75) && y_start <= yr2)
                ix1=1+365*(y_start-yr1);
                ix2=365*(1+y_end-yr1);
                nnans=sum(isnan(temps(ix1:ix2)));
                frac=1-(nnans/(1+ix2-ix1));
                if (isempty(yr_50) && frac >= .5 && sum(dv(:,iyr))>0)
                    yr_50 = y_start;
                    frac_50=100*frac;
                    yr_50_txt = sprintf('%04d',yr_50);
                elseif (isempty(yr_75) && frac >= .75 && sum(dv(:,iyr))>0)
                    yr_75=y_start;
                    frac_75=100*frac;
                    yr_75_txt = sprintf('%04d',yr_75);
                else
                    y_start=y_start+1;
                    iyr = iyr+1;
                end
            end
            if (isempty(yr_50))
                yr_50 = [];
                frac_50 = 0;
                yr_50_txt = ' -  ';
            end
            if (isempty(yr_75))
                yr_75 = [];
                frac_75 = 0;
                yr_75_txt = ' -  ';
            end                

            fprintf('%-12s %4d  %4d  %4s %5.1f%% %5s %5.1f%%  ', site_id, syr, eyr, yr_50_txt, frac_50, yr_75_txt, frac_75);
                    % output the monthly valid counts
            fprintf('%4d ', valids);
            fprintf(' %1s %6s  %s \n', ord_txt, validFlag, site);

            if (~doSkipFig || saveFigs || savePDFs)
                if (firstfig)
                    if (isempty(h) || ~(saveFigs || savePDFs))
                        if (doSkipFig)
                            h=figure('visible','off');
                        else
                            h=figure();
                        end
                    else
                        clf();
                    end
                    hh{figcount+1} = h;
                end
                [fn,pdfn, figcount, firstfig] = plot_temps(h, figcount, firstfig, i, nsites, temps, yr1, yr2, yr_50, yr_75, y_end, stnInfo, stns_per_plot, vname, doSkipFig, savePDFs, saveFigs, figdir);
                if (~isempty(fn))
                    fignames{figcount} = fn;
                end
                if (~isempty(pdfn))
                    pdfnames{figcount} = pdfn;
                end
            end
                
        catch me
            fprintf('--------%s: caught exception----------\n', mfilename);
            msgtext=getReport(me);
            fprintf('%s\n', msgtext);
            fprintf('------\n');
        end
    end
    
    if (~isempty(logname))
        diary off;
        fprintf('\nresults saved to %s\n\n',logname);
    end
    
    fprintf('\n');
    for i=1:length(fignames)
        if (~isempty(fignames{i}))
            fprintf('fig written: %s\n', fignames{i});
        end
    end
    
    fprintf('\n');
    for i=1:length(pdfnames)
        if (~isempty(pdfnames{i}))
            fprintf('pdf written: %s\n', pdfnames{i});
        end
    end
    fprintf('\n');
    
            % redisplay the figures
    if (length(fignames) < 15)
        if (~doSkipFig)
            for i=1:figcount
                if (~isempty(fignames{i}))
                    if (ishandle(hh{i})), close(hh{i}); end
                    h=openfig(fignames{i});
                else
                    h=hh{i};
                end
                set(h,'Units','Normalized');
                pos = [.01+.01*(i),.95-.005*i,.8,.8];
                set(h,'Position',pos);
            end
        end       
    end
    if (~usejava('desktop'))               % quit if running nodesktop and not displaying figures.
        if (doSkipFig)
            quit(); 
        else
            fprintf('\n\n----------done-----------\n\n\t\t*****     type the command "quit" to close figures and exit matlab     *****\n\n');
        end
    end
end

function [figname, pdfname, figcount, firstfig] = plot_temps(h, figcount, firstfig, i, nsites, temps, yr1, yr2, yr_50, yr_75, y_end, stnInfo, stns_per_plot, vname, doSkipFig, savePDFs, saveFigs, figdir)

    site_id=stnInfo.stnID(1);
    site=clean_string(stnInfo.stnName(1), 16);
    
    dnums=stnInfo.Properties.UserData.dates;
    day1=stnInfo.Properties.UserData.day1;
    dvecs = datevec(dnums);
    x = dvecs(:,1);
    y=mod(dnums-day1,365)+1;
    
    nr = floor(sqrt(stns_per_plot));
    nc = ceil(stns_per_plot/nr);
    
    if (sum(isnan(temps)) == length(temps))
        minval=0;
        maxval=1;
    else
    
        minval = min(0, min(temps(:)));
        minval1=minval+.1;
        minval2=minval1+.1;
        maxval = floor(max(temps(:))/10)*10 + 10;
    end
    
    mv1=maxval-20;
    mv2=maxval-.01;
    
    xlims = [ceil(-1+yr1/5)*5,floor(yr2/5+1)*5];
    ylims = [-10,380];
    zlims = [minval,maxval];
    
            % calculate the long term trend
    npts = length(temps);
    nyrs = npts/365;
    all_yrs = yr1:yr2;
    t2 = reshape(temps, 365, nyrs);
    good_pts = ~isnan(t2);
    yr_means  = nanmean(t2);
    yr_fracs  = nansum(good_pts) / 365;
    keepers = yr_fracs > .20;
    good_yrs = all_yrs(keepers);
    any_yr_count = sum(yr_fracs > 0);
    day_fracs = nansum(good_pts,2)/any_yr_count;
    
            % get the long term trend
    if (sum(keepers) > 5)

        [p,S,mu] = polyfit(good_yrs,yr_means(keepers),3);
        yr_fit = polyval(p, good_yrs,S,mu);
        do_yr_fit = true;
    else
        do_yr_fit = false;
    end
    
            % get the climatology
    try
        clim = climatology(temps, 5, 2);
        do_clim = true;
    catch
        do_clim = false;
    end
    
    subnum=mod(i-1,stns_per_plot)+1;
    if (stns_per_plot == 1)
        pos = [.01,.95,.6,.6];
        spotsize=4;
    elseif (stns_per_plot <= 4)
        pos = [.01,.95,.8,.8];
        spotsize=2;
    else
        pos = [.01,.95,.8,.8];
        spotsize=1;
    end
    if (firstfig)
        clf;
        if (~doSkipFig)
            set(h,'Units','Normalized');
            set(h,'Position',pos);
        end
        firstfig = false;
        figcount=figcount+1;
    end
    subplot(nr,nc,subnum);
%         surf(x,y,reshape(temps,365,nyrs),'EdgeColor','none');
    scatter3(x,y,temps,spotsize,temps);
    hold on;
    dv = ones(size(temps)) * minval;
    dv(isnan(temps))=nan;
    scatter3(x,y,dv,1,[0,0,0]);
    dv = ones(size(temps)) * minval;
    dv(~isnan(temps))=nan;
    scatter3(x,y,dv,1,[1,.75,.75]);
    
    if (~isempty(yr_50))
        plot3([yr_50,yr_50],[-5,370],[minval2,minval2],'-','linewidth',2,'color',[0,.85,0]);
        yr_50_txt = sprintf('%4d',yr_50);
    else
        yr_50_txt = ' -  ';
    end
    if (~isempty(yr_75))
        plot3([yr_75,yr_75],[-5,370],[minval2,minval2],'r-','linewidth',2);
        yr_75_txt = sprintf('%4d',yr_75);
    else
        yr_75_txt = ' -  ';
    end
    
            % plot long-term trend, fraction of good data
    left_wall = ones(1,nyrs)*ylims(2);
    left_wall2 = ones(1,length(good_yrs))*ylims(2);
    right_wall = ones(1,365)*xlims(2);
    plot3(yr1:yr2, left_wall, ones(1,nyrs)*mv1, 'b-');
    plot3(yr1:yr2, left_wall, ones(1,nyrs)*mv2, 'b-');
    plot3(right_wall, 1:365,  ones(1,365)*mv1,  'b-');
    plot3(right_wall, 1:365,  ones(1,365)*mv2,  'b-');
    plot3(yr1:yr2, left_wall, yr_fracs * 20 + mv1, 'b-');
    plot3(right_wall, 1:365, day_fracs * 20 + mv1, 'b-');
    if (do_yr_fit)
        plot3(good_yrs, left_wall2, yr_fit,'.r-','linewidth',2);
    end
    if (do_clim)
        plot3(right_wall, 1:365, clim,'.r-','linewidth',2);
    end

%             plot3([y_end+1,y_end+1],[-5,370],[tmax,tmax],'k-','linewidth',2);
    hold off;
    grid on;
    xlim(xlims);
    ylim(ylims);
    zlim(zlims);
    
%     if (~isempty(yr_75))
%         view([90,-90]);
%     end
    title(sprintf('%s %s %s %4s %4s %d',site_id, site, vname, yr_50_txt,yr_75_txt, y_end),'Interpreter','none');
    if (subnum == stns_per_plot || i==nsites)
        if (saveFigs || savePDFs)
            if (stns_per_plot < 4)
                set(h,'PaperSize',[16,12],'PaperPosition',[0,0,16,12]);
            else
                set(h,'PaperSize',[22,11.25],'PaperPosition',[0,0,22,11.25]);
            end
            if (stns_per_plot == 1)
                figname = fullfile(figdir, sprintf('%s_%d_%02d_%s.%s.fig',vname,stns_per_plot,figcount, site_id, site));
                pdfname = fullfile(figdir, sprintf('%s_%d_%02d_%s.%s.pdf',vname,stns_per_plot,figcount, site_id, site));
            else
                figname = fullfile(figdir, sprintf('%s_%d_%02d.fig',vname,stns_per_plot,figcount));
                pdfname = fullfile(figdir, sprintf('%s_%d_%02d.pdf',vname,stns_per_plot,figcount));
            end
            if (saveFigs)
                savefig(h,figname);
            else
                figname=[];
            end
            if (savePDFs)
                saveas (h,pdfname);
            else
                pdfname=[];
            end
        else
%            pause(.25);
            figname=[];
            pdfname=[];
        end
        firstfig = true;
    else
        figname=[];
        pdfname=[];
    end
    drawnow();
end

% these are from my ARRM_V2 code library, but I've pasted them in here so this code is standalone.

 function [clim] = climatology(y, nterms, sig_term)        
% clim = climatology(y, SIG) 
%
% calculates the average daily value (365x1) for the data in y  (365*nyrs x 1)
% then low-pass filters it with a filter defined by nterms and sig_term.
%
%   Note that the average here is a straight average for each day of the year (not gaussian-weighted), which is
%   then smoothed by circular convolution along the 365 days.
%
% Returns a 365x1 vector representing the smoothed daily average of the data in y.
%   Inputs:
%       y           data to work with.  must be of multiple of 365 days long
%       nterms      equivalent # of frequency terms to retain
%       sig_term    sigma for gaussian to smoothe rectangular filter.
%                       good values to use for nterms & sig_term are 5 & 2.0
%                       this gives the equivalent filtering of an ideal
%                       filter keeping 1st 6 terms of fft (up to 6 cycles
%                       per year)
%                       for pure gaussian filtering, use nterms=0, sig_term set to desired gaussian sigma.  
%frequency domain sigma to smooth with.  See "math notes" for equations to go from
%                   either time-domain gaussian sigma or idea-filter-equivalent sigma.

        % if keeping all terms, just do straight average:
    nyrs = length(y)/365;    
    clim = nanmean(reshape(y,365,nyrs),2);     % get average daily value for data range.

    if (sum(isnan(clim))>0)
        throw(MException('ICSF:CLIMNANS', sprintf('climatology:  Error:  %d days have no valid readings', sum(isnan(clim)))));
    end
    
        % Create the filter, and low-pass filter the daily averages
    FILT = calc_filter(nterms, sig_term, 1);
    [clim, ~] = lpf_FILT(clim, FILT);

end
function [FILT, filt] = calc_filter(nterms, sig_term, myrs)
%FILT = calc_filter(nterms, sig_term, myrs)
%
%   Returns fourier domain filter for an ideal rectangular filter with nterms
%   convolved with gaussian of sigterm, for myrs*365 days.
%   Filter is based on a 365-day year, so an nterms of 5 with sig_term of 3 means start from
%   an ideal filter passing frequencies of up to 5 cycles per 365-day year,
%   then convolve it with a guassian with a sigma of 3.
%
%   Inputs:
%       nterms          # of terms per year to keep 
%                           (positive terms only, not including DC or negative terms.)
%                           nterms=5 starts with a rect. filter of length 11
%                           nterms is relative to 1 year.  Will be scaled up by mterms.
%       sig_term        sigma for gaussian to convolve with rect. filter
%                           s/b specified relative to 1 year.  Will be scaled up by mterms
%       myrs            # of years to create filter for. 
%                           (1, 41, 141, etc., for 1 year, or 41 yrs of obs data, or 141 yrs of model data...)
%
%   Returns:
%       FILT            Fourier domain filter of length 365*myrs.
%       filt            circular spatial domain filter of length 365*myrs.
%                           note:  filt is centered at 1, intended for use
%                           with cconv(...), circular convolution.
%                           Centering at 1 means filter does not shift
%                           signal.
%                           filt is not calculated if only 1 return
%                           variable is detected.


        % the rectangular (ideal) low-pass filter
    len = 365*myrs;
    nterms = floor(nterms*myrs);
    sig = sig_term*myrs; 
    Frect = zeros(1,len);

    Frect([1:(nterms+1), (end-nterms+1):end]) = 1;
    
                % the gaussian
    midpt = ceil((len+1)/2);        % 11: -> 6.  12:  ->7.
    GF=gauss(len,sig, midpt);
    GF=circshift(GF,1-midpt,2);     % shift peak back to GF(1).
                % now convolve Frect with GF to make filter
    FILT = cconv(Frect,GF,len);
    FILT = FILT/FILT(1);          % make sure that DC is passed at full power.
    
    if (nargout > 1)
        filt = real(ifft(FILT));    % *should* be purely real, but may have tiny fractional imag. part due to computational limitations...~ 10^-15.
    else
        filt = nan;
    end
end

function [yout, YOUT] = lpf_FILT(y, FILT)

%   returns y circularly filtered with fourier-domain filter FILT
%   also returns Y_OUT, the fourier transform of y_out.
%

        % make sure y & filt are same shape...
    if ((iscolumn(y) && isrow(FILT)) || (isrow(y) && iscolumn(FILT)))
        FILT = FILT';
    end
    
    Y=fft(y);
    YOUT = Y .* FILT;
    yout = real(ifft(YOUT));
end

function valids = valid_counts(vals)

    st_day = [1, 32, 62, 92, 122, 153, 183, 214, 244, 274, 305, 335,366];

    valids = zeros(1,12);
    d = ~isnan(vals);
    nyrs = length(d)/365;
    d = reshape(d, 365, nyrs);
    for j=1:12
        s = st_day(j);
        e = st_day(j+1)-1;
        valids(j) = sum(sum(d(s:e,:)));
    end
end


function [ok, valids, ord_txt, validFlag] = check_valid(temps, yr1, yr2)
                % check to see if this site fails for either NAs or trend order.
    ok = true;
    valids  = valid_counts(temps);
    ord     = ARRM_V2_calc_order([yr1,yr2], temps, .5);
    validFlag = '';
    if (any(valids<600))
        validFlag = 'NA****';
        ok=false;
    end
    if (isnan(ord))
        ord_txt = '-';
        ok=false;
        if (isempty(validFlag))
            validFlag = '***ORD';
        else
            validFlag(4:6)='ORD';
        end
    else
        ord_txt = sprintf('%d',ord);
    end
    if (ok)
        validFlag = '  ok  ';
    end
end

function stns = make_stn_fail_list(stnTbl, yr1, yr2)
% returns a list of stations which fail for either NAs or trend order.

    allStns=stnTbl.stnID;
    keepers = true(size(allStns));
    
    for i=1:length(allStns)
        stn=allStns(i);
        stnInfo = QC_get_site_dadta(stnTbl, yr1,yr2, "stnID", stn, "removeLeaps",true, "searchType","stnID");
        stnInfo = stnInfo(1,:);
        temps = stnInfo.data(1,:);
%         sdate = datevec(stnInfo.startDate);
%         edate = datevec(stnInfo.endDate);
%         syr = max(yr1, sdate(1));
%         eyr = min(yr2, edate(1));
        ok = check_valid(temps, yr1, yr2);
        keepers(i) = ~ok;
    end
    stns = allStns(keepers);
end
