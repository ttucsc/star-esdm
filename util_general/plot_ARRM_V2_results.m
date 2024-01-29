function [fignum, nplots] = plot_ARRM_V2_results(plotFlags, fignum, DA_obs, DA_hist, DA_mdl, DA_results, runLbl)
%[fignum, nplots] = plot_ARRM_V2_results(plotFlags, fignum, raw_data, results, runLbl)
%
%   Plots data from ARRM_V2_SignalDisaggration's
%
%   plotFlags:   
%       0       don't do any plots
%       1       (scalar):  do all plots
%
%       (array of 7 boolean flags):
%         Flag
%           1   plot or raw_data, raw time series
%           2   Climatologies (obs, hist & mdl)
%           3   Running Climatology as time series
%           4   Running Climatology as surface
%           5   Anomaly Surface
%           6   PDF surfaces
%           7   PDF surfaces
%
%   This has been kludged around several times, and should be rewritten to work properly with DA objects, Ian!
%
   
    nplots=0;
    if (isempty(plotFlags) || ~any(plotFlags)), return; end
    
    if (length(plotFlags)==1), plotFlags = true(1,7); end
    
    yrlen = DA_results.RP.yrlen;
%    mapped = results.mapped;

    
    obsyrs = length(DA_obs.raw_data)/yrlen;
%   histyrs= length(raw_data{Hist}) /yrlen;
    mdlyrs = length(DA_mdl.raw_data)/yrlen;
    xo=(1:length(DA_obs.raw_data))  /yrlen + DA_obs.DP.base_yrs(1);
    xm=(1:length(DA_mdl.raw_data))  /yrlen + DA_mdl.DP.rolling_yrs(1);
    
%   mclim_obs  = reshape(moving_clim{Obs},yrlen,obsyrs);
%   mclim_hist = reshape(moving_clim{Hist},yrlen,histyrs);
    mclim_mdl  = reshape(DA_mdl.moving_clim, yrlen, mdlyrs);
%   mdl_bins = DA_mdl.RP.bins;
    data_yrs = DA_mdl.DP.data_yrs;
    data_yrs = data_yrs(1):data_yrs(2);
    nyrs = data_yrs(end)-data_yrs(1)+1;
    
        % in case any of the data is in Kelvin, rather than Celcius

    
    if (nanmean(DA_obs.raw_data) >150),  obs_off  = 273.15; else, obs_off = 0;  end
    if (nanmean(DA_hist.raw_data)>150),  hist_off = 273.15; else, hist_off = 0; end
    if (nanmean(DA_mdl.raw_data) >150),  mdl_off  = 273.15; else, mdl_off = 0;  end
    
    if (plotFlags(1))

        figure(fignum);
        nplots=nplots+1;
        fignum=fignum+1;
        
            % plot obs & model data.
        
        subplot(2,1,1);
        plot(xo, DA_obs.raw_data-obs_off, xm,DA_mdl.raw_data - mdl_off);
        grid on;
        legend('obs','model')
        title(sprintf('raw time series, %s', runLbl))
        subplot(2,1,2);
        rd1 = reshape(DA_obs.raw_data,yrlen,obsyrs);
        surf(rd1,'edgecolor','none');
        xlabel('year');
        ylabel('day of year')
        title(sprintf('Observations as surface, %s',runLbl));
        
    end        
    
    if (plotFlags(2))

        figure(fignum);
        nplots=nplots+1;
        fignum=fignum+1;
        
        % get the hist raw data
        [hstart, hend] = DA_hist.using_range('base_yrs');
        hist_raw = DA_hist.raw_data(hstart:hend);
        
        % Calculate the phase shift between Hist and Obs, and a good place to plot it.
        
        h_phase_shift = DA_hist.phase_clim - DA_obs.phase_clim;
%       fprintf("phase shift %.2f\n", h_phase_shift);
        % plot the climatology

        clmean_obs = nanmean(reshape(DA_obs.raw_data,yrlen,[]),2);
        clmean_hist = nanmean(reshape(hist_raw,yrlen,[]),2);
        plot(1:yrlen,clmean_obs-obs_off,'c-');
        hold on;
        grid on;
        plot(1:yrlen, clmean_hist-hist_off,'-','color',[1,.5,.5])
        hh = plot(1:yrlen,DA_obs.base_clim-obs_off, 'b-',1:yrlen,DA_hist.base_clim-hist_off,'r-', "linewidth",2); %,1:yrlen),DA_mdl.base_clim-mdl_off);
        yl = [min(DA_obs.base_clim-obs_off), max(DA_obs.base_clim-obs_off)];
        dyl = range(yl);
        yl1 = yl + .05*dyl*[-1,1];
        plot([DA_obs.phase_clim,DA_obs.phase_clim], yl1,'b-');
        hh2 = plot([DA_hist.phase_clim, DA_hist.phase_clim],yl1, 'r-');
        legend([hh;hh2], 'obs clim','mdl(hist) clim',sprintf("phase shift %.1f days", h_phase_shift)); %,'model clim')
                
        hold off;
        title(sprintf('Climatologies, %s', runLbl));
                
    end
    

    if (plotFlags(3))

        figure(fignum);
        nplots=nplots+1;
        fignum=fignum+1;
        
        obs_clim = DA_obs.moving_clim_1D(:) + DA_obs.trend;
        mdl_clim = DA_mdl.moving_clim_1D(:) + DA_mdl.trend;
        
%       plot(xo,DA_obs.moving_clim(:), xo, DA_obs.trend);
        hh = plot(xo,obs_clim, 'b-', xm, mdl_clim,'r-');
        hold on;
%       plot(xm, DA_mdl.moving_clim(:), xm, DA_mdl.trend);
        plot(xo, DA_obs.trend, 'b-', xm, DA_mdl.trend),'r-';
        hold off;
        grid on;
        legend(hh, 'obs','mdl');
        title(sprintf('running climatology, %s', runLbl));

    end
    
    if (plotFlags(4))

        figure(fignum);
        nplots=nplots+1;
        fignum=fignum+1;
        
        surf(data_yrs,1:yrlen, mclim_mdl,'edgecolor','none');
        hold on;
        mz = repmat(max(mclim_mdl(:)), nyrs,1);
        plot3(data_yrs, DA_mdl.phase_moving,mz,'r-');
        hold off;
        lights_on;
        xlabel('year');  
        ylabel('day of year');
        title(sprintf('Model climatology surface, %s and phase', runLbl));
                
    end
    
    if (plotFlags(5))

        figure(fignum);
        nplots=nplots+1;
        fignum=fignum+1;
        
        pmax = max([max(DA_obs.pdf_base(:)), max(DA_hist.pdf_base(:)), max(DA_mdl.pdf_base(:)), max(DA_mdl.pdf_rolling(:))]);

        subplot(3,1,1)
        surf(DA_obs.pdf_base,'edgecolor','none');
        zlim([0, pmax]);
        lights_on;
        xlabel('day of year');
        ylabel('anomaly (deg C)');
        title(sprintf('Obs Anomalies, %s', runLbl));
        subplot(3,1,2);
        surf(DA_hist.pdf_base,'edgecolor','none');
        zlim([0, pmax]);
        lights_on;
        xlabel('day of year');
        ylabel('anomaly (deg C)');
        title('Hist Anomalies');
        subplot(3,1,3);
        for i=1:size(DA_mdl.pdf_rolling,3)
            surf(DA_mdl.pdf_rolling(:,:,i),'edgecolor','none');
            zlim([0, pmax]);
            lights_on;
            xlabel('day of year');
            ylabel('anomaly (deg C)');
            title('Model Anomalies');
            pause(.25);
        end
%       lights_on;
    end
    
        % these need updating, Ian.
    
%     if (plotFlags(6) || plotFlags(7))
%         DA_results.disaggregate()
%     end
%     
%     if (plotFlags(6))
% 
%         figure(fignum);
%         nplots=nplots+1;
%         fignum=fignum+1;
%         
%         pdf_mapped = jc_kdm_anom(anom_mapped,   yedges, sigrange);
%         figure(6);
%         subplot(4,1,1);
%         surf(1:yrlen, yvals, pdf{Hist},'edgecolor','none');
%         lights_on;
%         xlabel('day of year');
%         ylabel('anomaly (deg C)');
%         title(sprintf('Hist, %s', runLbl));
%         subplot(4,1,4);
%         surf(1:yrlen, yvals, pdf_mapped,'edgecolor','none');
%         lights_on;
%         xlabel('day of year');
%         ylabel('anomaly (deg C)');
%         for i=1:size(pdf{Model},3)
%             subplot(4,1,2);
%             surf(1:yrlen, yvals, pdf{Model}(:,:,i),'edgecolor','none');
%             lights_on;
%             xlabel('day of year');
%             ylabel('anomaly (deg C)');
%             title('Model')
%             subplot(4,1,3);
%             surf(1:yrlen, yvals, pdf{Model}(:,:,i) - pdf{Hist},'edgecolor','none');
%             lights_on;
%             xlabel('day of year');
%             ylabel('anomaly (deg C)');
%             title('Model - Hist')
%             pause(.25);
%         end
% 
%     end
%     
%     if (plotFlags(7))
% 
%         figure(fignum);
%         nplots=nplots+1;
%         fignum=fignum+1;
%         
%         figure(7);
%         subplot(3,1,1);
%         surf(1:yrlen, yvals, pdf{Obs},'edgecolor','none');
%         lights_on;
%         xlabel('day of year');
%         ylabel('anomaly (deg C)');
%         title(sprintf('Obs, %s', runLbl));
%         
%         subplot(3,1,2);
%         surf(1:yrlen, yvals, pdf_mapped,'edgecolor','none');
%         lights_on;
%         xlabel('day of year');
%         ylabel('anomaly (deg C)');
%         title('Model');
%         
%         subplot(3,1,3);
%         surf(1:yrlen, yvals, pdf{Obs}-pdf_mapped,'edgecolor','none');
%         lights_on;
%         xlabel('day of year');
%         ylabel('anomaly (deg C)');
%         title('Difference (Obs - mapped)')   
%     end
%     

end

