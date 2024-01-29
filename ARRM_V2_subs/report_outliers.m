function report_outliers(DAs, yr_types, lbls, DSP, runLbl)
%   Reports outliers and far_outliers for 1 or more DAs.
%
%   DAs         1 or cell array of DAs
%   yr_types    1 or cell array of strings "base_yrs", "hist_yrs", "rolling_yrs", 1 for each DA
%   lbls        string labels to put on each outler
%   DSP         DP or DSP for run disaggregation
%   runLbl      addition label
%                             report_outliers({DA_obs_out,                          DA_results}, {'base_yrs',                         'rolling_yrs'}, ["obs","map"],DSP_base, runLbl);


    if (~iscell(DAs)),        DAs      = {DAs};             end
    if (~isstring(yr_types)), yr_types = string(yr_types);  end
    if (~isstring(lbls)),     lbls     = string(lbls);      end
    
    if (DAs{1}.DP.isPrecipRun)
        max_outliers = 3000;
    else
        max_outliers = 300;
    end
    
    for k=1:length(DAs)
        if (isempty(DAs{k})), continue; end
        if (strcmp(yr_types{k},'hist_yrs') && ~DSP.separate_hist), continue; end  % skip hist outliers if hist is same data as model
        outliers = find(DAs{k}.outlier_map(:));
        far_outliers = find(DAs{k}.far_outlier_map(:));     
        nouts = length(outliers);
        nfars = length(far_outliers);
        if (DAs{k}.RP.isPrecipRun)
            nvalid = numel(DAs{k}.anoms) - sum(DAs{k}.zeros_map);
        else
            nvalid = sum(~isnan(DAs{k}.anoms));
        end
        far_outflags = false(length(outliers),1);
        if (nfars > 0)
            for i=1:nfars
                ix=find(outliers == far_outliers(i));
                far_outflags(ix) = true; %#ok<FNDSB>
            end
        end
        if (~isempty(DAs{k}.orig_anoms))
            orig_anoms = DAs{k}.orig_anoms(outliers);  
        else
            orig_anoms = DAs{k}.anoms(outliers);
        end
        if (~strcmp(yr_types{k},'rolling_yrs'))     % obs & hist
            anom_probs = DAs{k}.anom_probs_base(outliers);
            sdev_equiv = DAs{k}.anom_sdevs_base(outliers);
        else
            anom_probs = DAs{k}.anom_probs_rolling(outliers);
            sdev_equiv = DAs{k}.anom_sdevs_rolling(outliers);
        end
        yrlen = DAs{k}.DP.yrlen;
        nyrs = length(outliers)/yrlen;
        [doy, yr] = ind2sub([yrlen, nyrs], outliers);
        DSP.print_log( "%s %s Outliers:  %d outliers of %d %.2f %%, %d far outliers %.2f %%\n", runLbl, DAs{k}.DAType, nouts, nvalid, nouts/nvalid*100, nfars, nfars/nvalid*100);
        if (nouts > max_outliers)
            DSP.warn_log("*****---Warning: report_outliers(): %s : Too many outliers to report: %d -----*****\n", DSP.stnName, nouts);
        else
            if (~strcmp(lbls{k},'map'))
                report_outliers_sub(DSP, DAs{k}.DAType, lbls(k), yr, doy, outliers, far_outflags, orig_anoms, anom_probs, sdev_equiv, [],     [],        [],              max_outliers, nvalid);
            else
                mapped = DAs{k}.mapped_anoms(outliers);
                mapped_im = DAs{k}.mapped_anoms_im(outliers);
                adjusted_mapped = DAs{k}.outlier_mapped_anoms(outliers);
                report_outliers_sub(DSP, DAs{k}.DAType, lbls(k), yr, doy, outliers, far_outflags, orig_anoms, anom_probs, sdev_equiv, mapped, mapped_im, adjusted_mapped, max_outliers, nvalid);
            end
        end
        DSP.print_log('\n');
    end
end

function report_outliers_sub(DSP, DAType,         lbl_in, yr, doy, outliers, far_outflags, orig_anoms, anom_probs, sdev_equiv, mapped, mapped_im, adjusted_mapped, max_outliers, nvalid)
    
%        report_outliers_sub(DSP, DAs{k}.DAType, lbls(k), yr, doy, outliers, far_outflags, anoms, anom_probs, sdev_equiv, mapped, mapped_im, adjusted_mapped, max_outliers, nvalid);

    nouts = length(outliers);
    nfars = sum(far_outflags); 
    DSP.print_log( "%s Outliers:  %d outliers of %d %.2f %%, %d far outliers %5.2f %%\n", DAType, nouts, nvalid, nouts/nvalid*100, nfars, nfars/nvalid*100);
    if (nouts > max_outliers)
        DSP.warn_log("*****---Warning: report_outliers_sub(): %s : Too many outliers to report: %d-----*****\n", DSP.stnName, nouts);
    else
        if (~exist('adjusted_mapped','var') || isempty(adjusted_mapped))            
            DSP.print_log( " Outliers         Year  Day  index  probability sdev_eq anomaly\n");
            ismapped = false;
        else
            ismapped = true;
            DSP.print_log( " Outliers         Year  Day  index  probability sdev_eq orig_anom mapped_to   diff     mapval adj_mapval    diff\n");
            difs1 = mapped - orig_anoms;
            difs2 = adjusted_mapped - mapped_im;
        end
        for i=1:nouts
            if (far_outflags(i))
                lbl = sprintf("far_%s_outlier", lbl_in);
            else
                lbl = sprintf(    "%s_outlier", lbl_in);
            end
            DSP.print_log('%16s: %4d  %3d %6d  %11.8f  %6.3f  %6.2f ', lbl, yr(i), doy(i), outliers(i), anom_probs(i), sdev_equiv(i), orig_anoms(i));
            if (ismapped)
                DSP.print_log( '%10.2f %7.2f %10.2f %10.2f %7.3f\n', mapped(i), difs1(i), mapped_im(i), adjusted_mapped(i), difs2(i));
            else
                DSP.print_log( '\n');
            end

        end
    end
end

