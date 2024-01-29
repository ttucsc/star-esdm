function report_counts(DAs, yr_types, DSP,runLbl)

    for k=1:length(DAs)
        if (isempty(DAs{k})), continue; end
        [ok,  month_ok, npts, nas_in_counts, outlier_counts, far_outlier_counts, valid_counts, yr_range] = DAs{k}.check_nans(yr_types{k});

                % output the monthly valid counts, NA_counts and added_na counts.
        DSP.print_log('\nmonthly counts: %11s %-11s           Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec  total', DAs{k}.DAType, yr_types{k});
        
         if (~ok)
            DSP.print_log(' *** insufficient data to downscale\n');
         else
            DSP.print_log( '\n');
         end

         nsets = size(npts,2);

         for j=1:nsets
             month_flags = repmat('*',12,1);
             month_flags(month_ok) = ' ';
             pct_valid = valid_counts(:,j)./npts(:,j)*100;
             tot_pct_valid = sum(valid_counts(:,j))/sum(npts(:,j))*100;
             DSP.print_log( 'nas_in:       %11s %11s %4d %4d ', DAs{k}.DAType, yr_types{k}, yr_range);
             for i=1:12
                DSP.print_log( ' %4d %c', nas_in_counts(i,j), month_flags(i));
             end
             DSP.print_log( '%5d    %s\n', sum(nas_in_counts(:,j)), runLbl);
             DSP.print_log( 'all_outliers: %11s %11s %4d %4d ', DAs{k}.DAType, yr_types{k}, yr_range);
             DSP.print_log( '%5d  ', outlier_counts(:,j));            
             DSP.print_log( '%5d    %s\n', sum(outlier_counts(:,j)), runLbl);
             DSP.print_log( 'far_outliers: %11s %11s %4d %4d ', DAs{k}.DAType, yr_types{k}, yr_range);
             DSP.print_log( '%5d  ', far_outlier_counts(:,j));
             DSP.print_log( '%5d    %s\n', sum(far_outlier_counts(:,j)), runLbl);
             DSP.print_log( 'npts:         %11s %11s %4d %4d ', DAs{k}.DAType, yr_types{k}, yr_range);
             DSP.print_log( '%5d  ', npts(:,j));
             DSP.print_log( '%5d    %s\n', sum(npts(:,j)), runLbl);
             DSP.print_log( 'pct_valid:    %11s %11s %4d %4d ', DAs{k}.DAType, yr_types{k}, yr_range);
             DSP.print_log( '%5.1f  ', pct_valid);
             DSP.print_log( '%5.1f    %s\n', tot_pct_valid, runLbl);
             DSP.print_log( '\n');
         end
    end  
end

