function report_trends(DAs, DSP, runLbl)

    ndas = length(DAs);
    DSP.print_log( 'trends:\n');
    for k=1:ndas
        if (isempty(DAs{k})), continue; end
        nused = sum(DAs{k}.RP.trend_yr_flags);
        nyrs  = length(DAs{k}.RP.trend_yr_flags);
        pct   = nused / nyrs * 100.0;
        if (isnan(DAs{k}.RP.trend_order))
            DSP.print_log('\t%-40s %-11s order TREND_FAIL : only %3d of %3d yrs ( %5.1f%% )\n', runLbl, DAs{k}.DA_yrlbl,                         nused, nyrs, pct);
        else  
            DSP.print_log('\t%-40s %-11s order       %3d  : used %3d of %3d yrs ( %5.1f%% )\n', runLbl, DAs{k}.DA_yrlbl, DAs{k}.RP.trend_order,  nused, nyrs, pct);
        end
    end

    DSP.print_log( '\n');
end

