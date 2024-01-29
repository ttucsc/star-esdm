function RP = ARRM_V2_RunParams_TTUCSC(runType)

    RP = struct;
    
    if (any(contains(string(runType),'precip')))
        RP.do_fixNAs = false;
        RP.trend_order = [];
        RP.trend_yr_flags = [];
        RP.trend_thresh = [];
    elseif (any(contains(string(runType),'temp')))
        
    end
    
    if (any(contains(string(runType),'station')))
        
    end
    
end

