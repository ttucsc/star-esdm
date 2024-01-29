function NCAT = ARRM_V2_NC_Attributes_TTUCSC(runType)
%
%       file to provide general metadata for the netcdf output files
%       processed on this system.
%
%       words in upper case, such as VARNAME  MDLSTARTYEAR etc. are
%       keywords and will be replaced with the appropriate info at run
%       time.
%
%     NC_Attributes = struct;
    
    [~, longname, hostname] = getusername();
    NCAT.hostname = hostname;
    NCAT.CreatedBy  = longname;
    NCAT.institution  = 'Texas Tech Climate Science Center';
    NCAT.CreationDate = string(datestr(now,'yyyy-mm-dd HH:MM:SS'));

    if (any(contains(runType,'station')))
        NCAT.long_title    = "ARRM V2 station Downscaling, VARNAME, MDLSTARTYEAR - MDLENDYEAR, STNID : STNNAME, ( LAT , LON )";
        NCAT.comments      = "Texas Tech Climate Science Center downscaling run, VARNAME, MDLSTARTYEAR - MDLENDYEAR, STNID : STNNAME  ( LAT , LON )";
    else
        NCAT.title         = "ARRM_V2_gridded_downscaling_VARNAME";
        NCAT.long_title    = "ARRM V2 gridded Downscaling, VARNAME, MDLSTARTYEAR - MDLENDYEAR ";
        NCAT.comments      = "Texas Tech Climate Science Center downscaling run, VARNAME, MDLSTARTYEAR - MDLENDYEAR, STNID : STNNAME  ( LAT1, LON1) to ( LATEND, LONEND)";
    end
    if (any(contains(runType,'precip')))
        % add comments specific to precip runs here.
    end
end

