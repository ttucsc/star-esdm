function DSP = ARRM_V2_DownscalingParams_TTUCSC(runType)

    [~,hostname]=system('hostname');
    hostname = strtrim(hostname);

    on_laptop      = contains(hostname,"jcsf") || contains(hostname,"jpro");
    on_dev_system  = contains(hostname,"neys");
    on_hpcc_system = contains(hostname,"compute");
    
%..................Parameters that depend on station/gridded and on which system we're running on

    if (any(contains('station',lower(runType))))
        if (on_dev_system)
            DSP = downscaling_params_station_dev(on_laptop);
        elseif (on_hpcc_system)            
            DSP = downscaling_params_station_hpcc();
        else
            DSP = downscaling_params_station_KMQ();
        end
        isStation=true;
    else
        if (on_dev_system)
            DSP = downscaling_params_gridded_dev(on_laptop);
        elseif (on_hpcc_system)            
            DSP = downscaling_params_gridded_hpcc();
        else
            DSP = downscaling_params_gridded_KMQ();
        end
        isStation=false;
    end
    DSP.NC_atts = NC_attributes(isStation);

end

function NC_atts = NC_attributes(isStation)

%          CF Conventions Global Attributes. 
%          Modify as appropriate.  You can override these in 

%   Required Meta-data for output file:
    if (~isStation)
        NC_atts.institution   = "Texas Tech Climate Center, Texas Tech University, Lubbock, TX";

        NC_atts.title         = "STAR-ESDM_gridded_downscaling_VARNAME";
        NC_atts.long_title    = "STAR-ESDM gridded Downscaling, VARNAME, MDLSTARTYEAR - MDLENDYEAR ";
        NC_atts.comments      = "Texas Tech Climate Science Center downscaling run, VARNAME, MDLSTARTYEAR - MDLENDYEAR";
        NC_atts.date_range    = "MDLSTARTYEAR - MDLENDYEAR";
        NC_atts.creation_date = "CURDATE";
    %   NC_attributes.references    = "references to any papers or documentation";   
        NC_atts.station_data  = "no";
        NC_atts.gridded_data  = "yes";
        NC_atts.downscaled    = "yes";
        NC_atts.references    = "";
        NC_atts.history       = "";   % this will appended to later.  Add any info you need to start the history text
        

    else    
        NC_atts.institution   = "Texas Tech Climate Center, Texas Tech University, Lubbock, TX";

        NC_atts.title         = "STAR-ESDM_station_downscaling_VARNAME";
        NC_atts.long_title    = "STAR-ESDM station Downscaling, VARNAME, MDLSTARTYEAR - MDLENDYEAR";
        NC_atts.comments      = "Texas Tech Climate Science Center Station downscaling run, VARNAME, MDLSTARTYEAR - MDLENDYEAR";
        NC_atts.date_range    = "MDLSTARTYEAR - MDLENDYEAR";
        NC_atts.creation_date = "CURDATE";
    %   NC_attributes.references    = "references to any papers or documentation";   
        NC_atts.station_data  = "yes";
        NC_atts.gridded_data  = "no";
        NC_atts.downscaled    = "yes";
        NC_atts.references    = "";
        NC_atts.history       = "";   % this will appended to later.  Add any info you need to start the history text
    end    
end

function [DSP] =  downscaling_params_station_dev(on_laptop)

    
%       Data parameters
%
%           this is the default. 
    DSP.region = "ncamerica";     % default

%           Data directories
    
    % output file & folder 
%   Set outdir in DataParams_TTUCSC, not here.


    if (~on_laptop)
        DSP.obsdir       = "/Volumes/lacie_1/data/obs/stations_netcdf";
        DSP.histdir      = "/Volumes/lacie_1/data/gcm_cmip5/rotated_netcdf";
        DSP.mdldir       = "/Volumes/lacie_1/data/gcm_cmip5/rotated_netcdf";
    else
        DSP.obsdir       = "/Volumes/2018_1/data";
        DSP.histdir      = "/Volumes/2018_1/data";
        DSP.mdldir       = "/Volumes/2018_1/data";        
    end    
%           Filenames 
%   DataParams.obsnames    = "stations.OBSVNAME.REGION.1850.2017.nc";
    DSP.obsnames    = []; % DataParams.stnInfo.Properties.UserData.ncName;
    DSP.default_obsnames = "stations.OBSVNAME.REGION.1850.2018.20_yrs_min.nc";
    DSP.histnames = [];
    DSP.mdlnames  = [];
    DSP.default_mdlnames = ["MODEL.ENSEMBLE.VARNAME.hist.day.1900.2005.llt.nc", "MODEL.ENSEMBLE.VARNAME.SCENARIO.day.2006.2100.llt.nc"];

    % output file  
%   outname:  set this in DataParams_TTUCSC, not here.

end

%----------------------------------------------------

function [DSP] =  downscaling_params_gridded_dev(on_laptop)

%          CF Conventions Global Attributes. 
%          Modify as appropriate.  You can override these in 

%       Data parameters
%
%           this is the default. 
    DSP.region = "ncamerica";     % default

%           Data directories
    
    % output file & folder 
%   Set outdir in DataParams_TTUCSC, not here.

    if (~on_laptop)
        DSP.obsdir       = "/Volumes/lacie_1/data/gcm_cmip5/rotated_netcdf";
        DSP.histdir      = "/Volumes/lacie_1/data/gcm_cmip5/rotated_netcdf";
        DSP.mdldir       = "/Volumes/lacie_1/data/gcm_cmip5/rotated_netcdf";
    else
        DSP.obsdir       = '/Volumes/2018_1/data';
        DSP.histdir      = "/Volumes/2018_1/data";
        DSP.mdldir       = "/Volumes/2018_1/data";        
    end    
    
%           Filenames 
    DSP.obsnames     = [];
    DSP.default_obsnames = "MODEL.ENSEMBLE.VARNAME.hist.day.1900.2005.llt.nc";
    DSP.histnames    = [];
    DSP.mdlnames     = [];
    DSP.default_mdlnames = ["MODEL.ENSEMBLE.VARNAME.hist.day.1900.2005.llt.nc", "MODEL.ENSEMBLE.VARNAME.SCENARIO.day.2006.2100.llt.nc"];

    % output file  
%   outname:  set this in DataParams_TTUCSC, not here.
   
end

%----------------------------------------------------

function [DSP] =  downscaling_params_station_hpcc()

%          CF Conventions Global Attributes. 
%          Modify as appropriate.  You can override these in 
    
%       Data parameters
%           this is the default. 
    DSP.region = "ncamerica";     % this needs to be passed in.
%           Data directories
    
    % output file & folder 
%   Set outdir in DataParams_TTUCSC, not here.

    DSP.obsdir       = "";   % LEAVE THIS BLANK!  the filename in stnInfo contains pathname.
    DSP.obsdir       = "/lustre/scratch/iscottfl/obs/stations_netcdf";
    DSP.histdir      = "/lustre/scratch/iscottfl/cmip5_rotated";
    DSP.mdldir       = "/lustre/scratch/iscottfl/cmip5_rotated";
 
%           Filenames 
%   DataParams.obsnames    = "stations.OBSVNAME.REGION.1850.2017.nc";
    DSP.obsnames    = []; % DataParams.stnInfo.Properties.UserData.ncName;
    DSP.default_obsnames = "stations.OBSVNAME.REGION.1850.2018.20_yrs_min.nc";
    DSP.histnames   = [];
    DSP.mdlnames     = [];
    DSP.default_mdlnames = ["MODEL.ENSEMBLE.VARNAME.hist.day.1900.2005.llt.nc", "MODEL.ENSEMBLE.VARNAME.SCENARIO.day.2006.2100.llt.nc"];
    % output file  
%   outname:  set this in DataParams_TTUCSC, not here.

end

%----------------------------------------------------

function [DSP] =  downscaling_params_gridded_hpcc()

%          CF Conventions Global Attributes. 
%          Modify as appropriate.  You can override these in 

%       Data parameters
%
%           Data directories
    
    % output file & folder 
%   Set outdir in DataParams_TTUCSC, not here.

    DSP.obsdir       = "/lustre/scratch/iscottfl/cmip5_rotated";   % this needs fixing, Ian!  this is just so I can test with one model against another.
    DSP.histdir      = "/lustre/scratch/iscottfl/cmip5_rotated";
    DSP.mdldir       = "/lustre/scratch/iscottfl/cmip5_rotated";
    
%           Filenames 
    DSP.obsnames     = "CCSM4.r1i1p1.tasmax.hist.day.1900.2005.llt.nc";  % this is just hardcoded to use CCSM4 as my "obs" data for testing.
    DSP.histnames    = [];
    DSP.mdlnames     = [];
    DSP.default_mdlnames = ["MODEL.ENSEMBLE.VARNAME.hist.day.1900.2005.llt.nc", "MODEL.ENSEMBLE.VARNAME.SCENARIO.day.2006.2100.llt.nc"];


    % output file  
%   outname:  set this in DataParams_TTUCSC, not here.

end

%----------------------------------------------------

function [DSP] =  downscaling_params_station_KMQ()


%       Data parameters
%           Data directories
    
    % output file & folder 
%   Set outdir in DataParams_TTUCSC, not here.

    DSP.obsdir       = "/data/obs/stations_netcdf";   % LEAVE THIS BLANK!  the filename in stnInfo contains pathname.
    DSP.histdir      = "/data/gcm_cmip5/daily.1900.2100/VARNAME/";
    DSP.mdldir       = "/data/gcm_cmip5/daily.1900.2100/VARNAME/";

    
%           Filenames 
%   DataParams.obsnames     = "stations.OBSVNAME.REGION.1850.2017.nc";
    DSP.obsnames     = []; % DataParams.stnInfo.Properties.UserData.ncName;
    DSP.default_obsnames = "stations.OBSVNAME.REGION.1850.2018.20_yrs_min.nc";
    DSP.histnames    = [];
    DSP.mdlnames     = [];
    DSP.default_mdlnames = ["MODEL.ENSEMBLE.VARNAME.hist.day.1900.2005.llt.nc", "MODEL.ENSEMBLE.VARNAME.SCENARIO.day.2006.2100.llt.nc"];

    % output file  
%   outname:  set this in DataParams_TTUCSC, not here.

end
   
%----------------------------------------------------

function [DSP] =  downscaling_params_gridded_KMQ()

%          CF Conventions Global Attributes. 
%          Modify as appropriate.  You can override these in 

%       Data parameters
%
%           Data directories
    
    % output file & folder 
%   Set outdir in DataParams_TTUCSC, not here.

    DSP.obsdir       = "/data/gcm_cmip5/daily.1900.2100/VARNAME/hist";   % this needs fixing, Ian!  this is just so I can test with one model against another.
    DSP.histdir      = [];
    DSP.mdldir       = "/data/gcm_cmip5/daily.1900.2100/VARNAME/";

    
%           Filenames 
    DSP.obsnames     = "CCSM4.r1i1p1.tasmax.hist.day.1900.2005.llt.nc";  % this is just hardcoded to use CCSM4 as my "obs" data for testing.
    DSP.histnames    = [];
    DSP.mdlnames     = [];
    DSP.default_mdlnames = ["MODEL.ENSEMBLE.VARNAME.hist.day.1900.2005.llt.nc", "MODEL.ENSEMBLE.VARNAME.SCENARIO.day.2006.2100.llt.nc"];
   
    % output file  
%   outname:  set this in DataParams_TTUCSC, not here.
end

%----------------------------------------------------


