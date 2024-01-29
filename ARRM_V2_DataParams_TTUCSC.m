function DP = ARRM_V2_DataParams_TTUCSC(runType)

    [~,hostname]=system('hostname');
    hostname = strtrim(hostname);

    on_jmac        = contains(hostname,"jmac");
    on_laptop      = contains(hostname,"jcsf") || contains(hostname,"jpro");
    on_dev_system  = contains(hostname,"neys") || on_laptop || on_jmac;
    on_hpcc_system = contains(hostname,"compute");
    
%..................Parameters that depend on station/gridded and on which system we're running on

    if (any(contains('station',lower(runType))))
        if (on_dev_system)
            DP = data_params_station_dev(on_laptop, on_jmac);
        elseif (on_hpcc_system)            
            DP = data_params_station_hpcc();
        else
            DP = data_params_station_KMQ();
        end       
    else
        if (on_dev_system)
            DP = data_params_gridded_dev(on_laptop, on_jmac);
        elseif (on_hpcc_system)            
            DP = data_params_gridded_hpcc();
        else
            DP = data_params_gridded_KMQ();
        end
    end
    
end


function [DP] =  data_params_station_dev(on_laptop, on_jmac)

    
%       Data parameters
%
%           this is the default. 
    DP.region = "ncamerica";     % default

%           Data directories
    
    if (on_jmac)
        DP.datadir      = "/Volumes/jcsf_data/data/obs/";
        DP.outdir       = "/Volumes/jcsf_data/data/downscaled_stations/REGION";        
    elseif (on_laptop)
        DP.datadir      = "/Volumes/2018_1/data";
        DP.outdir       = ".";
    else
        DP.datadir      = "/Volumes/lacie_1/data/obs/";       
        DP.outdir       = "/Volumes/lacie_1/data/downscaled/ARRM_V2/station/REGION";        
    end    
%           Filenames 
    DP.outname   = "downscaled.stations.MODEL.ENSEMBLE.VARNAME.SCENARIO.REGION.MDLSTARTYEAR.MDLENDYEAR.LLGRID.nc";

end

%----------------------------------------------------

function [DP] =  data_params_gridded_dev(on_laptop, on_jmac)

%          CF Conventions Global Attributes. 
%          Modify as appropriate.  You can override these in 

%       Data parameters
%
%           this is the default. 
    DP.region = "ncamerica";     % default

%           Data directories
    
    if (on_jmac)
        DP.datadir      = "/Volumes/jcsf_data/data/gcm_cmip5/daily.1900.2100/rotated";
        DP.outdir       = "/Volumes/jcsf_data/data/downscaled/arrm_v2";    
    elseif (on_laptop)
        DP.datadir      = "/Volumes/2018_1/data";
        DP.outdir       = ".";
    else
        DP.datadir      = "/Volumes/lacie_1/data/gcm_cmip5/daily.1900.2100/rotated";
        DP.outdir       = "/Volumes/lacie_1/data/downscaled/arrm_v2";            
    end    
    
%           Filenames 
    DP.fnames       = "";  % this is just hardcoded to use CCSM4 as my "obs" data for testing.

    % output file & folder 
%   Note:  outdir must exist and be writeable...code does not create the folder.    
    %DP.outname   = "downscaled.MODEL.ENSEMBLE.VARNAME.SCENARIO.MDLSTARTYEAR.MDLENDYEAR.mini.RUNID.nc";
    DP.outname   = "downscaled.MODEL.ENSEMBLE.VARNAME.SCENARIO.REGION.MDLSTARTYEAR.MDLENDYEAR.LLGRID.nc";

   
end

%----------------------------------------------------

function [DP] =  data_params_station_hpcc()

%          CF Conventions Global Attributes. 
%          Modify as appropriate.  You can override these in 
    
%       Data parameters
%           this is the default. 
    DP.region = "ncamerica";     % this needs to be passed in.
%           Data directories
    
    DP.datadir      = "/lustre/scratch/iscottfl/obs/stations_netcdf";

%           Filenames 
    DP.fnames     = ["MODEL.ENSEMBLE.VARNAME.hist.day.1900.2005.llt.nc","MODEL.ENSEMBLE.VARNAME.SCENARIO.day.2006.2100.llt.nc"];
  
    DP.outdir       = "/lustre/scratch/iscottfl/downscaled_stations/ncamerica/REGION";    
    DP.outname      = "downscaled.stations.MODEL.ENSEMBLE.VARNAME.SCENARIO.REGION.MDLSTARTYEAR.MDLENDYEAR.LLGRID.nc";

end

%----------------------------------------------------

function [DP] =  data_params_gridded_hpcc()

%          CF Conventions Global Attributes. 
%          Modify as appropriate.  You can override these in 

%       Data parameters
%
%           Data directories
    
    DP.datadir      = "/lustre/scratch/iscottfl/cmip5_rotated";
    
%           Filenames 
    DP.fnames     = ["MODEL.ENSEMBLE.VARNAME.hist.day.1900.2005.llt.nc", "MODEL.ENSEMBLE.VARNAME.SCENARIO.day.2006.2100.llt.nc"];


    % data directories  
%   Note:  out_dir must exist and be writeable...code does not create the folder.
    DP.outdir       = "/lustre/scratch/iscottfl/downscaled_gridded/";   
%   DP.outname      = "downscaled.MODEL.ENSEMBLE.VARNAME.SCENARIO.MDLSTARTYEAR.MDLENDYEAR.mini.RUNID.nc";
    DP.outname      = "downscaled.MODEL.ENSEMBLE.VARNAME.SCENARIO.REGION.MDLSTARTYEAR.MDLENDYEAR.LLGRID.nc";

end

%----------------------------------------------------

function [DP] =  data_params_station_KMQ()


%       Data parameters
%           Data directories
    
    DP.datadir      = "/data/gcm_cmip5/daily.1900.2100/VARNAME/";
    
%           Filenames 
    DP.fnames     = ["hist/MODEL.ENSEMBLE.VARNAME.hist.day.1900.2005.llt.nc", "SCENARIO/MODEL.ENSEMBLE.VARNAME.SCENARIO.day.2006.2100.llt.nc"];

    % data directories  
%   Note:  outdir must exist and be writeable...code does not create the folder.
    DP.outdir       = "/data/downscaled_stations/REGION";
    DP.outname      = "downscaled.stations.MODEL.ENSEMBLE.VARNAME.SCENARIO.REGION.MDLSTARTYEAR.MDLENDYEAR.LLGRID.nc";

end
   
%----------------------------------------------------

function [DP] =  data_params_gridded_KMQ()

%          CF Conventions Global Attributes. 
%          Modify as appropriate.  You can override these in 

%       Data parameters
%
%           Data directories
    
    DP.datadir       = "/data/gcm_cmip5/daily.1900.2100/VARNAME/";

    
%           Filenames 
    DP.fnames     = ["hist/MODEL.ENSEMBLE.VARNAME.hist.day.1900.2005.llt.nc", "SCENARIO/MODEL.ENSEMBLE.VARNAME.SCENARIO.day.2006.2100.llt.nc"];

    % out  
%   Note:  outdir must exist and be writeable...code does not create the folder.
    DP.outdir       = "/data/downscaled/arrmv2/gridded/VARNAME/MODEL/SCENARIO";
%   DP.outname      = "downscaled.MODEL.ENSEMBLE.VARNAME.SCENARIO.MDLSTARTYEAR.MDLENDYEAR.mini.RUNID.nc";
    DP.outname      = "downscaled.MODEL.ENSEMBLE.VARNAME.SCENARIO.REGION.MDLSTARTYEAR.MDLENDYEAR.LLGRID.nc";

   
end

%----------------------------------------------------


