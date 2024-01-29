classdef ARRM_V2_DataParams
    %ARRM_V2 Data parameter container and parsing code.
%   Constructor creates, and update(...) modifies a DataParams object for use in ARRM_V2 programs
%
%   Both the constructor and update(...) take matlab name, value pairs to set any of the dataParams.
%
%       data members can be set directly, but using the parser will error-check the values, and in some cases update
%       related parameters.
%   
%       Use variable name in quotes, followed by desired value, as keyword/value pair.
%       Where value requires multiple values (only trend_yr_flags), provide vector of values.
%       Explanations and default values are given below.
%       Terms may be in any order;  any term not specified uses the defaults given below.
%       (varargin may also be a vararginStruct from matlab's parse(...) function.)
%
%   6/7/22 icsf Modified to make sure lats & lons are of type double.  Some models use single, which can cause problems with location precision.
%
%               Examples:
%           'varName', 'tasmax', ', ...
%
%       varargin Keywords       Description
%       _________________       ___________

    properties
                                                % can be set to outside range of actual data.  Code will work with the
                                                % available years.
        trend_yrs       = [];               % data years to use to calculate long term trend.  set to [] to skip calculating trend
        base_yrs        = [];               % data years to use for base (historical/observation) period.
        rolling_yrs     = [];               % year range [start,end] to use for rolling pdfs.  Set to [] for OBS data
        rolling_steps   = [];               % years where rolling steps should be calculated, relative to rolling_yrs(1).  
                                            %   0-based... steps are calculated at rolling_yrs(1) + rolling_steps.
                                            %   Normally calculated by calling obj.calc_rolling_steps(yrstep).
        data_yrs        = [];               % data years to pull (for disaggregation only).  If empty, will be calculated by  code 
                                            %  as total range of bas_yrs and rolling_yrs, extended by 1/2*pdf_yrs at both ends.  
        data_final_yrs  = [];               % data years to retain in results. if empty, is set to total range of base_yrs and rolling_yrs
        
            % data Specifications
            % filenames and data directories.
            %  provide either fnames & datadir or obs,mdl & hist. filenames and dirs.
                % for basic signal decomposition of one datastream
        fnames          = strings(0,0);
        default_fnames  = strings(0,0);     % provide this in localized params if desired.
        datadir         = strings(0,0); 
        outname         = strings(0,0);     % currently not used in regular DataParams.  Used in DownScalingParams subclass.
        outdir          = strings(0,0);     % currently not used in regular DataParams.  Used in DownScalingParams subclass.
        interp_method   = 'bilinear';
        
        yrlen           = 365;  % length of year for input and output data.  Will remove or pad leap days, etc. as needed
        file_lats       = [];   % These are for the lats & lons extracted from the netcdf metadata.
        file_lons       = [];   % leave blank unless you need to override the lats & lons in the file!!!!!
        lats            = [];   % these are for the lat & lon locations to calculate data for
        lons            = [];   % Set these in your input parameter.  NOTE:  lons should be in range 0-360, not -180:180.  
        stninfo         = strings(0,0);   % stninfo table, for stations run
        model           = strings(0,0);
        varname         = strings(0,0);   % varname for model data
        varlongname     = strings(0,0);
        obsvname        = strings(0,0);   % (needed for DownscalingParams, but needs to be here for parameter replacing in strings.)
        units           = strings(0,0);   % cannot be set here, but can be set after reading units from netcdf file.
        ensemble        = strings(0,0);   % if blank, takes from filename
        scenario        = strings(0,0);   % if blank, takes from filename
        scaling         = 'linear';
        cdf_append_pts  = [-3,3; -4.75, 4.75];  % point(s) at which to replace local CDF with basic_cdf_mdl. Should be 2x2.  
                                                % 1st row defines where empirical cdf gets attached;  2nd row where
                                                % CDF should start tapering exponentially to zero.  [-3, 3;  -4.75,4.75] is
                                                % a good set of values to use.  -4.4,4.4 are ~5e-6. Use 4.75 for 1e-6 
                                                %  Emprical CDFs seem to be good out to just above 1e-6.
        cdf_append_type = "normal";         % 'emipirical','normal' or 'sk_normal'
        
        prcp_min        = 0.1;         % min threshold for precip;  below this we call it a dry day ("trace")
                                        % this will be adjusted when trimming precip to match obs precip counts.
        prcp_median     = [];           % will be set later to daily median precip value so we can differentiate NAs for drizzle vs NAs for heavy precip.
%       prcp_scaling    = 1.0;          % actual scaling info is now stored in a struct in the DA object.
%                                       % that struct has info for reversing the scaling and edges for binning.
        
        runType         = 'normal';        
        monthly_valid_count = 600;      % require 20 years' valid data for OBS (and HIST, for DSPs) data.
        Unmatched       = [];
        
        runID           = strings(0,0);     % single-word string identifying run:  model_ensemble_varname_scenario_gridix_ngrids

        runLbl          = strings(0,0);     % string with model, varname, scenario, ensemble ( & stnID & stnName for single-points)
        stnID           = strings(0,0);     % station ID, or gridpoint x,y indices within gridbox.
        stnName         = strings(0,0);     % station name, or gridpoint (lat,lon)
        region          = 'ncamerica';
        llgrid_size     = [15,15];      % size of [lat,lon] grid for creating output minifiles.  Usually set to something like [15,15];
                                            % if not empty, ARRM_V2_wrapper will create one child process for each
                                            % subgrid, and process only the stations or gridded obs locations within the
                                            % llgrid area.  Starts from floor(llgrid_size) for entire area being
                                            % processed.  NOTE:  if not set to empty, then child process reads all model
                                            % date for the entire grid.  So don't set too large if memory is limited.
        llgrid_lbl      = "";           % string to identify lat/lon region of minifile.  Should begin w/ a '.' if used
                                            % as part of the output filename.
        
        ARRM_V2_version = "";           % ARRM_V2_version_info(true);       % this doesn't actually call the function;  it seems to store the info available when the file was last modified.
                                                                            % so be sure to call ARRM_V2_version_info(true) in your code somewhere.  
 
        ARRM_V2_dir     = strings(0,0);

        FillValue       = single(1e20);
        to_struct       = false;
        keep_all        = false;            % retain extended data for plotting
        extended        = false;            % save extended output.  Currently just probability values for each point, stored as 2-byte packed zvalues.
        zval_offset     = -6;               % scaling and offset used to store z-vals of probabilities.
        zval_scaling    = 12/65534;         % (earlier typo:  12/64434 !
                                            % This gives us a resolution of 

                                            % things set by main routine, like where to write log info, etc.
        do_log          = 1;            % controls where logfile info is written. default:  1
                                        %   0: No log.  write only error messages, and those will go to console.  
                                        %   1: write  log info to log file (overwrite), then copy to console at end 
                                        %   2: append log info to log file, then copy to console at end 
                                        %       generally, for production, or if using parfor, use 1.  
   
        logname         = strings(0,0);
        do_parfor       = false;        % flags whether running in parallel or serial mode.
                                        % if in serial mode, then calling DP.print_log(...) instead of fprintf will
                                        % write to console and to file.
                                        % calling DP.warn_log(...) will always write to stderr and to log file.
        fidlog          = 1;    % fid(s) for log file.  by default, code will write to console, but connect to a file to log output.
                                    % if fidlog is vector (as in [1, fid], then calling DP.print_log(fmt, vars) will write to multiple locations
        progname        = [];   % Can be filled in by initiating program to identify which main program is running.
        exit_on_error   = true;
        
                                % info for drawing figures, if desired.
        displayResults  = false;
        figbase         = 1000;
        figflags        = [];   % 6 boolean flags;  enables drawing of 6 different figures.  See jc_plot_data_decomposition.m
        figname         = [];   % no extension!  Will append fig_ID, and then extension.  
        figext          = [];   % can be string, or cell array {'mat','tif','jpg'...}.  Will save one of each.
                                %   on return, will be string array, each starting with a period.
        figpos          = [];   % [low-left-x, low-left-y,nx,ny].  can be 1 for each figure, or 1 for all figs.
                                    % must be 4 columns wide.  low_left_x, low_left_y, nx, nx].  Can 1, or 1 for each
                                    % figure identified in figflags.
                                    
        debug_flag      = [];           % settable flag for debugging.
        internal_plotflag   = int64(0); % controls internal plotting for debugging. See obj.intern_plotflag(...) and obj.set_internal_plotflag(...)
                                        % bit 0 is global on/off, and can be set by calling obj.set_internal_plotflag(true/false)
                                        % bits 1-n can be set or tested individually.  
                                        % obj.intern_plotflag(k) returns false if bit 0 is off, or returns
                                        %                         bit k's state if bit 0 is on.
                                        % See code in ARRM_V2_DisaggregateSignal.m for info on plots generated
                                        % (plots may change in the future as needed...)
                                        %
                                        %   example:  obj.set_internal_plotflag([0,2,5,7]);  or pass in
                                        %   "do_internal_plots",[0,1,2,5,7], ...
                                        %                   Ian:  you need to document what the different plot flags do!  you keep forgetting! 
                                        %       0       turns on global, to enable plots
                                        %       1,2,5,7 turns on plots 1,2,5 & 7.
                                        %                   at present, these turn on plots of generating pdfs for 
                                        %                   basic_yrs (1), or rolling_yrs  (1,2,...)
                                        %       [0,32] turns on plots for appending to pdfs with sknormal or empirical or normal.  
                                        %
                                        %  if figbase not set, figbase will be set to 1.
        
                                    
                                 % list of fields to exclude when saving Params to netCDF file
       ncExcludeList     = {'default_fnames', 'file_lats', 'file_lons', 'lats', 'lons', ...
                             'stninfo', 'runID', 'stnID', 'stnName', 'runLbl', 'to_struct', 'keep_all','extended', ...
                             'fidlog', 'displayResults','figbase', 'figflags',...
                             'figname', 'figext', 'figpos', 'internal_plotflag', 'Unmatched','ncExcludeList'};
    end
    
    properties(Dependent)
        isPrecipRun;
        isStationRun;
        isDownscalingRun;                            
        
    end
        
    methods
        function obj = ARRM_V2_DataParams(varargin)

                    % check for runType.  We need to pass it to the localized DataParams function.
%             ix=find(strcmp(varargin,"runType"),1,'last');         % comparison against all elements at once fails
            ix=[];
            for i=1:2:length(varargin); if (strcmpi('runType',varargin{i})), ix=i; end; end
            if (~isempty(ix))
                runType=string(varargin{ix+1});
%                 varargin(ix:ix+1) = [];         % don't want to remove this, because it is used in both RunParams and DataParams
            else
                runType="unknown"; 
            end
                    
                % apply any localized parameters if the function ARRM_V2_DataParams_localized( ) exists.
            if (exist('ARRM_V2_DataParams_TTUCSC.m','file'))
                obj = obj.DP_update(ARRM_V2_DataParams_TTUCSC(runType));
                if (~isempty(fieldnames(obj.Unmatched))), error('error:  unrecognized parameters from DataParams_TTUCSC'); end
            end
            if (exist('ARRM_V2_DataParams_localized.m','file'))
                obj = obj.DP_update(ARRM_V2_DataParams_localized(runType));
                if (~isempty(fieldnames(obj.Unmatched))), error('error:  unrecognized parameters from DataParams_localized'); end
            end
            
%                 % now, check for a DP in the varagins
%                 % This applies any DP arguments before we run update(...) on the remaining args, so we can override
%                 % any DP arguments.
%             obj = obj.parse_for_DP(args);                        
            
                % now apply any command-line params
            obj = DP_update(obj, varargin{:});
            
%           obj = obj.replace_keywords();            % now done in finalize(...)
            
        end
        
        function obj = update(obj, varargin)
            
            % Updates DataParams object.  Uses current values as defaults, and parses any key/value pairs
            % Call this when you want to set multiple parameters with a single function call.
            % To update a single parameter, just set the parameter directly.
                
            obj = obj.DP_update(varargin{:});
        end

        function obj = finalize(obj, varargin)
            % does final check on obj, filling in defaults, etc. as needed.
            % Can be passed additional parameters, & will update first, then finalize.
            
            if (exist('varargin','var') && ~isempty(varargin))
                obj = obj.update(varargin{:});
            end
            
            if (~isempty_s(obj.fnames))
                obj = obj.parse_fnames(obj.fnames, false);      % fill in model, varname, etc. from filenames if needed.
            end
            
            if ( ~obj.isDownscalingRun)
                obj = obj.make_fnames();
                
                if (obj.isStationRun)      % for disaggregation run, get varname from station netcdf file 
                            % possibly could assemble station filename here if we were smart, Ian...
                    if (ischar_s(obj.stninfo))
                        obj.stninfo = QC_get_site_table(obj.stninfo, "stnID", obj.stnID);
                    end
                    if (~isQCstntbl(obj.stninfo)), error("error:  stninfo table not read in"); end
                            % extract station info for stations of interest if not processing all stations in the site table.
                    if (length(string(obj.stnID)) ~=size(obj.stninfo,1))
                        obj.stninfo = QC_get_site_table(obj.stninfo, "stnID", obj.stnID);
                    end

                    if (isempty(obj.fnames))
                        obj.fnames = obj.stninfo.UserData.ncName;
                    end
                    if (isempty(obj.varname) && size(obj.stninfo.Properties.UserData.varName,1)==1)
                        obj.varname = obj.stninfo.Properties.UserData.varName;
                    end                    
                end
                obj.fnames=string(obj.fnames);            
                obj = obj.update_filenames('datadir','fnames');           % make sure we have fully absolute filenames, not relative ones.
                obj = obj.update_filenames('outdir','outname');
            end
            
                % and make sure lats & lons are of type double
            obj.lats = double(obj.lats);
            obj.lons = double(obj.lons);
            obj.file_lats = double(obj.file_lats);
            obj.file_lons = double(obj.file_lons);
            
            if (istable(obj.stninfo))
                obj.stninfo.lat = double(obj.stninfo.lat);
                obj.stninfo.lon = double(obj.stninfo.lon);
            end            
            
            obj.lons = mod(obj.lons,360);   % make sure lons are in range 0-360, not -180 to 180.            
            
            obj = obj.replace_keywords();
            
        end
        
        function yrs = rolling_step_years(obj, n, pdf_yrstep)
            % returns year range ([start_yr, end_yr]) for rolling step n.            
            yrs = obj.rolling_yrs;
            try
                yrs(1) = max(obj.data_yrs(1), yrs(1) + obj.rolling_steps(n) - floor(pdf_yrstep/2));
            catch
                oops();
            end
            if (n ~=length(obj.rolling_steps))
                yrs(2) = yrs(1) + pdf_yrstep-1;
            end
        end
        
       function [start_ix, end_ix, nyrs] = using_range(obj, yr_type, rolling_set, pdf_yrstep)
           % returns start and end indices of data to use for the specified year range, and number of years.  
           % yr_type can be string 'data_yrs','base_yrs','rolling_yrs', or [start_yr, end_yr];
            if (exist('rolling_set','var') && ~isempty(rolling_set))
                yrs = obj.rolling_step_years(rolling_set, pdf_yrstep);                
            elseif (isnumeric(yr_type))
                yrs = yr_type;
            else
                yrs = obj.(yr_type);
            end
            nyrs = yrs(2)-yrs(1)+1;
            start_ix = (yrs(1) - obj.data_yrs(1)  )*obj.yrlen+1;
            end_ix   = (yrs(2) - obj.data_yrs(1)+1)*obj.yrlen;
        end
                
                
        function obj =set_yr_limits(obj, file_yrs)
            % function to limit the various year ranges to the data available in the files.
            % file_yrs should be set to the range of data available in the file
            % rolling_pdf_yrs is the # of years for each pdf in the rolling pdfs;  s/b RP.pdf_yrs.
            
            if(isempty(obj.base_yrs) && isempty(obj.rolling_yrs)), error("error:  DSP::set_yr_limits():  neither base_yrs nor rolling_yrs are set"); end

                    % limit base_yrs to data available in files
            if (~isempty(obj.base_yrs))
                obj.base_yrs = min_yr_range(obj.base_yrs, file_yrs);
            end
                    % limit rolling yrs to data available in files
                    
                    % limit the trend years
            if (~isempty(obj.trend_yrs))
                obj.trend_yrs = min_yr_range(obj.trend_yrs, file_yrs); 
            end
            
                    % limit the rolling years
            if (~isempty(obj.rolling_yrs))
                obj.rolling_yrs = min_yr_range(obj.rolling_yrs, file_yrs); 
            end
            
                    % limit the data_final_yrs
            if (~isempty(obj.data_final_yrs))
                obj.data_final_yrs = min_yr_range(obj.data_final_yrs, file_yrs);
            end
            if (~isempty(obj.data_yrs))                    
                obj.data_yrs = min_yr_range(obj.data_yrs, file_yrs);
            end
            obj.data_yrs  = max_yr_range(obj.data_yrs, obj.trend_yrs, obj.base_yrs, obj.rolling_yrs, obj.data_final_yrs);
        end
        
        function out = replace_keywords(obj, kwds, instring)
            %   replaces keywords either in instring or in all the properties of obj.
            %   if instring is provided, then all keywords in instring are replaced, and new string returned.
            %   Otherwise, all instances of the keywords are replaced in all the objects' properties.
            %
            %   NOTE:  if replacing LAT1 & LATEND, or LON1 and LONEND, be sure to specify LATEND or LONEND before LAT or LON 
            %   
            if (~exist('kwds','var') || isempty_s(kwds))
                kwds=["MODEL","ENSEMBLE","VARNAME","OBSVNAME","SCENARIO","STNID","RUNID","REGION",...
                      "MDLSTARTYEAR","MDLENDYEAR","LATEND","LONEND","LAT","LON","LAT1","LON1","LLGRID", "CURDATE"];
            else
                kwds = string(kwds);    % make sure keywords are strings, not cell array of chars...
            end
            if (exist('instring','var'))
                out = instring;
                if (ischar(out))
                    was_char = true;
                    out = string(instring);
                else
                    was_char = false;
                end
                for k=1:length(out)
                   for i=1:length(kwds)
                        kwd=kwds(i);
                        if (iscell(out))
                            out{k} = obj.do_replace(out{k},kwd);
                        else
                            out(k) = obj.do_replace(out(k),kwd);
                        end
                   end
                end
                if (was_char), out = char(out); end     % if it was a single char array, convert it back to chars.
            else
                for i=1:length(kwds)
                    kwd=kwds(i);
                    props=properties(obj);
                    for j=1:length(props)
                        prop=props{j};
                        if (strncmp(prop,'default_',8)), continue; end  % skip replacing keywords in any default parameter.
                        if (isempty(obj.(prop))), continue; end
                        if (isstruct(obj.(prop)))
                            flds = fieldnames(obj.(prop));
                            for k=1:length(flds)
                                fld = flds{k};
                                if (ischar_s(obj.(prop).(fld)))
                                    obj.(prop).(fld) = obj.do_replace(obj.(prop).(fld),kwd);
                                end
                            end
                        elseif (ischar_s(obj.(prop)))
                            try
                                obj.(prop) = obj.do_replace(obj.(prop),kwd);
                            catch
                                oops();
                            end
                        elseif (isstring(obj.(prop)) && ~strcmp(prop,'stnID'))
                            for k=1:length(obj.(prop))
                                obj.(prop)(k) = obj.do_replace(obj.(prop)(k),kwd);
                            end
                        elseif (iscell(obj.(prop)))
                            for k=1:length(obj.(prop))
                                if (ischar_s(obj.(prop){k}))
                                    obj.(prop){k} = obj.do_replace(obj.(prop){k},kwd);
                                end
                            end                           
                        end
                    end
                end
                out=obj;
            end            
        end
        
        function yesno = get.isPrecipRun(obj)
            yesno = any(contains(["pr","prcp","precip","prec"], lower(obj.varname))) || any(contains(["pr","prcp","precip","prec","precipitation"],lower(obj.runType)));
        end
                
        function yesno = get.isDownscalingRun(obj)
            yesno = isprop(obj, "obsnames");
        end
                
        function yesno = get.isStationRun(obj)
            yesno = ~isempty(obj.stninfo);
        end
                
        function DP_struct = toStruct(obj)
            props = properties(obj);
            DP_struct = struct();
            for i=1:length(props)
                prop = props{i};
                DP_struct.(prop) = obj.(prop);
            end
        end
            
        function roll_steps = calc_rolling_steps(obj, yrstep)
            if (isempty(obj.rolling_yrs) || isempty(yrstep) || yrstep<=0)
                roll_steps = [];
            else
                nyrs = obj.rolling_yrs(2)-obj.rolling_yrs(1)+1;
                roll_steps = floor(yrstep/2):yrstep:nyrs;
            end
        end
        
        function [ok, badnames] = check_files_exist(obj, props)
            if (~exist('props','var'))
                props="fnames";
            end
            
            ok = true;
            badnames="";
            nbad = 0;
            for j=1:length(props)
                prop=string(props(j));
                myfnames = string(obj.(prop));
                for i=1:length(myfnames)
                    if (~exist(myfnames(i),'file'))
                        ok = false;
                        if (nargout == 0)
                            error("error:  %s file %s does not exist", prop, myfnames(i));
                        elseif (nargout == 1)
                            return
                        end
                        nbad=nbad+1;
                        badnames(nbad) = myfnames(i);
                     end
                end
            end
        end    
        
        function obj = parse_fnames(obj, fnames, do_force)
                % looks at filenames, extracts model, ensemble, variable, etc. and updates the DP fields.
            if (exist('fnames','var') && ~isempty_s(fnames))
                obj.fnames = string(fnames);
            end
            if (nargin < 3), do_force = false; end
            run_info = ARRM_V2_parse_netcdf_filenames(obj.fnames);
            fields=fieldnames(run_info);
            for i=1:length(fields)
                f=fields{i};
                if (isempty(obj.(f)) || do_force)
                    obj.(f) = run_info.(f);
                end
            end
        end
        
        function obj = make_fnames(obj, prop)
                % replaces keywords in fnames, using the default fnames field if fnames is empty.
                % This function is called by finalize().  If you don't call finalize() yourself, then
                % call this function later if you want to assemble the filenames on the fly from the
                % keywords (model, varname, ensemble, scenario, etc.).
                % To do this, supply a default filename in the ARRM_V2_DataParams_localized(...) function
                % You'll probably want to do this separately for station runs and gridded runs.
            if (~exist('prop','var') || isempty(prop))
                prop = 'fnames';
            end
            if (isempty_s(obj.(prop)))
                defprop = sprintf('default_%s', prop);
                obj.(prop) = obj.(defprop);
            end
            if (isempty_s(obj.(prop))), return; end   % nothing to do!  bail out.
            
            obj.(prop) = obj.replace_keywords([],obj.(prop));
        end
                      % These next 2 functions are used for generating debugging plots.
                      % These plots can be used to generate figures to help explain the disaggregation process.
                                        % See code in ARRM_V2_DisaggregateSignal.m for info on plots generated
                                        % (plots may change in the future as needed...)
                      % by default, only 1 64-bit flag (flagnum==1 or missing).  But if set_internal_plotflag(...) can 
                      % append additional 64-bit flags for other purposes.  See notes above for info on plots generated
          function yn = intern_plotflag(obj, k, flagnum)
            if (nargin < 3)
                flag = obj.internal_plotflag(1);
            else
                flag = obj.internal_plotflag(flagnum);
            end
                
            if (nargin == 1 || k == 0) 
                yn = bitget(flag,1);
            else
                yn = bitget(flag,1) & bitget(flag,k+1);
            end
        end
        
        function [obj, oldflag] = set_internal_plotflag(obj, k, newstate, flagnum)
            % if k is true/flase, turns on or off global flag to enable or disable internal plots
            % if k is numeric (1-63), turns on or off (newstate==false) flags 1-63.
            % k can be an array of flags.  Be sure to set flag 0 to enable plots if turning on an internal plot.
            if (nargin < 4), flagnum = 1; end
            if (islogical(k))
                oldflag = bitget(obj.internal_plotflag(flagnum),1);
                obj.internal_plotflag(flagnum) = bitset(obj.internal_plotflag(flagnum),1,k);
            else
                if (nargin < 3), newstate = true; end
                if (iscell(k))
                    oldflag = zeros(length(k),1);
                    for i=1:length(k)
                        [obj, oldflag(k)] = set_internal_plotflag(obj,k{i},newstate);
                    end
                else
                    nflags = numel(k);                    
                    for kk=1:nflags
                        obj.internal_plotflag(flagnum) = bitset(obj.internal_plotflag(flagnum),k(kk)+1, newstate);
                    end
                end
            end
            if (isempty(obj.figbase)), obj.figbase = 1; end
        end
        
        function print_log(obj, fmt, varargin)
                % if not running parallel, does fprintf of varargins to console and to file
                % else just does fprintf to file.
            if (~obj.do_log), return; end
            fids = obj.fidlog;
            if (~obj.do_parfor)
                if (~any(fids==1)), fids(end+1)=1; end
            end
            
%             fprintf("print_log:  do_parfor:  %d nfids: %d: ", obj.do_parfor, length(fids));
%             for i=1:length(fids)
%                 fprintf("fid %3d: %3d\n", i, fids(i));
%             end
%             dbstack(1); 
%             
%             if (contains(fmt,"xxx")) fprintf("\n***** %d xxx \n", i); end
%             if (contains(fmt,"yyy")) fprintf("\n***** %d yyy \n", i); end
%             if (contains(fmt,"zzz")) fprintf("\n***** %d zzz \n", i); end
%             fprintf("logging:  %s\n", fmt);            
% 
%             fprintf("^^^\n");
%             pause(.025);
            for i=1:length(fids)
                fprintf(fids(i), fmt, varargin{:});
            end
        end
        function warn_log(obj, fmt, varargin)
                % does fprintf of varargins to stderr and to log file.
                % (skips to console if console is an element of fid_log)
            fids = obj.fidlog;
            fids(end+1)=2;
            for i=1:length(fids)
                if (fids(i)==1), continue; end          % skip console, because we've added stderr.
                fprintf(fids(i), fmt, varargin{:});
            end
        end
        function error_log(obj, fmt, varargin)
                % does fprintf of varargins to stderr and to log file.
                % (skips console if console is an element of fid_log)
                % then throws the same message as an error.
                
            st=dbstack();
            if (length(st)>1)
                mycaller=st(2).name;
                myline=st(2).line;
            end
            fids = obj.fidlog;
            fids(end+1)=2;
            for i=1:length(fids)
                if (fids(i)==1), continue; end          % skip console, because we've added stderr.
                if (exist("mycaller","var")), fprintf(fids(i), "%s:%d : ", mycaller, myline); end
                fprintf(fids(i), fmt, varargin{:});
            end
            error(fmt, varargin{:});
        end
    end
        
    methods(Access = protected)
        % prepends directory to filenames if name provided is not an absolute path.
        %   does not parse the filename for model, ensemble, etc. or update those fields.
        %   Use obj.parse_fnames(...) for that.
        
        function obj = update_filenames(obj, dirprop, namesprop)
            if (isempty(obj.(dirprop)) || isempty(obj.(namesprop))), return; end
            if (~isAbsolute(obj.(dirprop)))
                obj.(dirprop)=fullfile(pwd(),obj.(dirprop));
            end                
            if (isstring(obj.(namesprop)))
                for i=1:length(obj.(namesprop))
                    if (~isAbsolute(obj.(namesprop)(i)))
                        obj.(namesprop)(i) = fullfile(obj.(dirprop),obj.(namesprop)(i));
                    end
                end
            elseif (iscell(obj.(namesprop)))
                for i=1:length(obj.(namesprop))
                    if (~isAbsolute(obj.(namesprop){i}))
                        obj.(namesprop){i} = fullfile(obj.(dirprop),obj.(namesprop){i});
                    end
                end
            elseif (~isAbsolute(obj.(namesprop)))
                obj.(namesprop) = fullfile(obj.(dirprop),obj.(namesprop));
            end
        end
        
        function outstring = do_replace(obj, instring, kwd)
            % replaces all occurrences of kwd with appropriate value
            % This is for replacing keywords in various things, like filenames, netcdf metadata, etc.

            outstring = instring;

            if (~ischar_s(instring) || ~contains(instring,kwd))
                return;
            end

            if (ischar(instring))
                waschar = true; 
            else
                waschar = false; 
            end
            
            if     (kwd == "MODEL" && ischar_s(obj.model) && ~isempty_s(obj.model))
                outstring=strrep(instring,kwd, obj.model);
            elseif (kwd == "ENSEMBLE" && ischar_s(obj.ensemble) && ~isempty_s(obj.ensemble))
                outstring=strrep(instring,kwd, obj.ensemble);
            elseif (kwd == "OBSVNAME" && ischar_s(obj.obsvname) && ~isempty_s(obj.obsvname))
                outstring=strrep(instring,kwd, obj.obsvname);
            elseif (kwd == "VARNAME" && ischar_s(obj.varname) && ~isempty_s(obj.varname))
                outstring=strrep(instring,kwd, obj.varname);
            elseif (kwd == "SCENARIO" && ischar_s(obj.scenario) && ~isempty_s(obj.scenario))
                outstring=strrep(instring,kwd, obj.scenario);
            elseif (kwd == "STNID" && ischar_s(obj.stnID) && ~isempty_s(obj.stnID))
                outstring=strrep(instring,kwd, obj.stnID);
            elseif (kwd == "RUNID" && ischar_s(obj.runID) && ~isempty_s(obj.runID))
                outstring=strrep(instring,kwd, obj.runID);
%             elseif (kwd == "RUNLBL" && ischar_s(obj.runLbl) && ~isempty_s(obj.runLbl))
%                 outstring=strrep(instring,kwd, obj.runlbl);
            elseif (kwd == "STNNAME" && ischar_s(obj.stnName) && ~isempty_s(obj.stnName))
                outstring=strrep(instring,kwd, obj.stnName);
            elseif (kwd == "LLGRID" && ischar_s(obj.llgrid_lbl) && ~isempty_s(obj.llgrid_lbl)) 
                outstring=strrep(instring,kwd, obj.llgrid_lbl);    %note underscore.  If used for filename, make sure llgrid_lbl includes period at start.
            elseif (kwd == "REGION" && ischar_s(obj.region) && ~isempty_s(obj.region))
                outstring=strrep(instring,kwd, obj.region);
            elseif (kwd == "MDLSTARTYEAR" && ~isempty(obj.rolling_yrs))
                outstring=strrep(instring,kwd, sprintf("%04d", obj.rolling_yrs(1)));
            elseif (kwd == "MDLENDYEAR" &&   ~isempty(obj.rolling_yrs))
                outstring=strrep(instring,kwd, sprintf("%04d", obj.rolling_yrs(2)));
            elseif ((kwd== "LAT" || kwd=="LAT1") && ~isempty(obj.lats))
                outstring=strrep(instring,kwd, sprintf("%.4f", obj.lats(1)));
            elseif ((kwd== "LON" || kwd=="LON1") && ~isempty(obj.lons))
                outstring=strrep(instring,kwd, sprintf("%.4f", obj.lons(1)));
            elseif (kwd== "LATEND" && ~isempty(obj.lats))
                outstring=strrep(instring,kwd, sprintf("%.4f", obj.lats(end)));
            elseif (kwd== "LONEND" && ~isempty(obj.lons))
                outstring=strrep(instring,kwd, sprintf("%.4f", obj.lons(end)));
            elseif (kwd=="CURDATE")
                outstring=datestr(now,'yyyy-mm-dd HH:MM:SS');
            end
            if (waschar), outstring = char(outstring); end
        end    
    end

    methods(Access=private)
        
        function obj = DP_update(obj, varargin)
            % I need to add verification to these still, Ian!
            

                    % first, check for a DP in the varargins.
                        % parse input for DP
            ix=find(strcmp(varargin,"DP"));
            if (~isempty(ix))
                for j=1:length(varargin)
                    DP=varargin{ix{j}+1}; 
                    if (~isa(DP,"ARRM_V2_DataParams")), error("DP argument is not an ARRM_V2_DataParams object"); end
                    props = properties(DP);
                    for i=1:length(props)
                        try
                            obj.(props{i}) = DP.(props{i});     % in a try block so loop doesn't abort on dependent properties.
                        catch
                        end
                    end
                end            
                for i=length(ix):-1:1
                    varargin{ix:ix+1} = [];
                end
            end
                % and parse the remaining arguments.
            p = inputParser;
            p.StructExpand = true;
            p.KeepUnmatched = true;
    
            addParameter(p,'trend_yrs',obj.trend_yrs);
            addParameter(p,'base_yrs',obj.base_yrs);
            addParameter(p,'rolling_yrs',obj.rolling_yrs);
            addParameter(p,'data_yrs',obj.data_yrs);
            addParameter(p,'data_final_yrs',obj.data_final_yrs);
            
                % data run params
            addParameter(p, 'fnames',obj.fnames);
            addParameter(p, 'default_fnames',obj.default_fnames);
            addParameter(p, 'datadir',obj.datadir);
            addParameter(p, 'outname',obj.outname);
            addParameter(p, 'outdir',obj.outdir);

            addParameter(p,'interp_method', obj.interp_method);
            
            addParameter(p,'yrlen', obj.yrlen);
            addParameter(p,'file_lats', obj.file_lats);
            addParameter(p,'file_lons', obj.file_lons);
            addParameter(p,'lats', obj.lats);
            addParameter(p,'lons', obj.lons);
            addParameter(p,'stninfo',obj.stninfo);
            addParameter(p,'model', obj.model);
            addParameter(p,'varname', obj.varname);
            addParameter(p,'obsvname', obj.obsvname);
            addParameter(p,'varlongname', obj.varlongname);
            addParameter(p,'ensemble', obj.ensemble);
            addParameter(p,'scenario', obj.scenario);
            addParameter(p,'scaling', obj.scaling);
            addParameter(p,'cdf_append_pts', obj.cdf_append_pts);
            addParameter(p,'cdf_append_type',obj.cdf_append_type, @(s) strlength(s)==0 || any(strcmpi(s,["sk_normal","empirical","normal"])));
            
            
            addParameter(p,'runType',obj.runType);
%           addParameter(p,'Unmatched',[]);
            addParameter(p,'monthly_valid_count',obj.monthly_valid_count);

            
            addParameter(p,'prcp_min',obj.prcp_min);
                        
            
            addParameter(p,'runID',obj.runID);
            addParameter(p,'runLbl',obj.runLbl);
            addParameter(p,'stnID',obj.stnID);
            addParameter(p,'stnName',obj.stnName);
            addParameter(p,'region',obj.region);
            addParameter(p,'llgrid_size',obj.llgrid_size);
            addParameter(p,'llgrid_lbl',obj.llgrid_lbl);
            
            addParameter(p,'ARRM_V2_version',obj.ARRM_V2_version);
            addParameter(p,'ARRM_V2_dir',obj.ARRM_V2_dir);
            addParameter(p,'FillValue',obj.FillValue);
            addParameter(p,'to_struct',obj.to_struct);          % set to true to get output as structs rather than classes.  Used for disaggregation code.         
            addParameter(p,'keep_all',obj.keep_all);            % set to true to retain trimmed data when calling DA.trim_data_yrs(...)  Used for disaggregation code.         
            addParameter(p,'extended',obj.extended);            % set to true to save probabilities for each downscaled output point.
            addParameter(p,'zval_offset',obj.zval_offset);      
            addParameter(p,'zval_scaling',obj.zval_scaling);    

            addParameter(p,'do_log',obj.do_log);
            addParameter(p,'logname', obj.logname);
            addParameter(p,'do_parfor', obj.do_parfor);
            addParameter(p,'fidlog', obj.fidlog);

            addParameter(p,'progname', obj.progname);
            addParameter(p,'exit_on_error', obj.exit_on_error);
            
            addParameter(p,'displayResults', obj.displayResults);
            addParameter(p,'figbase', obj.figbase);
            addParameter(p,'figflags', obj.figflags);
            addParameter(p,'figname', obj.figname);
            addParameter(p,'figext', obj.figext);
            addParameter(p,'figpos', obj.figpos, @(x) isempty(x) || (isnumeric(x) && size(x,2)==4));  
            
            addParameter(p, 'debug_flag', obj.debug_flag);
            
            addParameter(p,'internal_plotflag',obj.internal_plotflag);      % sets internal_plotflag directly with bit pattern.   make flag odd to turn on global!
            addParameter(p,'do_internal_plots',[]);                         % pass in list of flags to set.
                                                                            % bit 0 is global on/off
                                                                            %
                                                                            % bit 1-n turns on plots of  creating PDF
                                                                            % surfaces for basic years or rolling year
                                                                            % sets. 1 for basic years, 1,2,3... for
                                                                            % rolling years set 1,2, 3...
                                                                            %
                                                                            % bits [0,32] turns on plots for appending
                                                                            % extended pdf shapes beyond 1e-4.
            parse(p, varargin{:});
            
            
                % copy parser output into data fields.
            f=fieldnames(p.Results);
            for i=1:length(f)
                if (isprop(obj, f{i}))
                    res=p.Results.(f{i});
                        % convert anything passed in as cell array of chars to array of strings                    
                    if (iscell(res) && ischar(res{1}))
                        res=string(res);
                    end
%                     if (strcmp(f{i},'outdir'))
%                         fprintf("here!\n");
%                     end
                    obj.(f{i}) = res;
                end
            end
            obj.Unmatched = p.Unmatched;    
            
            if (isempty(obj.ARRM_V2_dir))
                [mydir,~] = fileparts(mfilename);
                 obj.ARRM_V2_dir = mydir;
            end
                        
             if (~isempty(p.Results.do_internal_plots))
                obj = obj.set_internal_plotflag(p.Results.do_internal_plots); % see notes above in internal_plotflag
             end
            
            if (~isempty(obj.figext))
                obj.figext = obj.check_fig_extensions(obj.figext);
            end
            
                % make sure we keep trend_yrs empty for precip data.
            if (obj.isPrecipRun)
                obj.trend_yrs = [];
            end
            
                % and make sure lats & lons are of type double
            obj.lats = double(obj.lats);
            obj.lons = double(obj.lons);
            obj.file_lats = double(obj.file_lats);
            obj.file_lons = double(obj.file_lons);
            
            if (istable(obj.stninfo))
                obj.stninfo.lat = double(obj.stninfo.lat);
                obj.stninfo.lon = double(obj.stninfo.lon);
            end

        end        
    end
    
    methods(Access=private,Static)
         function exts = check_fig_extensions(exts)
            if (isempty(exts)), return; end
            exts = string(exts);
            for i=1:length(exts)
                if (extractBefore(exts(i),2) == '.')
                    exts(i) = extractAfter(exts(i),1);
                end
            end
        end
        
    end       
        
end

