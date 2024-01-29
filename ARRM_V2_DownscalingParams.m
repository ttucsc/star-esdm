classdef ARRM_V2_DownscalingParams < ARRM_V2_DataParams
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
                % for downscaling obs, historical and model data
        obsnames        = [];
        default_obsnames= [];
        mdlnames        = [];
        default_mdlnames= [];
        histnames       = [];   % if empty, uses mdlnames
        default_histnames=[];
        obsdir          = [];
        mdldir          = [];
        histdir         = [];   % if empty, uses mdldir
        hostname        = [];
                
        globatts_name   = [];
        NC_atts         = [];
        
        recalc_pdfs     = false;    % if true, recalculate PDFs & CDFs of anomalies (& poss. problines) after downscaling.
                                    % results will be slightly different than if you re-disaggregate the mapped output
                                    % because this uses the adjusted climatology and moving climatology, rather than
                                    % generating new climatologies based on the mapped output.
%       obsvname        = [];   % variable name for observations, if different from varname.  Located in DataParams.

        obs_interp_method = 'closest';  % main DSP's interp_method is used for Model & Hist.
                                        % this defines how to interpolate the obs data.
        
        hist_yrs        = [];   % years of data to use for hist base period.       
%                                 % Specify if different from base_yrs, which is used for obs base period.
%                                 % NOTE: This parameter should be removed. elsewhere I'm forcing use of base_yrs for the historical period.
        
        dspType         = '';   % 'obs','hist','model', 'downscaling','global' or ''.  Flags what DSP is being used for.
        saveMetaData    = true;
        plotFlags       = false(1,7);
        fignum          = [];
        gridbox         = [];       % gridbox to be processed. [lat1, lon1, lat2, lon2] This should probably be in DataParams, along with llgrid_str.
        lastgrid        = false;
        gridpt_loc      = [];       % location of location being downscaled

                                % list of fields to exclude when saving Params to netCDF file
        ncExcludes      = {'fnames', 'datadir', 'default_obsnames', 'default_mdlnames', 'default_histnames', 'globatts_name', 'NC_atts', ...
                             'saveMetaData', 'do_log', 'exit_on_error', 'plotFlags', 'fignum', ...
                             'ncExcludes','UsingDefaults', 'lastgrid'};
    end
    
    properties(Dependent)
                
    end
        
    properties(Access = private)
        
        UsingDefaults = [];
        
    end
        
    methods
        function obj = ARRM_V2_DownscalingParams(varargin)
            
                % Note for changes:
                %   I should parse varargin for DSP-specific params to strip them out...then apply them afterwards.
                % This will reduce the confusion over having to reapply changes.  
                %   Order should be:
                %       1.  peel off the DSP-specific varargins w/ parser
                %               per MATLAB documentation, this should be OK as long as we don't modify the obj yet.
                %       2.  call DP constructor with remaining params
                %       3.  call localized/TTUCSC function and apply
                %       4.  apply DSP-specific params from parsing
                %       5.  replace keywords.
            ix=[];    
            for i=1:2:length(varargin); if (strcmpi('runType',varargin{i})), ix=i; end; end
            if (~isempty(ix))
                runType=string(varargin{ix+1});
%                 varargin(ix:ix+1) = [];         % don't want to remove this, because it is used in both RunParams and DataParams
            else
                runType="unknown"; 
            end

            obj@ARRM_V2_DataParams(varargin{:});    % call superclass constructor to initialize all the DP stuff
            
            args = obj.Unmatched;                   % and remember any unmatched parameters for later.
            obj.Unmatched = [];
            
                                                    % get any localized parameters
            if (exist('ARRM_V2_DownscalingParams_TTUCSC.m','file'))
                local_parms = ARRM_V2_DownscalingParams_TTUCSC(runType);
            elseif (exist('ARRM_V2_DownscalingParams_localized.m','file'))
                local_parms = ARRM_V2_DownscalingParams_localized(runType);
            else
                local_parms = [];
            end
            
                                                    % apply any DP-specific local params to the DP data members.
                                                    % We do this so the DSP's local file can override DP params if needed.
                                                    % Only DSP-specific local params will be in obj.Unmatched now.
            if (~isempty(local_parms))
                obj = obj.update(local_parms);
            end
            
            obj = obj.DSP_update(args);             % and finally, apply any DSP-specific command-line params.
            
%             
%                 % apply any remaining unmatched params, which should all now be DSP-specific.
%             obj = obj.DSP_update(obj.Unmatched);
%             if (~isempty(fieldnames(obj.Unmatched))), error('error:  unrecognized parameters from ARRM_V2_DownscalingParams_TTUCSC'); end
%             
%             if (exist('ARRM_V2_DownscalingParams_TTUCSC.m','file'))
%                 local_parms = ARRM_V2_DownscalingParams_TTUCSC(obj.runType);
%                 obj = obj.Update(local_parms);     %Use Update here, so local_params can override 
% %               obj = obj.DSP_update(obj.Unmatched);                
%                 obj = obj.DSP_update(args);                
%             elseif (exist('ARRM_V2_DownscalingParams_localized.m','file'))
%                 local_parms = ARRM_V2_DownscalingParams_localized(obj.runType);
%                 obj = obj.DSP_update(local_parms);
% %               obj = obj.DSP_update(obj.Unmatched);                
%                 obj = obj.DSP_update(args);                
%                 if (~isempty(fieldnames(obj.Unmatched))), error('error:  unrecognized parameters from ARRM_V2_DownscalingParams_localized'); end
%             end
                % now apply any command-line params (which are in Unmatched);
                
                % first, check for a DSP in the varagins
%           obj = obj.parse_for_DSP(args); 
                        
%             obj = obj.update(args, obj.Unmatched,varargin{:});    % also reapply any fields that might have been overridden by parameters in ARRM_V2_localized(...);
%             
%             obj = obj.replace_keywords();      % and apply any keyword replacements   % now done in DSP_update(...)

            for i=1:length(obj.ncExcludes)      % add ncExcludes to ncExcludeList in parent class.
                prop = obj.ncExcludes{i};
                if (~ismember(prop, obj.ncExcludeList))
                    obj.ncExcludeList{end+1}=prop;
                end
            end
        end
        
       function obj = update(obj, varargin)
            
            
            % Updates DownscalingParams object.  Uses current values as defaults, and parses any key/value pairs
            % Call this when you want to set multiple parameters with a single function call.
            % To update a single parameter, just set the parameter directly.
    
            % I need to add verification to these still, Ian!

            obj = update@ARRM_V2_DataParams(obj, varargin{:});
                        
            obj = obj.DSP_update(obj.Unmatched);
       end
       
       function obj = finalize(obj, varargin)
            % does final check on obj, filling in defaults, etc. as needed.
            % Can be passed additional parameters, & will update first, then finalize.

            
            if (exist('varargin','var') && ~isempty(varargin))
                obj = obj.update(varargin{:});
            end
            
            obj = finalize@ARRM_V2_DataParams(obj);
           
            obj.obsnames=string(obj.obsnames);
            obj.mdlnames=string(obj.mdlnames);
            if (~isempty(obj.histnames)),  obj.histnames=string(obj.histnames); end
            
            obj = obj.make_fnames();
            obj.hostname = get_hostname();
            
            if (~isempty(obj.stninfo) || ~isempty(obj.stnID))   % can't use obj.isStationRun yet, because we might have stnIDs but not stnInfo.
                if (isempty(obj.stninfo)), obj.stninfo = obj.obsnames; end  % get filename from obsnames if stninfo is empty.
                if (ischar_s(obj.stninfo))
                    obj.stninfo = QC_get_site_table(obj.stninfo, "stnID", obj.stnID);
                end
                if (isempty(obj.obsnames) && isprop(obj.stninfo, "UserData"))
                    obj.obsnames = obj.stninfo.Properties.UserData.ncName;
                end
                if (isempty(obj.obsvname) && isprop(obj.stninfo, "UserData") && size(obj.stninfo.Properties.UserData.varName,1)==1)
                    obj.obsvname = obj.stninfo.Properties.UserData.varName;
                end
                if (isempty(obj.obsvname))
                    obj.obsvname = QC_station_varname(obj.varname);
                end
            end
                            % make filenames absolute, rather than relative.
            obj = obj.update_filenames('obsdir','obsnames');
            obj = obj.update_filenames('mdldir','mdlnames');
            obj = obj.update_filenames('histdir','histnames');
            obj = obj.update_filenames('outdir','outname');
            
            if (~obj.isPrecipRun && isempty(obj.trend_yrs))
                trend_yrs = max_yr_range(obj.base_yrs, obj.hist_yrs, obj.rolling_yrs);
%               trend_yrs = max_yr_range(obj.base_yrs,               obj.rolling_yrs);  % hist_yrs now removed.  icsf 2/2021.  
%                                                                                         reinstated icsf 4/21, but ignored because forced to match base_yrs
                obj = obj.update('trend_yrs',trend_yrs);
            end
    
                    % if user specified obsnames, histnames or mdlnames, update the fnames field.

            if     (strcmp(obj.dspType,"obs")),   obj.fnames = obj.obsnames;  
            elseif (strcmp(obj.dspType,"hist")),  obj.fnames = obj.histnames; 
            elseif (strcmp(obj.dspType,"model")), obj.fnames = obj.mdlnames;  
            end
            
            if (isempty(obj.varname) || isempty(obj.model) || isempty(obj.ensemble) || isempty(obj.scenario))
                obj = obj.parse_fnames();
            end
                                    
            if (isempty_s(obj.obsvname)), obj.obsvname=obj.varname; end
            
            obj = obj.replace_keywords();
        end
            
        function [ok, badnames] = check_files_exist(obj, props)
            if (~exist('props','var'))
                props=["obsnames","histnames","mdlnames"];
            end
            badnames=[];
            
            if (nargout == 0)
                check_files_exist@ARRM_V2_DataParams(obj, props);       % will throw an exception if any file doesn't exist
                ok=true;
                return;
            elseif (nargout == 1)
                ok = check_files_exist@ARRM_V2_DataParams(obj, props);  % returns false as soon as it finds a file doesn't exist
                return;
            else
                [ok, badnames] = check_files_exist@ARRM_V2_DataParams(obj, props);  % checks all files for existence and returns list of missing as well
            end
            
        end
        
        function obj = make_fnames(obj, prop)
                % replaces keywords in obsnames, mdlnames and histnames, using the default_names as the template if ?name is empty.
                % This function is called by finalize().  If you don't call finalize() yourself, then
                % call this function later if you want to assemble the filenames on the fly from the
                % keywords (model, varname, ensemble, scenario, etc.).
                % To do this, supply a default filename in the ARRM_V2_DataParams_localized(...) function
                % You'll probably want to do this separately for station runs and gridded runs.
           if (~exist('prop','var') || isempty(prop))
                prop = ["obsnames","mdlnames","histnames"];
            else
                prop = string(prop);
            end
            for i=1:length(prop)
                obj = make_fnames@ARRM_V2_DataParams(obj, prop(i));
            end
        end
        
    end     % end methods, public
    
    methods(Access = private)
        
        function obj = DSP_update(obj, varargin)
            
%                     % first, check for a DSP in the varargins.
%                         % parse input for DSP
%                         
%             obj = obj.parse_for_DSP(varargin{:});
%             args = obj.Unmatched;
%             if (isempty(fieldnames(args))), return; end
%             
                % and parse the remaining arguments.
            p = inputParser;
            p.StructExpand = true;
            p.KeepUnmatched = true;
    
                % data run params
            addParameter(p, 'obsnames',obj.obsnames);
            addParameter(p, 'mdlnames',obj.mdlnames);
            addParameter(p, 'histnames',obj.histnames);
            addParameter(p, 'default_obsnames',obj.default_obsnames);
            addParameter(p, 'default_mdlnames',obj.default_mdlnames);
            addParameter(p, 'default_histnames',obj.default_histnames);
            addParameter(p, 'obsdir',obj.obsdir);
            addParameter(p, 'mdldir',obj.mdldir);
            addParameter(p, 'histdir',obj.histdir);
            
            addParameter(p, 'globatts_name',obj.globatts_name);
            addParameter(p, 'NC_atts',obj.NC_atts);
            addParameter(p, 'recalc_pdfs',obj.recalc_pdfs);

%           addParameter(p,'obsvname', obj.obsvname);   % moved back to DataParams.
            
            addParameter(p, 'hist_yrs',obj.hist_yrs);
            addParameter(p, 'dspType',obj.dspType, @(x) ismember(x,{'obs','hist','model','downscaling','global',''}));
            addParameter(p, 'obs_interp_method',obj.obs_interp_method);
            
            addParameter(p, 'do_log',obj.do_log);
            addParameter(p, 'exit_on_error',obj.exit_on_error);
            addParameter(p, 'plotFlags',obj.plotFlags);
            addParameter(p, 'fignum',obj.fignum);
            addParameter(p, 'saveMetaData',obj.saveMetaData);
            addParameter(p, 'gridbox',obj.gridbox);
            addParameter(p, 'lastgrid',obj.lastgrid);
            addParameter(p, 'gridpt_loc',obj.gridpt_loc);
            
            addParameter(p,'runType',obj.runType);
%           addParameter(p,'Unmatched',[]);
            
%             addParameter(p,'monthly_valid_count',obj.monthly_valid_count);  % moved to DataParams.  icsf.
                        
%             parse(p, args);
            parse(p, varargin{:});
            
            f=fieldnames(p.Results);
            for i=1:length(f)
                if (isprop(obj, f{i}))
                    obj.(f{i}) = p.Results.(f{i});
                end
            end
            obj.Unmatched = p.Unmatched;
            
            obj.UsingDefaults = p.UsingDefaults;
            
            if (any(strcmp(obj.dspType,{'global','obs'})) && ((isempty(obj.stninfo) || ~isempty(obj.stnID))))   % can't use obj.isStationRun yet, because we might have stnIDs but not stnInfo.
                if (isempty(obj.obsnames))
                    obj = obj.make_fnames('obsnames');
                    obj = obj.update_filenames('obsdir','obsnames');
                end
                try
                    if (~isempty(obj.stninfo) && (~istable(obj.stninfo) || size(obj.stninfo,1) ~= length(string(obj.stnID)) || any(~strcmp(obj.stnID, obj.stninfo.stnID))))
                        obj.stninfo = QC_get_site_table(obj.obsnames,"stnID", obj.stnID);
                    end
                catch
                    obj.DP.warn_log(oops());
                end
            end            
            
        end
      
                    
%        function obj = parse_for_DSP(obj, varargin)
%             p = inputParser;
%             p.StructExpand = true;
%             p.KeepUnmatched = true;    
%             addParameter(p,'DSP',[]);
%             parse(p, varargin{:});
%             DSP = p.Results.DSP;
%             if (~isempty(DSP) && ~isa(DSP,'ARRM_V2_DownscalingParams'))
%                 error("error:  DSP parameter is not an ARRM_V2_DownscalingParams object");
%             end           
%                     % copy fields into self.
%             props = properties(DSP);
%             for i=1:length(props)
%                 try
%                     obj.(props{i}) = DSP.(props{i});     % in a try block so loop doesn't abort on dependent properties.
%                 catch
%                     fprintf('oops\n');
%                 end
%             end
%             obj.Unmatched = p.Unmatched;        % save rest of input params for later
%         end
        
    end         % end methods, private
   
    
end

