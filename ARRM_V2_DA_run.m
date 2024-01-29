function DA_results = ARRM_V2_DA_run(varargin)
%   Front-end for ARRM_V2_DisaggregateSignal(...) on one or more stations or lat/lon pairs of model data.
%
%   Program to run all latitudes for a single longitude through ARRM_V2
%   Can be run standalone, or started via ARRM_V2_wrapper(...)
%
%   Inputs (all inputs except varargin are required):     
%                   NOTE:  variable name, ensemble, etc. can be picked up from the netcdf filename or the file for both 
%                           station and gridded data.  
%                           Current list of variables is found in the call to function ncfind_varinfo_ic(...)  
%                           below, which finds the first variable matching any name in the list.
%                   
% 
%       varargin:   optional keyword/value pairs
%
%       varargin must specify either fnames or model, varname, ensemble, & scenario
%
%       this needs updating, Ian!
%                       any of the keyword/value pairs listed in the ARRM_RunParams class.  
%                       Type "help ARRM_RunParams" for more info.
%                   plus
%                       'rolling_yrs', [yr1,yr2]    years over which to calculate moving climatology, pdfs & CDFs
%                                                       set start year to 1 or true to start at beginning of 1st year of valid data for site
%                                                       set end year to 1 or true to end at end of last year of valid  data for site
%                                                       example:  [1950,2100].  default:  [true,true]
%                       'base_yrs', [yr1,yr2]       to specify the years to use for the base years (usually historical period) [1950,2005].  
%                                                       Used to calculate the basic climatology, pdfs and CDFs,
%                                                       and the basic_probability isolines.
%                       'varname', varname          specify variable name in netcdf file to process
%                                                       if not specified, will extract from filename or use 1st
%                                                       non-dimension variable in file.                                                       
%                                                       Current list of variables is found in the call to function ncfind_varinfo_ic(...)  
%                                                       below, which finds the first variable matching any name in the list.
%                       'interp_method', 'method'   interpolation method for gridded data.  'method' can be any of:
%                                                       'closest'           use closest gridcell only
%                                                       'inverse_distance'  weight gridcells by inverse distance-squared
%                                                       'bilinear_weights'  use bilinear-interpolation weights
%                                                       'interpolate'       interpolate to single data series using bilinear interpolation
%                                                                           default:  'bilinear_weights'
%                       'scaling', scaling          scale variable data with one of the following
%                                                       'linear' (default), 'log','log10' or numeric
%                                                       if numeric, data is scaled  data^(1/pwr)  i.e. pwr defines the ROOT
%                                                       to which the data is raised.  Examples:  2 for sqrt,  3 for cube
%                                                       root, etc.  must be positive, can be fractional  (e.g., 2.5) .
%                       'figbase', fignum           base figure number, to produce graphics of each result.
%                                                       if figbase is provided, 4 figures are generated for each station or lat/lon
%                                                       default:  none (no figures)
%                       'figflags', figflags        vector of 6 boolean flags selecting which figures to draw.
%                                                       if missing, and figbase given, draws all 5 figures.
%                                                              Figures are:
%                                                                   1.  raw data  && raw data surface
%                                                                   2.  Static (annual) climatology
%                                                                   3.  Moving climatology, moving climatology surface and diff from annual climatology
%                                                                   4.  Anomaly PDF & CDF surfaces
%                                                                           if multiple PDFs & CDFS done, then surfaces are annimated over time
%                                                                   5.  Overall PDF & CDF with probability lines drawn in.
%
%                   If data is not temperature data, then you should specify
%                       'pdfEdgeSpec', [start,step,end]  or vector of edges
%                                   edges (bottom values) specification to use for histogramming to create the pdfs.
%                                   Given as [minval, stepsize, maxval], as in [-50,0.1,50];
%                                   default:  [-50:(calculated):50]
%       
%   Output:
%       results         results struct from  jc_signal_decomposition(...), if only 1 location specified, or
%                                   cell array of results structs, one for each location specified
%
%                       NOTE:  if data is scaled (scaling is anything other than 'linear'), results are in the
%                               scaled data space.  This should be fixed at some point, Ian...
% 
%-----------------------------------------------------------------------------------

    if (isempty(varargin))
        help(mfilename('fullpath'));
        return;
    end
    
       % make sure the path is set to find all the code we need.      
     ARRM_V2_setpath(mfilename('fullpath')); 

        % parse input parameters and return as individual variables ( matlab's parse(...) returns them as a struct.
  
    [RP, DP, DA_title] = initParams(varargin{:}); 
        
    if (DP.isStationRun)       % NOTE:  Station runs Currently only tested for a single station..
%        DP.stnInfo = get_site_info(DP.stations, DP.fnames);
        DP.print_log( 'ARRM_V2_DA_run:  Station run\n');
        DP.print_log( 'input files: \n');
        DP.print_log( '\t%s\n', basename(DP.fnames));
        if (length(DP.lats) <= 3)
            DP.print_log( 'lats:  ');
            DP.print_log( '%9.4f ', DP.lats);
            DP.print_log( '\nlons;  ');
            DP.print_log( '%9.4f ', DP.lons);
        else
            DP.print_log( 'lats range:  ');
            DP.print_log( '%9.4f - %9.4f', min(DP.lats), max(DP.lats));
            DP.print_log( '\nlons range;  ');
            DP.print_log( '%9.4f - %9.4f', min(DP.lons), max(DP.lons));
        end
        DP.print_log( '\n');

        DP.units = DP.stninfo.Properties.UserData.units;
        
        nstns = size(DP.stninfo,1);
        DA_results = cell(nstns,1);
        figbase=DP.figbase;
        for i=1:nstns
            DPr = DP.update('stnID',DP.stninfo.stnID(i),'stnName',DP.stninfo.stnName(i), 'lats',DP.stninfo.lat(i),'lons',DP.stninfo.lon(i),'figbase',figbase);
            lbl = sprintf('%s %s', DPr.stnID, DPr.stnName);
            DP.print_log( 'Station %d: %s\n', i, lbl);
            svec = datevec(DPr.stninfo.startDate(1));        % note:  DPr.stninfo's start and end dates are really index from beginning of file's date range, not the datenum of the actual start/end date
            evec = datevec(DPr.stninfo.endDate(1));
            DPr.data_yrs = min_yr_range(DPr.data_yrs,[svec; evec]);      % will update base_yrs, etc.

            sInfo = QC_get_site_table(DPr.stninfo, DPr.data_yrs(1), DPr.data_yrs(2), "stnID",DPr.stnID, "removeLeaps",true);
            raw_data = scale_data(sInfo.data(1,:), DPr.scaling, 'forward');
            RP.weights = 1;
            if (isempty(DA_title)), DA_title = sprintf('DA station: %s %s', DPr.stnID, DPr.stnName); end
            results = ARRM_V2_DisaggregateSignal(raw_data,RP, DPr,DA_title);  
            if (~isempty(figbase))
                nfigs = jc_plot_data_decomposition(results, lbl);
                figbase = figbase + nfigs;                
                % I should scale_data(...'reverse') here on the anomalies, Ian!
            end
            if (DP.to_struct)
                DA_results{i} = results.toStruct();
            else
                DA_results{i} = results.clone();
            end
            results.clean();

        end
        
    else
                % get start, end year from 1st data file if startYr or endYr not specified. 
        DP.print_log( 'ARRM_V2_DA_run:  model %s, variable %s, ensemble %s, scenario %s\n', DP.model, DP.varname, DP.ensemble, DP.scenario);
        DP.print_log( 'input files: \n');
        DP.print_log( '\t%s\n', basename(DP.fnames));
        if (length(DP.lats) <= 3)
            DP.print_log( 'lats:  ');
            DP.print_log( '%9.4f ', DP.lats);
            DP.print_log( '\nlons;  ');
            DP.print_log( '%9.4f ', DP.lons);
        else
            DP.print_log( 'lats range:  ');
            DP.print_log( '%9.4f - %9.4f', min(DP.lats), max(DP_lats));
            DP.print_log( '\nlons range;  ');
            DP.print_log( '%9.4f - %9.4f', min(DP.lons), max(DP.lons));
        end
        DP.print_log( '\n');

        if (DP.base_yrs(1)==1 || DP.base_yrs(2)==1)
            [timevals, calendar, units] = ncget_timeinfo_ic(fnames(1));
%           [from_vec, timescale, isUTC] = nc_parse_date_str(units);
            from_vec                     = nc_parse_date_str(units);
            dnums = datenum_cal(from_vec, calendar) + timevals;
            svec = datevec_cal(dnums(1),calendar);
            evec = datevec_cal(dnums(2),calendar);
            if (DP.base_yrs(1) == 1), DP.base_yrs(1) = svec(1); end
            if (DP.base_yrs(2) == 1), DP.base_yrs(2) = evec(1); end
        end
        
        % read netcdf files.  extracts lat/lon region & joins data from multiple netcdf files into a single ncdf object.
                    % historical and future data is put together into a single data stream for each lat/lon point, and  calendar to 365-day.
        [nc_in, DP] = read_all_nc_data(DP);  

        nlons = length(DP.lons);
        nlats = length(DP.lats);
        
        DA_results = cell(nlons,nlats);
        figbase = DP.figbase;
        for ilon = 1:nlons
            for ilat = 1:nlats
                lat = DP.lats(ilat);
                lon = DP.lons(ilon);
%                lbl = sprintf('(%8.4f,%8.4f)',lat,lon);
                lbl = DP.runLbl;
                DP.print_log('Location:  %s  (%8.4f,%8.4f)\n', lbl, lat, lon);
                [raw_data, weights] = extract_data(nc_in, DP.varname, lat,lon, DP.interp_method);  % extracts 4 closest gridcells
                if (~strcmpi(DP.scaling,'linear'))
                    raw_data = scale_data(raw_data, DP.scaling, 'forward');
                end
                [npts,nr,nc] = size(raw_data);
                nsets = nr*nc;
                raw_data = reshape(raw_data, npts,nsets);
%                 DP.print_log('model %s: nsets %d   weights:',DP.model, nsets);
%                 DP.print_log('%.10f ', weights);
%                 DP.print_log('\n');
                if (nsets == 1)
                    results = ARRM_V2_DisaggregateSignal(raw_data, RP, DP.update('lats',lat,'lons',lon), sprintf('DA gridded: %.4f %.4f',lat,lon));

                else
                    DAs = cell(nsets,1);
                    minedge = zeros(nsets,1);
                    maxedge = zeros(nsets,1);
                    nbins    = zeros(nsets,1);
                    for i=1:nsets
                        DAs{i} = ARRM_V2_DisaggregateSignal(raw_data(:,i), RP, DP, 'DA_run gridcell',true, false, false, false, false, false);  % do anomalies & calc binning;  we'll do the other steps later. 
                        minedge(i) = DAs{i}.RP.edges(1);
                        maxedge(i) = DAs{i}.RP.edges(end);
                        nbins(i) = length(DAs{i}.RP.edges);
                    end
                            % now find max bin range and number of bins, so all sets use the same binning.
                    for i=1:nsets
                        RP.edges = linspace(min(minedge),max(maxedge), max(nbins));  % slight bug here, Ian.  Should recalculate nbins using smallest delta-bin...
                        RP.do_calc_binning = false;     % make sure we don't redo the binning
                        DAs{i}.RP.edges = RP.edges;
                        DAs{i}.calc_hists();        
                    end
                    
                    
                    if (isempty(DA_title)), DA_title = sprintf('DA weighted gridcells: %.4f %.4f', lat,lon); end
                    results = ARRM_V2_merge_DAs(DAs, RP.update('weights',weights), DP.update('lats',lat,'lons',lon,'figbase',figbase), RP.do_pdfs, DA_title);
%                    for i=1:nsets; DAs{i}.clean(); end  % release memory to avoid memory leaks.
                    keep_all = ~isempty(DP.data_final_yrs) && ~isempty(figbase) && figbase > 0; 
                    if (~isempty(DP.data_final_yrs))
                        results.trim_data_yrs(DP.data_final_yrs, keep_all);    % trims data to just DP.data_final_yrs.  note:  DA_results is a handle object, so no need to assign back to DA_results.
                    end
                end
                if (~isempty(figbase) && any(DP.figflags))
                    nfigs = jc_plot_data_decomposition(results, lbl);
                    figbase = figbase + nfigs;
                end
                        % If data was scaled earlier, I should rescale the data if it was scaled earlier. 
                        % For now, I just return all the results in the scaled data space.
                        % unscaling will take some thought if we've removed trend and climatology...
                        
%                 if (~strcmp(scaling,'linear') && isnan(results.trend_order))
%                     results = unscale_results(results,scaling);
%                 end
                if (results.DP.to_struct)
                    DA_results{ilon,ilat} = results.toStruct();
                else
                    DA_results{ilon,ilat} = results.clone();
                end
                results.clean();  % release memory to avoid memory leaks.
            end
        end
    end
    if (numel(DA_results) == 1)
        DA_results = DA_results{1};
    end
end
        
function [nc_in, DP] = read_all_nc_data(DP)  % reads netcdf files and interpolates data to specified lats and lon
        
    if (ischar(DP.fnames)), DP.fnames=string(DP.fnames); end
    if (isempty(DP.varname))
        DP.varname = ncfind_varinfo_ic(DP.fnames{1}, ["Tmax","Tmin","Prec",'tasmax','tasmin','pr', "temp_F","rh_F","rhsmax","rhsmin","relhum"]);
    end
                                                %                                                  V make this numeric, ian!
                                                %   will use 2 as seed for RNG for 360-day to 365-day expansion.
    nc_in =  ncdf_read_files(DP.fnames, DP.varname,   DP.lats, DP.lons,  DP.data_yrs,   '365-day', 2, DP.varname);      
    
    ncvar = nc_in.get(DP.varname);
    try
        DP.units = ncvar.getattvalue('units');
    catch
        DP.units = '';
    end
    
    try
        DP.varlongname = ncvar.getattvalue('long_name');
    catch
        DP.varlongname = '';
    end
 
end

function [raw_data,weights] = extract_data(ncobj, varname, lat_pt, lon_pt, interp_method)

    lats = ncobj.getvardata('lat');
    lons = ncobj.getvardata('lon');
    data = ncobj.getvardata(varname);
    
    if (strncmpi(interp_method,'clo',3))        % use closest gridcell
        [latix, lonix, ~, ~] = closest(lat_pt, lon_pt, lats, lons, 1, 1);
        raw_data = data(:,latix(1),lonix(1));
        weights=1;
    else
        [latix, lonix, latout, lonout] = closest(lat_pt, lon_pt, lats, lons, 2, 2);     % find closest 4 gridcells

        if (strncmpi(interp_method,'bil',3))              % 'bilinear_sampling':  use standard ARRM_V2 sampling with bilinear weights.
            raw_data = data(:,latix,lonix);
            weights = bilinear_weights_ic(lat_pt, lon_pt, latout, lonout);            
            
                %  special case:  latpt,lonpt is exactly on one of the gridcell centers.
                %  Set weight for that gridcell to almost 1, and set tiny weights for the others based on inverse
                %  distance. This allows us to use alternate sampling when main gridpoint is NA.
                %   This is same as zero-distance case in inverse_distance weighting below.
            if (any(weights==1.0))
                lat_pts = repmat(lat_pt,4,1);
                lon_pts = repmat(lon_pt,4,1);
                dists = distance(lats, lons, lat_pts, lon_pts);
                jx=find(dists==0,1);
                dists(jx) = nan;
                dinv = 1./(dists.^2);   % weights based on 1/square(distance)
                weights = .000001*dinv/nansum(dinv);
                weights(jx) = .999999;
            end

        elseif (strncmpi(interp_method,'inv',3))            % 'inverse_distance':  use standard ARRM_V2 sampling with inverse_distance weights.
            raw_data = data(:,latix,lonix);
            lats  = [latout(1); latout(1); latout(2); latout(2)];
            lons  = [lonout(1); lonout(2); lonout(1); lonout(2)];

            lat_pts = repmat(lat_pt,4,1);
            lon_pts = repmat(lon_pt,4,1);
            dists = distance(lats, lons, lat_pts, lon_pts);
                %  special case:  latpt,lonpt is exactly on one of the gridcell centers.
                %  Set weight for that gridcell to almost 1, and set tiny weights for the others
                %   This allows us to use alternate sampling when main gridpoint is NA.
            if (any(dists==0))  
                jx=find(dists==0,1);
                dists(jx) = nan;
                dinv = 1./(dists.^2);   % weights based on 1/square(distance)
                weights = .000001*dinv/nansum(dinv);
                weights(jx) = .999999;
            else
                dinv = 1./(dists.^2);   % weights based on 1/square(distance)
                weights = dinv/sum(dinv);
            end
        elseif (strncmpi(interp_method,'int',3))    % use bilinear interpolation to get single data series.
            Tlat1lon1 = data(:,latix(1),lonix(1));
            Tlat1lon2 = data(:,latix(1),lonix(2));
            Tlat2lon1 = data(:,latix(2),lonix(1));
            Tlat2lon2 = data(:,latix(2),lonix(2));

            raw_data = interp_bilinear_ic(lat_pt, lon_pt, latout(1), lonout(1), latout(2),lonout(2), Tlat1lon1, Tlat1lon2, Tlat2lon1, Tlat2lon2);
            weights=1;
        else
            error('error:  %s:  unknown sampling method: %s', mfilename, interp_method);
        end
            
    end
end

function [scaled, data_nans] = scale_data(raw_data, pwr, dir)
%   if dir is 'forward':  returns log(pr) or pr^(1/pwr)
%       zeros returned as NAs
%       data_zeros is logical array flagging where original data was 0
%       pr_nans  is logical array flagging where original data was NA
%   if dir is 'reverse':  returns exp(pr) or pr^(pwr)
%       NAs in pr are set to 0, except for locations flagged by was_nans
%       was_nans is logical array of booleans flagging which data to reset to NAs,   

    data_nans = isnan(raw_data);

    if (strcmpi(pwr,'linear'))
        scaled = raw_data;
        return;
    end
    if (strcmp(dir,'forward'))
        if (~isempty(pwr) && isnumeric(pwr))
            scaled = (raw_data).^(1/pwr);
        elseif (strcmpi(pwr,'log10'))
            scaled = log10(raw_data);
        elseif (strcmpi(pwr,'log'))
            scaled = log(raw_data);
        else
            error('error:  bad scaling info');
        end                      
    elseif (strcmp(dir,'reverse'))
        if (~isempty(pwr) && isnumeric(pwr))
            scaled = (raw_data.^(pwr));
        elseif (strcmpi(pwr,'log10'))
            scaled = 10.^(raw_data);
        elseif (strcmpi(pwr,'log'))
            scaled = exp(raw_data);
        else
            error('error:  bad scaling info');
        end
    else
        error('error:  bad scaling info');
    end   
end

function [RP,DP, DA_title] = initParams(varargin)

    % returns an ARRM_V2_RunParams and an ARRM_V2_DataParams object with settings from input arguments.
        
                    % parse input for DA_title
    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = true;    
    addParameter(p,'DA_title',[]);
    parse(p, varargin{:});
    DA_title = p.Results.DA_title;
    Unmatched = p.Unmatched;        % save rest of input params for later
    
    if (isfield(Unmatched, 'RP'))
        RP = Unmatched.RP;
        Unmatched = rmfield(Unmatched,'RP');
        RP = RP.update(Unmatched);
    else
        RP = ARRM_V2_RunParams("temperature", Unmatched);   % we're assuming this is running on temp only.  will check this below.
    end
    Unmatched = RP.Unmatched;
    
    if (isfield(Unmatched,'DP'))
        DP = Unmatched.DP;
        Unmatched = rmfield(Unmatched,'DP');
        DP = DP.update(Unmatched);
    elseif (isfield(Unmatched,'DSP'))
        DP = Unmatched.DSP;
        Unmatched = rmfield(Unmatched,'DSP');
        DP = DP.update(Unmatched);
    else
        DP = ARRM_V2_DownscalingParams(Unmatched);
    end
    Unmatched = DP.Unmatched;
    
    if (DP.isPrecipRun), error("Error.  cannot run on precip currently."); end
    
    if (~isempty(fieldnames(Unmatched)))
        error("unexpected input parameters");
    end
    
    if (DP.isStationRun)
%         file_range=nan(1,2);
%         for i=1:length(DP.fnames)
%             [DP.stninfo, DP.lats, DP.lons] = get_site_info(DP.stations, DP.fnames);
            startvecs = datevec(DP.stninfo.startDate);
            endvecs   = datevec(DP.stninfo.endDate);
%           yr_range = max_yr_range(file_range, [startvecs;endvecs]);
%         end
        yr_range = [startvecs(1,1), endvecs(1,1)];
    else
        file_ranges = ncdf_get_date_range(DP.fnames);
        yr_range = [file_ranges(1,1), file_ranges(2,1)];
    end   
            % limit the range of years to pull and process to what's available in the files
            
    DP = DP.set_yr_limits(yr_range);

    DP.rolling_steps = DP.calc_rolling_steps(RP.pdf_yrstep);
    
            % report any unused input parameters.
            
    if (~isempty(DP.Unmatched))
        fields = fieldnames(DP.Unmatched);
        if (~isempty(fields))
            DP.warn_log( "Warning:  unexpected input parameters:  \n");
            disp(DP.Unmatched);
            error("please correct and rerun");
        end
    end
    DP.lons = mod(DP.lons, 360);  % make sure lons are in range 0-360!
    
end

% function results = unscale_results(results,scaling)
% end
% 
