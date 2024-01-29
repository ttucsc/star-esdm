function elev = gdem2_elevation(lat,lon, verbose, in_feet, gdem_dir)
% elev = gdem2_elevation(lat,lon, verbose, in_feet, gdem_dir)
%
%   Returns elevation(s) for lat,lon, and flag(s) for locations that are not on a continent or an island.
%   Elevations are taken from ASTER GDEM2 database, which has 1-degree gridcells with 1-arcsecond elevations for all 
%   land locations.  Lookup is into 1-degree gridcell file, and returns the elevation of closest pixel for each lat/lon.
%   This should be within approx. 30 m (100 ft) of any location on the earth.  Comparing with Google Maps lat/lons,
%   the locations agree within +/-1 pixel/30m/100ft.  GDEM2 elevations over ocean is 0.
%
%   GDEM2 accuracy is reported as 14.14 m. (~50 ft)
%
%   Not all gridcells are present in the GDEM database.  It covers latitudes -83 to +83, and does not have a gridcell
%   file for any gridcell with no islands or continent area.  
%   If no gridcell is available, for the lat/lon, an elevation of 0  is returned.

%   GDEM values over large bodies of water, such as the great lakes, are not reliable.
%
%   NOTE:  elevations for the northern tip of Greenland, and the center of Antarctica, will return elevations of 0.
%
%   This program expects the GDEM2 database of geotiff files to be in the folder given by variable gdem_dir.
%
%   Inputs:  
%       lat, lon    same-length arrays giving lat & lon of location(s) to find elevations for.
%       in_feet     boolean.  If true, results returned in feet.  else in meters.  default:  false
%       gdem_dir    directory of where to find GDEM2 geotiff files.  default:  /Volumes/jcsf_data/GDEM2/dems
%       verbose     false/true/0/1/2    falseor 0: no warnings, no progress bar
%                                       true or 1: progress bar
%                                       2:         warnings/error messages and progress bar 
%
%   Outputs
%       elev        single value or column matrix of elevations, one for each lat/lon input pair
%                       any location without a valid elevation is set to 0.  
%                       This is correct for oceans, but not correct for places like the great lakes.
%                   
%
%   Author:  Ian Scott-Fleming, TTU Climate Center, Texas Tech University, Lubbock, TX
%

    if (~exist('verbose', 'var') || isempty(verbose)),  verbose = false; end
    if (~exist('in_feet', 'var') || isempty(in_feet)),  in_feet = false; end
    if (~exist('gdem-dir','var') || isempty(gdem_dir)), gdem_dir = "/Volumes/lacie_1/jcsf_data/GDEM2/dems"; end
    
    npts    = length(lat);
    fnames  = strings(npts,1);
    latchar = repmat('S',npts,1);
    lonchar = repmat('W',npts,1);
    elev    = zeros(npts,1); 
    
    lon = mod(lon+180,360)-180;     % make sure lons are in range [-180.180)
    latdeg = floor(lat);            % latitude degree, used to select filename
    latdd  = lat-latdeg;            % fraction of degree, used to calculate pixel
    londeg = floor(lon);
    londd  = lon-londeg;
    
    latrow = round((1-latdd)*3600)+1;   % pixel in file
    latrow = max(1,min(3601,latrow));   % make sure it's in range...
    loncol = round(londd*3600)+1;        
    loncol = max(1,min(3601,loncol));
    
    latchar(lat>=0) = 'N';
    lonchar(lon>=0) = 'E';
    
    latdeg = abs(latdeg);
    londeg = abs(londeg);

    for i=1:npts
        fnames(i) = fullfile(gdem_dir,sprintf("ASTGTM2_%c%02d%c%03d_dem.tif", latchar(i),latdeg(i),lonchar(i),londeg(i)));
    end
    
     
    if (npts==1)
        try
            tif=geotiffread(fnames(1));
            elev = tif(latrow,loncol);
        catch
            elev = nan(3601,3601);
        end
    else                
        ffnames = unique(sort(fnames));
        nffn = length(ffnames);
        
        if (verbose)
            fprintf("reading %d elevation files to find %d elevations\n", nffn, npts);
        end
        
        for j=1:nffn
            ix = find(strcmp(fnames,ffnames(j)));    % find all locations in this gridbox
            try
%               fprintf("processing sector %4d of %4d:  %s\n", j, nffn, ffnames(j));
                tif = geotiffread(ffnames(j));
            catch
                if (verbose==2)
                    fprintf("\n")
                    fprintf("error reading geotiff file %s: ", ffnames(j));
                    for i=1:min(20,length(ix))
                        iix = ix(i);
                        fprintf('%9.4f %9.4f ; ',  lat(iix),lon(iix));                    
                    end
                    if (length(ix) > 20)
                        fprintf('and %d more\n', length(ix)-20);
                    else
                        fprintf('\n');
                    end                
                end
                tif = nan(3601,3601);
            end
%             elev(ix) = tif(latrow(ix),loncol(ix));
            for i=1:length(ix)
                iix = ix(i);
                elev(iix) = tif(latrow(iix),loncol(iix));
            end
            if (verbose), show_progress(j, nffn); end
        end

        
    end   
    
    if (in_feet)
        meters2feet = 100/2.54/12;     % 1m / (2.54 cm/inch)
        elev = double(elev) * meters2feet;
    end
    elev = round(elev);
end

