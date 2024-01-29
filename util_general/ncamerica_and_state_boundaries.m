function [nalat, nalon, statelat,statelon] = ncamerica_and_state_boundaries(fignum, do_holdon, linespecs, latbnds, lonbnds) 
% Returns Coastline lats & lons for North America, and US state lats & lons
%   coastlines are from Matlab's coastlines.mat file
%   states are from Matlab's usastatelo.shp shape file
%
%   Unfortunately, the coastlines are lower resolution than the states, so the lines are not identical along the US
%   coasts.  
%
%   Optional inputs:
%       fignum                  if present and not empty, draw the boundaries in the figure.
%       do_holdon               logical.  If true to add to existing map.  
%                                         If false, create worldmap for region identified by latbnds, lonbnds and draw
%                                         boundaries into map
%       linespecs               linespecs for NA boundaries and state boundaries.  ["k-"]
%                                         if only 1 linespec, use for both sets of boundaries.
%       latbnds, lonbnds        row vectors or matrices of lat/lon boundaries to include
%
% ---
%   To create your own map & projection, create the map with one of Matlab's mapping commands.  A reasonable approach is
%   to use 
%       worldmap([lat1,lat2],[lon1,lon2]) ;     % will select a reasonable projection for the lat/lon region specified.
%   and then draw the lines with 
%       [nalat,nalon,statelat,statelon] = nc_america_and_state_boundaries();
%       plotm(nalat,nalon,'k-');
%       plotm(statelat,statelon,'k-');
% ---


    load("coastlines","coastlat","coastlon");
    
    if (~exist('fignum','var') || isempty(fignum))
        fignum = false;
    end
    if (~exist('do_holdon','var') || isempty(do_holdon))
        do_holdon = false;
    end    
    if (~exist('linespecs','var') || isempty(linespecs))
        linespecs = "k-";
    else
        linespecs = string(linespecs);
    end
    us = min(2, length(linespecs));     % indexes for us & canada linespecs.
    can = min(3, length(linespecs));
    
    latbnds = [0, 85];
    lonbnds = [-180, -10];
    
    if (~exist('latbnds','var') || isempty(latbnds))
        latbnds = [7,15; 15,45; 45, 85];
    elseif (~exist("lonbnds","var"))
        error("missing input:  lonbnds");
    end
    if (~exist('lonbnds','var') || isempty(lonbnds))
        lonbnds = [-180,-76.45; -180,-50; -180, -10];
    else
        lonbnds = mod(lonbnds+180,360)-180;
    end
    
    keepers = false(size(coastlat)); %#ok<NODEF>

    for i=1:size(latbnds,1)
        keepers(coastlat >= latbnds(i,1) & coastlat <= latbnds(i,2) & coastlon >= lonbnds(i,1) & coastlon <= lonbnds(i,2)) = true; %#ok<NODEF>
    end
    coastlat(~keepers) = nan;
    coastlon(~keepers) = nan;
    edges = false(size(keepers));
    for i=1:(length(keepers)-1)
        if ((keepers(i) && ~keepers(i+1)))
            edges(i+1) = true;
        end
    end
    keepers = keepers | edges;
    
%     dupnas = false(size(keepers));
%     for i=2:length(keepers)-1
%         if (~keepers(i-1) && ~keepers(i) && ~keepers(i+1))
%            dupnas(i) = true;
%         end
%     end
% 
%     difk = ~logical(diff(keepers));       % will be true when bother keepers(i) && keepers(i+1) are false.  (i.e., when coastlat(i) && coastlat(i+1) are both NAs)
%     dupnas = [difk;false] & ~keepers;     % we don't want to keep repeating NAs, so we flag them and remove them.
%   keepers = keepers & ~dupnas;
    nalat = coastlat(keepers);
    nalon = coastlon(keepers);
    
    states=shaperead("usastatelo","UseGeoCoords",true);
    
    statelat=[];
    statelon=[];
    for i=1:length(states)
        statelat=cat(2,statelat,states(i).Lat);
        statelat(end+1)=nan; %#ok<AGROW>
        statelon = cat(2,statelon,states(i).Lon);
        statelon(end+1)=nan; %#ok<AGROW>
    end
    
    load("canada_provinces.mat", "all_provinces");
        
    if (fignum || do_holdon)
        if (fignum), figure(fignum); end
        if (do_holdon)
            hold on;
        else
            ltbnds = [min(latbnds(:,1)),max(latbnds(:,2))];
            lnbnds = [min(lonbnds(:,1)),max(lonbnds(:,2))];
            worldmap(ltbnds,lnbnds);
        end
        plotm(nalat,nalon,linespecs(1));
%       plotm(coastlat,coastlon,linespecs(1));
        plotm(statelat,statelon,linespecs(us));
        plotm(all_provinces.lat, all_provinces.lon, linespecs(can));
    end
    
    if (nargout == 0)
        nalat = [];
        nalon = [];
        statelat = [];
        statelon = [];
    else
            % make these column vectors.  They're row vectors at the moment...
        statelat = to_column(statelat);
        statelon = to_column(statelon);
    end
end

