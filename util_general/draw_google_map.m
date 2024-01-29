function h = draw_google_map(lats, lons, connected, fignum, apikey)
% function h = draw_google_map(lats, lons, fignum, connected)
%
%   Ian Scott-Fleming
%   Texas Tech University
%
%   function to plot one or more sets of lat/lon coordinates on top of a google map of the area.
%   Google map is scrollable and (to an extent) zoomable
%       Code uses plot_google_map(...) from Matlab File Exchange (see notes below)
%   
%   Inputs:
%       lats        vector of latitudes  (for points or single line) or 
%                       matrix or cell array of vectors (for multiple lines)
%       lons        vector of longitudes (for points or single line) or 
%                       matrix or cell array of vectors (for multiple lines)
%       connected   boolean (single or vector, 1 per line, controlling whether to join points or just mark them
%   
%       fignum      figure number.  If empty, opens new figure.  If absent, draws in last figure
%
%   Code uses plot_google_map(...) from Matlab File Exchange, by Zohar Bar-Yehuda
%   https://www.mathworks.com/matlabcentral/fileexchange/27627-zoharby-plot-google-map
%
%       Note:   I added the following code to plot_google_map's argument parsing case statement,
%               along with an if(do_back)  to conditionally move the map to the background.
%               to allow moving the map to the foreground 
%               with transparency, because at present Matlab's line plots do not have a transparency option.
%
%              case 'back'      % added icsf
%                 do_back = varargin{idx+1};        
%
%              if (do_back)    % if and end added, icsf.  (unconditional) uistack(...) was already there
%                  uistack(h,'bottom') % move map to bottom (so it doesn't hide previously drawn annotations)
%              end
%
%-------------------------------------------------------------------------------------

        % parameter defaults
    if (exist('fignum','var'))
        if (isempty(fignum))
            h = figure;
        else
            h = figure(fignum);
        end
    end
    if (~exist('connected','var') || isempty(connected))
        if (~iscell(lats) && length(lats)==2)
            connected=false;
        else
            connected = true;
        end
    end
    if (~exist('apikey','var') || isempty(apikey))
        apikey = get_api_key('maps');
    end
    
    if (iscell(lats))
        nsets = length(lats);
        mylats = lats;
        mylons = lons;
    else            % get lats & lons into cell arrays if they are in matrices
        if (~iscolumn(lats)),  lats = lats'; end
        if (~iscolumn(lons)),  lons = lons'; end
        nsets = size(lats,2);
        mylats = cell(nsets,1);
        mylons = cell(nsets,1);
        for i=1:nsets
            mylats{i} = lats(:,i);
            mylons{i} = lons(:,i);
        end
    end
        
    if (length(connected) == 1 && nsets > 1)
        connected = repmat(connected, nsets,1);    % duplicate so one for each set of points
    end
    
    if (sum(connected)==0)  % if points are just disconnected, scatter the points and put map in background 
        for i=1:nsets
            mylons{i}=mod(mylons{i}+540,360)-180;       % google maps wants lons between -180 and + 180...
            scatter(mylons{i}, mylats{i},36,'filled');
            hold on;
        end
        hold off;
        plot_google_map_ic('MapType','map', 'APIkey',apikey,'verbose',true);
        
    else                     % if any points are connected, draw lines, then put map on top w/ transparency
        for i=1:nsets
            mylons{i}=mod(mylons{i}+540,360)-180;       % google maps wants lons between -180 and + 180...
            if (nsets == 1)
                plot(mylons{i}, mylats{i}, 'y-', 'linewidth',6);
            else
                plot(mylons{i}, mylats{i}, 'linewidth',6);
            end
            hold on;
            plot(mylons{i}, mylats{i}, 'k-');
            scatter(mylons{i}, mylats{i},36,'r','filled');
       end
                % note:  this needs the change to plot_google_map(...) as described above.
        plot_google_map_ic('MapType','map','alpha',.5,'back',false, 'ShowLabels', 1, 'APIkey',apikey,'verbose',true);
        hold off;
    end
    
    xlabel('lon');
    ylabel('lat');
end
