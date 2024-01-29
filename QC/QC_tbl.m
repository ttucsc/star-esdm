function tbl = QC_tbl(tbl, stns, search_type)
%   Displays or returns a QC station table, with dates as strings instead of datenums.
    %
    %   tbl             either a stnTbl from QC_get_site_table(...) or the name of a QC netcdf file.
    %   stns            indexes or list of stations to select (optional).
    %                       stns can be numeric row indexes, vector of stnIDs, or site names.  
    %   search_type     can be 'stnName', or 'any', for search_type if stns is not indexes or stnIDs
    %       can preface search type with 'like_' for partial match on name.  See QC_get_site_table

    if (~exist('stns', 'var')), stns = []; end
    if (~exist('search_type','var')), search_type = []; end
    
    if(isempty(stns))
        tbl =  QC_get_site_table(tbl, "showDates",true);
    else
        tbl = QC_get_site_table(tbl, "stnID", stns, "searchType", search_type,"showDates",true);
    end
    
    if (nargout == 0)
        disp(tbl);
        tbl=[];
    end
end

