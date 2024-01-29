function hist_att = ncdf_update_history(ncobj, histinfo, tstamp)
%
%   updates history attribute (or creates it) for ncdf object ncobj
%   tstamp is optional timestamp to add; should be datenum object, such as now();
%

    if (~isa(ncobj, "ncObj"))
        error("error: ncdf_update_history: first param is not an ncObj");
    end
    
    if (exist('tstamp','var') && ~isempty(tstamp))        
        timestr = datestr(tstamp, "yyyy-mm-dd HH:MM") + ": ";
    else
        timestr = "";
    end
        
    try
        hist_att = ncobj.get("/Attributes/history");
        hist_att.Value = sprintf("%s%s; \n%s",timestr, histinfo, hist_att.Value); 
    catch
        hist_att = Attribute("history",sprintf("%s: %s; \n%s",datestr(tstamp,'yyyy:mm:dd HH:MM'), histinfo));
    end
    ncobj.put("/Attributes", hist_att);

end

