function ncclose_ic( ncid )
 
    global ncids;   % keeps track of all ncids I've opened using ncopen_ic(...)
    
        % close all known open ncid's, if called with 'all')
    if (~exist('ncid','var') || isempty_s(ncid) || (ischar_s(ncid) && strcmp(ncid,'all')))
        for i=1:length(ncids)
            try
                netcdf.close(ncids(i));
            catch
            end
        end
        ncids=[];
    else
            % close a specific ncid:
        netcdf.close(ncid);
            % see if ncid is in global list of open ncids
        ix=find(ncids==ncid);
        if (~isempty(ix))
            ncids(ix)=[];
        end
    end
end

