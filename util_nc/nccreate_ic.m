function [ ncid ] = nccreate_ic( fname, cmode )
    
    global ncids;
    
    ncid = netcdf.create(fname, cmode);
    
    ncids(end+1)=ncid;
    
end

