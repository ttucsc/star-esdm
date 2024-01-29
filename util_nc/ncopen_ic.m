function [ ncid ] = ncopen_ic( fname, mode )
    
    global ncids;
    
    if (exist('mode','var'))
        ncid = netcdf.open(fname, mode);
    else
        ncid = netcdf.open(fname);
    end
    
    ncids(end+1)=ncid;
    
end

