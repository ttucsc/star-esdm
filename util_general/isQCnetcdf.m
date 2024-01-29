function [yn,nc] = isQCnetcdf(fname)
% returns true if fname is a QC netcdf file.  

    [yn,nc] = isnetcdf(fname, true);
    if (~yn), return; end
    if (~ismember("stnID",{nc.Variables.Name}) || ~ismember("lat",{nc.Variables.Name}) || ~ismember("lon",{nc.Variables.Name})), yn=false; return; end
    yn=true;
end