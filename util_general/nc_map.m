function [gavg,v] = nc_map(fname, varname, ttl, fignum, subpos, ming, maxg, cmap)

    if (~exist('fignum','var')), fignum = []; end
    if (~exist('subpos','var')), subpos = []; end
    if (ischar_s(fname))
        [~,~,ext] = fileparts(fname);
        if (~strcmpi(ext,".nc")),  error("%s is not a netcdf file", fname); end
        nc=ncdf(fname,'create',false);
    elseif (isa(fname,'ncdf'))
        nc = fname;
        fname = nc.Filename;
    else
       error("can't figure out what kind of thing fname is");
    end
    
    if (~exist('cmap','var') || isempty(cmap)), cmap = cmap_CSC(255); end
    
    if (~exist('ttl','var') || isempty(ttl))
        [~,fb,fext] = fileparts(fname); 
        ttl = sprintf("%s%s",fb,fext);
    end
    
    nc.loadvars([],true);
    [lats,lons, ~,~, latbnds, lonbnds] = ncdf_get_latlons(fname);
    nlats = length(lats);
    nlons = length(lons);
    
    
    if (~isempty(latbnds))
%         latbnds = latbnds{1};
%         lonbnds = lonbnds{1};
        R = make_refmat([latbnds(1,1), latbnds(end,end)], [lonbnds(1,1),lonbnds(end,end)], nlats,nlons);
    else
        dlat = lats(2)-lats(1);
        dlon = lons(2)-lons(1);
        R = make_refmat(lats,lons, nlats, nlons, dlat, dlon);
    end

    vbose = nc.verbose(true);
    v=nc.readvar(varname);
    nc.verbose(vbose);
    
    gavg = squeeze(mean(v));
    if (nargout < 2), v=[]; end     % free up memory if user doesn't want v back.
    if (nanmean(gavg(:)) > 150)
        gavg=gavg - 273.15;
    end
    
    mymax = max(gavg(:));
    mymin = min(gavg(:));

    if (~exist('ming','var') || isempty(ming))
        ming = mymin;
    end
    if (~exist('maxg','var') || isempty(maxg))
        maxg = mymax;
    end
    
    gavg(gavg>maxg) = maxg;
    gavg(gavg<ming) = ming;
    cmap = cmap_CSC(255)
    display_map(gavg, ttl, R, fignum, subpos, ming, maxg, true, cmap);
    
    if (nargout == 0), gavg = []; end     % free up memory if user doesn't want gavg back.
    
end