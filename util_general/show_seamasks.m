[~,hostname]=system('hostname');
onlaptop = contains(hostname, "icsf-jpro");
do_conv = true;

if (isfolder("/Volumes/lacie_1"))
    mydir = "/Volumes/lacie_1/killarney_data/gcm_cmip5/land_sea_mask/";
    [dinf,fnames] = dirmatch(fullfile(mydir, "sftlf*.nc"));
    nfn = numel(fnames);
    fnames = ["/Volumes/lacie_1/killarney_data/gcm_cmip5/hires/gfdl/c360_hiram_H2/N360_to_M45/atmos.gridinfo.N360_to_M45.nc";
              "/Volumes/lacie_1/killarney_data/gcm_cmip5/hires/gfdl/c360_hiram_H2/N360_amip/landfrac.N360_amip.nc";
              fnames];
    varnames=["land_mask"; "land_frac"; repmat("sftlf",nfn,1) ];
else
    mydir = "land_frac";
    [dinf,fnames] = dirmatch(fullfile(mydir, "*.nc"));
    bnames = basename(fnames);
    nfn = sum(strncmp(bnames,"sftlf",5));
    varnames=["land_frac"; "land_mask"; repmat("sftlf",nfn,1) ];
end
      

load ('coastlines','coastlat','coastlon');

% cmap1 = summer(8);
% cmap  = cmap1; % [cmap1(:,1),cmap1(:,3), cmap1(:,2)];
% for i=1:8
%     cmap(i,2) = .5*(1+cmap(i,2)); 
% end
% 

cmap = mycmap;

for i=1:length(fnames)
    
%   try
        nc = ncdf(fnames(i));
        [~,fn,ext] = fileparts(fnames(i));
        fname = sprintf('%s%s', fn, ext);
    %   nc.loadvars([],true);
    
        if (strcmp(fname,"sftlf_fx_BNU-ESM_historical_r0i0p0.nc"))
            lat_bnds = linspace(-87.8638,87.8638,64);
            lon_bnds = nc.readvar('lon_bnds');
            nlats=max(size(lat_bnds));
            nlons=max(size(lon_bnds));
            fprintf("%s:  generated lat_bnds.\n\tsize:  %d x %d\n", fname, nlats, nlons);
        else
            try
                lat_bnds=nc.readvar('lat_bnds');
                lon_bnds=nc.readvar('lon_bnds');
                nlats=max(size(lat_bnds));
                nlons=max(size(lon_bnds));
                fprintf("%s:\n\tsize:  %d x %d\n", fname, nlats,nlons);
            catch
                lats = nc.readvar('lat');
                lons = nc.readvar('lon');

                dlat=(lats(4)-lats(3))/2;
                dlon=(lons(4)-lons(3))/2;
                nlats=length(lats);
                nlons=length(lons);
                lat_bnds = round([to_row(lats)-dlat;to_row(lats)+dlat],4);
                lon_bnds = round([to_row(lons)-dlon;to_row(lons)+dlon],4);
                fprintf("%s:  using lat,lon.  \n\tsize:  %d x %d\n", fname, max(size(lat_bnds)),max(size(lon_bnds)));


            end

            latdifs = unique(round(diff(lat_bnds(1,:)),4));
            londifs = unique(round(diff(lon_bnds(1,:)),4));
            if (length(latdifs) > 1 || length(londifs) > 1)
                fprintf("warning:  unevenly spaced lat_bnds or lon_bnds: %s\n", fnames(i));
                fprintf("\tlatdifs:  "); 
                disp(latdifs);
                fprintf("\tlondifs:  ");
                disp(londifs);
            end
        end

        ttl = fnames(i);

%       d=nc.readvar('sftlf');
        d=nc.readvar(varnames(i));
        d_orig = d;
        land = d_orig == 1;
        if (do_conv)
            convsize = round(.025*nlats);
            krnl = ones(convsize,convsize)/(convsize*convsize);
            d = conv2(d, krnl,"same");
            mx = max(d(:));
            nz = d>0;
            d(nz) = .75*mx + .25*d(nz);
            d(land) = mx;   % reset land to 1
        end

        d=(d')/max(d(:));
        try
            R = mk_refmat(lat_bnds, lon_bnds);
        catch
            oops();
        end

        h=figure(i);
        clf;
        if (onlaptop)
            figpos = [5*i, 25-2*i, 800,650 ];
        else
            figpos = [200+5*i, 300-2*i, 2000,1200 ];
        end
        h.Position = figpos;
        axesm ('miller', 'Frame', 'on', 'Grid', 'on');
        geoshow(d, R, 'DisplayType', 'texturemap');
        plotm(coastlat, coastlon,'-','color','k','linewidth',1); %[1,.5,.25]
        title(sprintf("%s  size %d x %d", fname, nlats, nlons),'interpreter','none');
        colormap(cmap);
%     catch
%         fprintf("problem with %d %s\n", i, fnames(i));
%     end
%   display_map(d, ttl, R, i, [], 0, 100, 1, [], -90:30:90,0:30:360, figpos)
end

function R = mk_refmat(lat_bnds, lon_bnds)

    nlats = length(lat_bnds);
    nlons = length(lon_bnds);
    
    latrange = [lat_bnds(1,1),lat_bnds(end,end)];
    lonrange = [lon_bnds(1,1),lon_bnds(end,end)];
        
    R = georefcells(latrange, lonrange, [nlats, nlons]);
    
end

function cmap = mycmap
    cmap = zeros(8,3);
    cmap(:,1) = linspace(0,1,8);
    cmap(:,2) = linspace(.5,.4,8);
    cmap(:,3) = linspace(.75,1,8);
    cmap(8,:) = [1,1,.4];
end