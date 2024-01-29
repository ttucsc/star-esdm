function elev_temp(varname, e_or_c)

    season_lbl = ["winter", "spring", "summer", "fall", "full-year"];
    if (~exist('varname','var'))
        varname="tasmax";
        e_or_c = "cm3";
    end

%         hist = load(matname_hist, "varname", "e_or_c", "mean_fig", "trend_fig", "clim_mu", "clim_sig", "mclim_mu",  "mclim_sig", "prob_rms_figs", "prob_mean_figs", "probs","lats1","lons1", "latslow", "lonslow");
%         fut  = load(matname_fut,  "varname", "e_or_c", "mean_fig", "trend_fig", "clim_mu", "clim_sig", "mclim_mu",  "mclim_sig", "prob_rms_figs", "prob_mean_figs", "probs","lats1","lons1", "latslow", "lonslow");

    syr_fut  = 2086;
    eyr_fut  = 2115;
    syr_hist = 1979;
    eyr_hist = 2038;
    
    is_tasmax = strcmp(varname,"tasmax");
    is_cm3    = strcmp(e_or_c, "cm3");
    
    if (strcmp(e_or_c,"cm3"))
        fname_hi_fut  = sprintf("/Volumes/lacie_1/projects/ian/work/future/c360_hiram_cm3_rcp85_X1X3X4.hires.%s.20860101-21151231.llt.nc", varname);
        fname_hi_hist = sprintf("/Volumes/lacie_1/projects/ian/work/hist/%s_day_GFDL-HIRAM-C360_amip_19790101-20381231.H2H3.hires.llt.nc", varname);

        fname_low_fut  = sprintf("/Volumes/lacie_1/projects/ian/work/future/c360_hiram_cm3_rcp85_X1X3X4.lowres.%s.20860101-21151231.llt.nc", varname);
        fname_low_hist = sprintf("/Volumes/lacie_1/projects/ian/work/hist/%s_day_GFDL-HIRAM-C360_amip_19790101-20381231.H2H3.lowres.llt.nc", varname);
    else
        fname_hi_fut  = sprintf("/Volumes/lacie_1/projects/ian/work/future/c360_hiram_esm_rcp85_X1X2X3.hires.%s.20860101-21151231.llt.nc", varname);
        fname_hi_hist = sprintf("/Volumes/lacie_1/projects/ian/work/hist/%s_day_GFDL-HIRAM-C360_amip_19790101-20381231.H2H3.hires.llt.nc", varname);

        fname_low_fut  = sprintf("/Volumes/lacie_1/projects/ian/work/future/c360_hiram_esm_rcp85_X1X2X3.lowres.%s.20860101-21151231.llt.nc", varname);
        fname_low_hist = sprintf("/Volumes/lacie_1/projects/ian/work/hist/%s_day_GFDL-HIRAM-C360_amip_19790101-20381231.H2H3.lowres.llt.nc", varname);
    end 
    
        % get the difference between the downscaled future data and the hires futue data from the post-downscaling disaggregation
    DA_matname = sprintf("GFDL_PM.%s.%s.%d.%d.mat", varname, e_or_c, syr_fut, eyr_fut);
    DA_mat = load(DA_matname,"mean_fig");
    DA_mean_diff = - DA_mat.mean_fig;



    if (~isfile(fname_hi_fut)),   error("missing future hires  file: %s\n", fname_hi_fut); end
    if (~isfile(fname_hi_hist)),  error("missing hist   hires  file: %s\n", fname_hi_hist); end
    if (~isfile(fname_low_fut)),  error("missing future lowres file: %s\n", fname_low_fut); end
    if (~isfile(fname_low_hist)), error("missing hist   lowres file: %s\n", fname_low_hist); end

    latrange = [20,50];
    lonrange = [-125,-60];

    [clim_hi_fut,     lats,    lons] = get_mean_clim(fname_hi_fut,   varname, latrange, lonrange, 30, false);
    [clim_hi_hist,       ~,       ~] = get_mean_clim(fname_hi_hist,  varname, latrange, lonrange, 60, false);
    [clim_low_fut, latslow, lonslow] = get_mean_clim(fname_low_fut,  varname, latrange, lonrange, 30, true);
    [clim_low_hist,      ~,       ~] = get_mean_clim(fname_low_hist, varname, latrange, lonrange, 60, true);
    
    nlats_hi  = size(clim_hi_fut, 2);
    nlons_hi  = size(clim_hi_fut, 3);
%   nlats_low = size(clim_low_fut,2);
%   nlons_low = size(clim_low_fut,3);
    

    clim_delta_hi  = clim_hi_fut  - clim_hi_hist;
    clim_delta_low = clim_low_fut - clim_low_hist;

%       [elevs, avg_elevs, diff_elevs] = conus_elevations("conus_elevs.mat");

%     dmax = max([max(clim_hi_fut(:)), max(clim_hi_hist(:)), max(clim_low_fut(:)), max(clim_low_hist(:))]);
%     dmin = min([max(clim_hi_fut(:)), min(clim_hi_hist(:)), min(clim_low_fut(:)), min(clim_low_hist(:))]);
%     ddmax = max(max(clim_delta_hi(:)), max(clim_delta_low(:)));
%     ddmin = min(min(clim_delta_hi(:)), min(clim_delta_low(:)));

    dmin = -25;
    dmax =  50;

    ddmin = 0;
    ddmax = 12;
    
    delmin = -4;
    delmax = 5;

    R_hi  = georefcells(latrange, lonrange, [nlats_hi,  nlons_hi]);
%   R_low = georefcells(latrange, lonrange, [nlats_low, nlons_low]);
    
    cmap = jet(32);
    cmap2 = jet(18);
    
    cmap(15:17,:) = [.3333,1,1;.8,1,.8;1,1,1];
    cmap2(8:10,:) = [.3333,1,1;.9,1,1;1,1,1];

                % get elevation data for scatterplot of difference vs elevation.
%     e = load("conus_elevs.mat");
%     elevs = e.elevs;
%     diff_elevs = e.diff_elevs;

    subpos = [1,4,7,10,13; 2, 5, 8, 11, 14; 3, 6, 9, 12, 15];
    
    figbase = 500 + 10*is_cm3 + 20*is_tasmax;
    for issn = 1:5

        fut_hi      = squeeze(mean(get_season(clim_hi_fut,    issn)));
        hist_hi     = squeeze(mean(get_season(clim_hi_hist,   issn)));
        fut_low1    = squeeze(mean(get_season(clim_low_fut,   issn)));
        hist_low1   = squeeze(mean(get_season(clim_low_hist,  issn)));
        delta_hi    = squeeze(mean(get_season(clim_delta_hi,  issn)));
        delta_low1  = squeeze(mean(get_season(clim_delta_low, issn)));
        DA_diff     = squeeze(DA_mean_diff(issn,:,:));
        
        fut_low   = expand(fut_low1,   lats, lons, latslow, lonslow);
        hist_low  = expand(hist_low1,  lats, lons, latslow, lonslow);
        delta_low = expand(delta_low1, lats, lons, latslow, lonslow);
    
        

        
        [del_nearest, del_quadrant] = diff_temps(delta_hi, delta_low1, lats, lons, latslow, lonslow);
        DA_del = DA_diff - del_quadrant;

        fignum = figbase + mod(issn,5);
        h = figure(fignum);

        display_map(fignum, [3,2,1], fut_hi,   sprintf("%s %s %s hires Mean %d-%d",  season_lbl(issn), varname, e_or_c, syr_fut, eyr_fut),    dmin, dmax,  cmap, R_hi, latslow, lonslow);
        display_map(fignum, [3,2,3], hist_hi,  sprintf("%s %s %s hires Mean %d-%d",  season_lbl(issn), varname, e_or_c, syr_hist, eyr_hist),  dmin, dmax,  cmap, R_hi, latslow, lonslow);
        display_map(fignum, [3,2,5], delta_hi, sprintf("%s %s %s hires Mean Change", season_lbl(issn), varname, e_or_c),                     ddmin, ddmax, cmap, R_hi, latslow, lonslow);

        display_map(fignum, [3,2,2], fut_low,   sprintf("%s %s %s lowres Mean %d-%d",  season_lbl(issn), varname, e_or_c, syr_fut, eyr_fut),    dmin, dmax,  cmap, R_hi, latslow, lonslow);
        display_map(fignum, [3,2,4], hist_low,  sprintf("%s %s %s lowres Mean %d-%d",  season_lbl(issn), varname, e_or_c, syr_hist, eyr_hist),  dmin, dmax,  cmap, R_hi, latslow, lonslow);
        display_map(fignum, [3,2,6], delta_low, sprintf("%s %s %s lowres Mean Change", season_lbl(issn), varname, e_or_c),                     ddmin, ddmax, cmap, R_hi, latslow, lonslow);

        figname = sprintf("GFDL_pm_tempmap.%s.%s.%s.png", varname, e_or_c, season_lbl(issn));
        saveas(h, figname, "png");  
        
        display_map(figbase+5, [3,2,issn], del_nearest,  sprintf("%s %s %s hi - Nearest(low) ", season_lbl(issn), varname, e_or_c),  delmin, delmax, cmap2, R_hi, latslow, lonslow);
        display_map(figbase+6, [5,3,subpos(1,issn)], del_quadrant, sprintf("%s %s %s hi - Qudrant(low) ", season_lbl(issn), varname, e_or_c),  delmin, delmax, cmap2, R_hi, latslow, lonslow);
        display_map(figbase+6, [5,3,subpos(2,issn)], DA_diff,      sprintf("%s %s %s hi - Downscaled ", season_lbl(issn), varname, e_or_c),  delmin, delmax, cmap2, R_hi, latslow, lonslow);
        display_map(figbase+6, [5,3,subpos(3,issn)], DA_del,       sprintf("%s %s %s Downscaled difference ", season_lbl(issn), varname, e_or_c),  delmin, delmax, cmap2, R_hi, latslow, lonslow);
        
%       scatter_elevs(issn,figbase+7, del_quadrant, elevs, diff_elevs, season_lbl(issn), varname, e_or_c);
            
%       draw_elevs();
    end
    h = figure(figbase+5);
    figname = sprintf("GFDL_pm_tempmap_trend.%s.%s.%s.png", varname, e_or_c, "nearest");
    saveas(h, figname, "png");        
    h = figure(figbase+6);
    figname = sprintf("GFDL_pm_tempmap_trend.%s.%s.%s.png", varname, e_or_c, "quadrant");
    saveas(h, figname, "png");        
%     h = figure(figbase+7);
%     figname = sprintf("GFDL_pm_trend_vs_elev.%s.%s.%s.png", varname, e_or_c, "quadrant");
%     saveas(h, figname, "png");        
end

function data = get_season(data, issn)
    season_days={[334:365,1:60]; 61:151; 152:242; 243:333; 1:365};

    if (issn < 1 || issn > 4), return; end
    keepers = season_days{issn};
    data = data(keepers, :,:);
    
end

function display_map(fignum, pos, img, ttl, ming, maxg, cmap, R, latslow, lonslow)
    img = squeeze(img);
    myming = min(img(:));
    mymaxg = max(img(:));
    if (~exist('ming','var') || isempty(ming)), ming = myming; end
    if (~exist('maxg','var') || isempty(maxg)), maxg = mymaxg; end
    
    img(img<ming) = ming;
    img(img>maxg) = maxg;
    if (min(img(:)) > ming)
        img(1,1)=ming;
    end
    if (max(img(:)) < maxg)
        img(end,end) = maxg;
    end

    fig=figure(fignum);
    if (~isempty(pos)), subplot(pos(1), pos(2), pos(3)); end
    worldmap(img, R);
    
    meshm(img,R);
    ncamerica_and_state_boundaries([],true);
    if (exist('latslow','var'))

        for i=1:length(latslow)
            plotm(latslow(i)*ones(size(lonslow)), lonslow, '-','color',[.75,.75,.75]);
        end
        for i=1:length(lonslow)
            plotm(latslow, lonslow(i)*ones(size(latslow)), '-','color',[.75,.75,.75]);
        end
    end
    hold off;

%   imshow(flipud(img),[ming,maxg],'colormap',cmap);
    if (myming ~= ming || mymaxg ~= maxg)
        title(sprintf("%s range %.2f-%.2f", ttl, myming, mymaxg));
    else
        title(ttl);
    end
    colorbar('location','eastoutside');
    colormap(cmap);
    if (pos(2)==2)
        fig.Position = [200, 100, 1400, 1300];
    else
        fig.Position = [200, 100, 1800, 1300];
    end
    drawnow();
    pause(.25)

end

% function draw_elevs()
%     
%     e = load("conus_elevs.mat");
%     elevs = e.elevs;
%     avg_elevs = e.avg_elevs;
%     diff_elevs = e.diff_elevs;
% %     lats = e.lats;
% %     lons = e.lons;
%     latslow = e.latslow;
%     lonslow = e.lonslow;
%     
%     minhi = min(elevs(:));
%     maxhi = max(elevs(:));
%     minlow = min(avg_elevs(:));
%     maxlow = max(avg_elevs(:));
%     
%     mindiff = min(diff_elevs(:));
%     maxdiff = max(diff_elevs(:));
%     
%     dmax = max(max(abs(elevs(:))), max(abs(avg_elevs(:)))); 
%     elevs(1,1) = -dmax;
%     elevs(end,end) = dmax;
%     avg_elevs(1,1) = -dmax;
%     avg_elevs(end,end) = dmax;
%     
% %    dmax = max(abs(diff_elevs(:)));
%     dmax = 1500;
%     diff_elevs(1,1) = -dmax;
%     diff_elevs(end,end) = dmax;
%     diff_elevs(diff_elevs > dmax) = dmax;
%     diff_elevs(diff_elevs < -dmax) = dmax;
%     
%     R=georefcells([20,50],[-125,-60],[120,208]);
% 
%     subplot(3,3,3);
%     worldmap(elevs, R);
%     meshm(elevs,R);
%     ncamerica_and_state_boundaries([],true);
%     lowgrid(latslow, lonslow, [20,50],[-125,-60]);
%     colorbar("location","eastoutside");
%     title(sprintf("Elevations (m)  range: %.0f-%.0f", minhi, maxhi));
%     hold off;
%     
% 
%     subplot(3,3,6);
%     worldmap(elevs, R);
%     meshm(avg_elevs,R);
%     ncamerica_and_state_boundaries([],true);
%     lowgrid(latslow, lonslow, [20,50],[-125,-60]);
%     colorbar("location","eastoutside");
%     title(sprintf("Mean Elevation for gridcells (m)  range:  %.0f-%.0f", minlow, maxlow));
%     hold off;
%     
%     subplot(3,1,9);
%     worldmap(elevs, R);
%     meshm(diff_elevs,R);
%     ncamerica_and_state_boundaries([],true);
%     lowgrid(latslow, lonslow, [20,50],[-125,-60]);
%     colorbar("location","eastoutside");
%     title(sprintf("Elevation difference (m)  range:  %.0f-%.0f", mindiff, maxdiff));
%     hold off;
%     
%     colormap(jet(16));
%     
% end

function [clim, lats, lons] = get_mean_clim(fname,   varname, latrange, lonrange, nyrs, is_low)


    nc = ncdf(fname);
    nc.loadvars([], true);
    try
        lats = nc.getvardata("lat");
        lons = nc.getvardata("lon");
    catch
        lats = nc.getvardata("latitude");
        lons = nc.getvardata("longitude");
    end
    
    if (is_low)     % make sure we bracking the full range if this is for low-resolution file.
        dlat = lats(10)-lats(9);
        dlon = lons(10) - lons(9);
        latrange = latrange + [-dlat, dlat];
        lonrange = lonrange + [-dlon, +dlon];
    end
    
    latkeepers = lats >= latrange(1) & lats <= latrange(2);
    lats = lats(latkeepers);
    ilat1 = find(latkeepers,1);
    nlats = length(lats);
    
    lons = mod(lons+180, 360)-180;
    lonkeepers = lons >= lonrange(1) & lons <= lonrange(2);
    lons = lons(lonkeepers);
    ilon1 = find(lonkeepers, 1);
    nlons = length(lons);
    
    start = [1, ilat1, ilon1];
    count = [nyrs*365, nlats, nlons];
    
    d = nc.readvar(varname, start, count) - 273.15;
    
    d = reshape(d, 365, nyrs, nlats, nlons);
    
    clim = squeeze(mean(d,2));
end

% function lowgrid(latslow, lonslow, latrange, lonrange)
% 
%     lonslow = mod(lonslow+180, 360) - 180;
%     keepers = latslow >= latrange(1) & latslow <= latrange(2);
%     latslow = latslow(keepers);
%     keepers = lonslow >= lonrange(1) & lonslow <= lonrange(2);
%     lonslow = lonslow(keepers);
%     
%     for i=1:length(latslow)
%         plotm(latslow(i)*ones(size(lonslow)), lonslow, '-','color',[.75,.75,.75]);
%     end
%     for i=1:length(lonslow)
%         plotm(latslow, lonslow(i)*ones(size(latslow)), '-','color',[.75,.75,.75]);
%     end
% end    

function [del_nearest, del_quadrant] = diff_temps(delta_clim_hi, delta_clim_low, lats, lons, latslow, lonslow)

    [nr,nc] = size(delta_clim_hi);
    del_nearest = nan(nr,nc);
    del_quadrant = nan(nr,nc);
    for ilat=1:nr
        for ilon=1:nc
            [quad_latix, quad_lonix, near_latix, near_lonix] = find_ix(ilat, ilon, lats, lons, latslow, lonslow);           
            dl = delta_clim_low(quad_latix, quad_lonix);
            del_quadrant(ilat,ilon) = delta_clim_hi(ilat,ilon) - mean(dl(:));
            del_nearest(ilat,ilon) = delta_clim_hi(ilat,ilon) - delta_clim_low(near_latix, near_lonix);
        end
    end
end

function [quad_latix, quad_lonix, near_latix, near_lonix] = find_ix(ilat, ilon, lats, lons, latslow, lonslow)

    lat=lats(ilat);
    quad_latix = nan(1,2);
    quad_lonix = nan(1,2);
    
    ql = find(latslow < lat,1, 'last');
    if (isempty(ql))
        quad_latix(1) = 1; 
    elseif (ql == length(latslow))
        quad_latix(1) = ql-1; 
    else
        quad_latix(1) = ql;
    end
    quad_latix(2) = quad_latix(1)+1;
    
    if (abs(lat-latslow(quad_latix(1))) < abs(latslow(quad_latix(2))-lat))
        near_latix = quad_latix(1);
    else
        near_latix = quad_latix(2);
    end

    lon = lons(ilon);
    ql = find(lonslow < lon, 1, 'last');
    if (isempty(ql))
        quad_lonix(1) = 1; 
    elseif (ql == length(lonslow))
        quad_lonix(1) = ql-1; 
    else
        quad_lonix(1) = ql;
    end
    quad_lonix(2) = quad_lonix(1)+1;
    
    if (abs(lon-lonslow(quad_lonix(1))) < abs(lonslow(quad_lonix(2))-lon))
        near_lonix = quad_lonix(1);
    else
        near_lonix = quad_lonix(2);
    end                        
end

function [clim_quadrant, clim_nearest]  = expand(clim_low,  lats, lons, latslow, lonslow)

    nr = length(lats);
    nc = length(lons);
    clim_nearest  = nan(nr,nc);
    clim_quadrant = nan(nr,nc);
    
    for ilat = 1:nr
        for ilon = 1:nc
            [qlatix, qlonix, llatix, llonix] = find_ix(ilat, ilon, lats, lons, latslow, lonslow);
            clim_nearest(ilat, ilon) = clim_low(llatix, llonix);
            lc = clim_low(qlatix, qlonix);
            clim_quadrant(ilat,ilon) = interp_bilinear_ic(lats(ilat), lons(ilon), latslow(qlatix), lonslow(qlonix), lc(1,1), lc(2,1), lc(1,2), lc(2,2));
%           clim_quadrant(ilat, ilon) = mean(lc(:));
%           wts = bilinear_weights_ic(lats(ilat), lons(ilon), latslow(qlatix(1)), latslow(qlatix(2)), lonslow(qlonix(1)), lonslow(qlonix(2)));
        end
    end
end

% function         scatter_elevs(issn,fignum, del_quadrant, elevs, elev_diffs, season_lbl, varname, e_or_c)
% 
%     subpos = [1,2,5,6,9; 3,4,7,8,11];
%     fig = figure(fignum);
%     subplot(3,4,subpos(1,issn));
%     scatter(elevs(:), del_quadrant(:), 5, 'b','filled','o');
%     xlabel("elevation");
%     ylabel("trend difference, deg C");
%     ylim([-5,5]);
%     grid on;
%     title(sprintf("Hires - Lowres trend vs Elevation, %s %s %s", season_lbl, varname, e_or_c));
% 
%     subplot(3,4,subpos(2,issn));
%     scatter(elev_diffs(:), del_quadrant(:), 5,'b','filled','o');
%     xlabel("elevation");
%     ylabel("trend difference, deg C");
%     ylim([-5,5]);
%     grid on;
%     title(sprintf("Hires - Lowres trend vs Elevation Change, %s %s %s", season_lbl, varname, e_or_c));
%     fig.Position = [200 50, 1400, 1400 ];
%     pause(.25);
% 
% 
% end
