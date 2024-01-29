
yr0=2065;
nyrs=30;
ncfut_hi = "/Volumes/lacie_1/data/gcm_cmip5/hires/ian_work/future/c360_hiram_cm3_rcp85_X1X3X4.hires.pr.20660101-20951231.llt.nc";
prcp_lbb_hi=ncdf_read_nearest_single_location(ncfut_hi,"pr",33.5667, 258.1167,[yr0+1,yr0+nyrs],"365-day");
units = ncdf_getvar_info(ncfut_hi,"pr");
prcp_lbb_hi = jc_units_conversion(prcp_lbb_hi, units,"mm");  % convert from kg/m/sec to mm/day
nyrs = length(prcp_lbb_hi)/365;
yrs=yr0+1+(1:length(prcp_lbb_hi))/365;
for j=1:nyrs
    
    mu  = zeros(nyrs,1);
    sig = zeros(nyrs,1);
    jx1=(j-1)*365+1;
    jx2= j*365;
    yr1 = prcp_lbb_hi(jx1:jx2);  
    
    h=figure(30+j);  
    h.Position=[225 648 2033 677];
    subplot(2,1,1)
    plot(yrs, prcp_lbb_hi, yrs(jx1:jx2), yr1); 
    grid on;

    for i=1:nyrs
        yr2 = prcp_lbb_hi(((i-1)*365)+1:(i*365));
        mu(i)  = mean(yr1-yr2);
        sig(i) = std(yr1-yr2);
    end
    subplot(2,1,2);
    plot(yr0+(1:nyrs), mu, yr0+(1:nyrs), sig);
    legend("mean","sig");
    grid on;
    mu(j)=nan;
    sig(j)=nan;
    k1 = find(abs(mu)==min(abs(mu)),1,"omitnan");
    k2 = find(sig==min(sig),1,"omitnan");
    title(sprintf("%d %.5g %.4f %.4f", yr0+j-1, mu(k1), sig(k2), max(prcp_lbb_hi)))
    drawnow();
end

% nchist_hi = "/Volumes/lacie_1/data/gcm_cmip5/hires/ian_work/hist/pr_day_GFDL-HIRAM-C360_amip_19790101-20381231.H2H3.hires.llt.nc";
% prcp_lbb_hi=ncdf_read_nearest_single_location(nchist_hi,"pr",33.5667, 258.1167,[1979,2038],"365-day");
% units = ncdf_getvar_info(nchist_hi,"pr");
% prcp_lbb_hi = jc_units_conversion(prcp_lbb_hi, units,"mm");  % convert from kg/m/sec to mm/day
% nyrs = length(prcp_lbb_hi)/365;
% yrs=1979+(1:length(prcp_lbb_hi))/365;
% 
% for j=1:60
%     
%     mu  = zeros(nyrs,1);
%     sig = zeros(nyrs,1);
%     jx1=(j-1)*365+1;
%     jx2= j*365;
%     yr1 = prcp_lbb_hi(jx1:jx2);  
%     
%     h=figure(j);  
%     h.Position=[225 648 2033 677];
%     subplot(2,1,1)
%     plot(yrs, prcp_lbb_hi, yrs(jx1:jx2), yr1); 
%     grid on;
% 
%     for i=1:nyrs
%         yr2 = prcp_lbb_hi(((i-1)*365)+1:(i*365));
%         mu(i)  = mean(yr1-yr2);
%         sig(i) = std(yr1-yr2);
%     end
%     subplot(2,1,2);
%     plot(1978+(1:60), mu, 1978+(1:60), sig);
%     legend("mean","sig");
%     grid on;
%     mu(j)=nan;
%     sig(j)=nan;
%     k1 = find(abs(mu)==nanmin(abs(mu)),1);
%     k2 = find(sig==nanmin(sig),1);
%     subplot(2,1,1); 
%     title(sprintf("%d %.5g %.4f", 1979+j-1, mu(k1), sig(k2)))
%     drawnow();
% end
