% diary off;
  % try
%     QC_to_netcdf('Tmax','ncamerica',[1850,2017]);
%     !scp stations.Tmax.ncamerica.1850.2017.nc killarney:/data/obs/stations_netcdf
    QC_to_netcdf('Tmax','new_india',[1950,2017]);
    !scp stations.Tmax.new_india.1950.2017.nc killarney:/data/obs/stations_netcdf
% catch
% end
diary off;
% try
%     QC_to_netcdf('Tmin','ncamerica',[1850,2017]);
    QC_to_netcdf('Tmin','new_india',[1950,2017]);
    !scp stations.Tmin.new_india.1950.2017.nc killarney:/data/obs/stations_netcdf
% catch
% end
diary off;
% try
%     QC_to_netcdf('Prec','ncamerica',[1850,2017]);
    QC_to_netcdf('Prec','new_india',[1950,2017]);
    !scp stations.Prec.new_india.1950.2017.nc killarney:/data/obs/stations_netcdf
% catch
% end
diary off;