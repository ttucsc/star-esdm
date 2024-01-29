function arrm_put_globals( ncid, info )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    % Insert global attributes

    % These are just my stuff, you need to customize. But note attributes related to the convention. 
    % This is very useful for you netcdf file to read by anybody else than you.

    %   this is needed because putAtt won't take a string as the Varid.
    %   getConstant returns -1 for NC_GLOBAL, but better not to hard-code -1 here in case it changes in the future.
    NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');
    
%    if (~is_netcdf4),  netcdf.reDef(ncid); end;
    

    if (~isempty_s(info.title)),          netcdf.putAtt(ncid, NC_GLOBAL, 'title',             info.title);       end
    if (~isempty_s(info.long_title)),     netcdf.putAtt(ncid, NC_GLOBAL, 'long_title',        info.long_title);  end
    if (~isempty_s(info.comments)),       netcdf.putAtt(ncid, NC_GLOBAL, 'comments',          info.comments);    end
    if (~isempty_s(info.institution)),    netcdf.putAtt(ncid, NC_GLOBAL, 'institution',       info.institution); end
    if (~isempty_s(info.progname)),       netcdf.putAtt(ncid, NC_GLOBAL, 'source',            info.progname);    end
    if (~isempty_s(info.data_source)),    netcdf.putAtt(ncid, NC_GLOBAL, 'data_source',       info.data_source); end
    if (~isempty_s(info.history)),        netcdf.putAtt(ncid, NC_GLOBAL, 'history',           info.history);     end
    if (~isempty_s(info.references)),     netcdf.putAtt(ncid, NC_GLOBAL, 'references',        info.references);  end
    
    netcdf.putAtt(ncid, NC_GLOBAL, 'Conventions', 'CF-1.6');      % orig. was 1.4
    % This way I can always get back to fundamentals:
    netcdf.putAtt(ncid, NC_GLOBAL, 'Conventions_help', 'http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html'); 
    netcdf.putAtt(ncid, NC_GLOBAL, 'CreationDate', datestr(info.tstamp,'yyyy/mm/dd HH:MM:SS'));
    netcdf.putAtt(ncid, NC_GLOBAL, 'CreatedBy', getenv('LOGNAME'));
    
end

