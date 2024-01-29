function [sunrise,sunset,solarnoon, sr, ss, sn] = solar_info(dates, latitude, longitude_east, utcoff, secs)
%function [sunrise,sunset,solarnoon, sr, ss, sn] = solar_info(dates, latitude, longitude_east, utcoff, secs)
%   Returns sunrise, sunset & solarnoon for date(s) @ lat, lon.
%       lat, lon default to lubbock, tx
%
%   Inputs:
%       dates       [ yyyy, mm, dd] or datenum;  can be multiple rows for
%                                   multiple dates
%       lat, lon    decimal or [dd hh mm].  Note that lon is EAST of GMT
%                                   (i.e., use neg for NAmerica)
%       nsteps      # of timesteps (sets precision of sunrise, sunset, noon)
%                       288 = every 5 minutes; 1440 = every minute
%       utcoff      # # of hours between current loc & GMT (neg for US)
%
%       secs        accuracy desired (in secs) [30]
%
%   Outputs:
%       sunrise, sunset, solarnoon   times, in decimal hours
%       sr, ss, sd                   times in [hh mm ss]
%
%   Uses:  suncycle.m from http://mooring.ucsd.edu/software/matlab/doc/toolbox/geo/suncycle.html
%
%           note:  suncycle is not a good tool for this...doesn't appear to calculate
%           using  refraction, and generates a table of sun positions.
%           I find max value in solar zenith to get transit time.
%
%           I should use equations from NREL's "Solar Position Algorithm for
%           Solar Radiation Applications", or some other source, but didn't
%           have time to code them up.
%
%               ian
%
%---------------------------


    if (~exist('latitude','var') || isempty(latitude))
        latitude=33.578;
    end
    if (~exist('longitude_east','var') || isempty(longitude_east))
        longitude_east=-101.855;
    end
    if (~exist('secs','var') || isempty(secs))
        secs=30;     % every 5 minutes
    end
    
    if (~exist('utcoff','var') || isempty(utcoff))
        utcoff=-6;
    end
    
    if (size(dates,2)==1)       % input is datenums already
        dnums=dates;
    else
        dnums=datenum(dates);      % or not.  convert dates to datenums
    end
    ndnums = length(dnums);
    nsites = size(latitude,1);

    sunrise     = nan(nsites, ndnums,1);
    sunset      = nan(nsites, ndnums,1);
    solarnoon   = nan(nsites, ndnums,1);
    sr          = nan(nsites, ndnums,3);
    ss          = nan(nsites, ndnums,3);
    sn          = nan(nsites, ndnums,3);

    if (size(latitude,2)>1)     % convert [d,m,s] to deg
        signs=sign(latitude);
        latitude = abs(latitude);
        latitude=signs.*(latitude(:,1)+latitude(:,2)/60.0+latitude(:,3)/60.0/60.0);
    end
    if (size(longitude_east,2)>1)   % convert [d,m,s] to deg
        signs=sign(longitude_east);
        longitude_east = abs(longitude_east);
        longitude_east = signs .* (longitude_east(:,1)+longitude_east(:,2)/60.0+longitude_east(:,3)/60.0/60.0);
    end

    nsteps=floor(24*3600/secs);
    
    for s_ix=1:nsites
        
    %   [rs,t,d,z,a,r]=suncycle(latitude,longitude_east,dn,nsteps);
        [rs,t,~,z,~,~]=suncycle(latitude(s_ix),longitude_east(s_ix),dnums,nsteps);

        rs=mod(rs+24+utcoff(s_ix),24);

        for i=1:length(dnums)
            sunrise(s_ix,i)=rs(i,1);
            sunset(s_ix,i)=rs(i,2);
            [sr(s_ix,i,1),sr(s_ix,i,2),sr(s_ix,i,3)]=hms(rs(i,1));
            [ss(s_ix,i,1),ss(s_ix,i,2),ss(s_ix,i,3)]=hms(rs(i,2));
            ix=find(z(i,:)==max(z(i,:)));
    %        snd(i)=24*(ix-1)/nsteps;
            solarnoon(s_ix,i)=mod(mod(t(i,ix),1)*24+24+utcoff(s_ix),24); %#ok<FNDSB>
            [sn(s_ix,i,1),sn(s_ix,i,2),sn(s_ix,i,3)]=hms(solarnoon(s_ix,i));
    %         dv=datevec(dn(i));
    %         fprintf('%4d/%02d/%02d: %02d:%02d:%02d\t%02d:%02d:%02d\t%02d:%02d:%02d\n',dv(1),dv(2),dv(3),sr(i,1:3),sn(i,1:3),ss(i,1:3));
        end
    end

        % squeeze out the singleton dimension if only 1 lat/lon given.
    if (nsites == 1)
        sunrise     = squeeze(sunrise);
        sunset      = squeeze(sunset);
        solarnoon   = squeeze(solarnoon);
        sr          = squeeze(sr);
        ss          = squeeze(ss);
        sn          = squeeze(sn);
    end
        
end
    
function [h, m, s]=hms(t)

    sgn = sign(t);
    t = abs(t);
    h=sgn * floor(t);
    m=floor(mod(t*60,60));
    s=floor(mod(t*3600,60));
end
