    function s = llgrid_string(lat,lon, lat2, lon2, isdelta)
        % isdelta:  true if lat2 & lon2 are gridbox size, rather than ending lat & lon.
    
        if (lat < 0)
            latc='S';
        else
            latc='N';
        end
        
        if (mod(lat,1)==0)
            s=sprintf('%c%02d',latc,abs(lat));
        else
            s=sprintf('%c%05.2f',latc,abs(lat));
        end

        lon=mod(lon+180,360)-180;
        if (lon<0)
            lonc='W';
        else
            lonc='E';
        end
        if (mod(lon,1)==0)
            s=sprintf('%s%c%03d',s,lonc,abs(lon));
        else
            s=sprintf('%s%c%06.2f',s,lonc,abs(lon));
        end
        
        if (exist('lat2','var') && ~isempty(lat2))
            
             if (lat2 < 0)
                latc='S';
            else
                latc='N';
             end
           if (mod(lat2,1)==0)
               if (~isdelta)
                    s=sprintf('%s_%c%02d',s,latc,abs(lat2));
               else
                    s=sprintf('%s_%02d',s,abs(lat2));
               end
           else
                if (~isdelta)
                    s=sprintf('%s_%c%.2f',s,latc,abs(lat2));
                else
                    s=sprintf('%s_%.2f',s,abs(lat2));
                end                    
           end
            
            lon2=mod(lon2+180,360)-180;
            if (lon2<0)
                lonc='W';
            else
                lonc='E';
            end
            if (mod(lon2,1)==0)
                if (~isdelta)
                    s=sprintf('%s_%c%03d',s,lonc,abs(lon2));
                else
                    s=sprintf('%sx%03d',s,abs(lon2));
                end
            else
                if (~isdelta)
                    s=sprintf('%s%c%.2f',s,lonc,abs(lon2));
                else
                    s=sprintf('%sx%.2f',s,abs(lon2));
                end
            end
        end
    end
