function G = ghcn_station_info(fname, do_sort,do_table)

    % reads GHCN station info file (usually ghcnd-stations.txt), returns a struct w/ all fields.
    % Matlab's tables were terribly slow for large amounts of data, so default is struct.
    % (Matlab seems to have made them faster... icsf, 2019, so added do_table as option.)
    
    % inputs:
    %   fname       filename
    %   do_sort     if true, sort by stn_id before returning.   default: false
    %   do_table    if true, return table, rather than struct.  default: false
    
    % outputs:
    % struct with fields
    %   stn_id
    %   lat
    %   lon
    %   elev
    %   state
    %   stn_name
    %   gsn         bool flag, true if data is from GSN database
    %   hcn         bool, true if from HCN
    %   wmo_id
    
    % note:  code does not identify whether station is COCORAHS precip station, which could be deteremined from the
    % stn_id.  There are a huge number of COCORAHS precip stations, which only have precip data;  COCORAHs started in
    % 1998.
            
    if (~exist('do_sort', 'var') || isempty_s(do_sort)),  do_sort  = false; end
    if (~exist('do_table','var') || isempty_s(do_table)), do_table = false; end

    fid = fopen(fname);
    s=textscan(fid,'%11s%*1c%8f%*1c%9f%*1c%6f%*1c%2s%*1c%30s%*1c%3s%*1c%3s%*1c%5s','WhiteSpace','');
    fclose(fid);

    stn_id = string(strtrim(s{1}));
    lat=s{2};
    lon=s{3};
    elev=s{4};
    state=string(strtrim(s{5}));
    stn_name=string(strtrim(s{6}));
    gsn1 = string(strtrim(s{7}));
    hcn_crn=string(strtrim(s{8}));
    gsn= gsn1=='GSN';
    hcn = hcn_crn=="HCN";
    crn = hcn_crn=="CRN";
    wmo_id=string(strtrim(s{9}));
    
        % make all lons between -180.0 & 180.0;
    ix=lon>180.0;
    lon(ix)=lon(ix)-360.0;
    
    missing_lats=find(isnan(lat));
    missing_lons=find(isnan(lon));
    
    if (~isempty(missing_lats) || ~isempty(missing_lons))
        fprintf('%s:  missing lats:  \n', fname);
        for jj=1:length(missing_lats)
            j=missing_lats(jj);
            fprintf('%6d %-20s\t%.4f %.4f\n', j, stn_id(j), lat(j),lon(j));
        end

        fprintf('%s: missing lons:  \n', fname);
        for jj=1:length(missing_lons)
            j=missing_lons(jj);
            fprintf('%6d %-20s\t%.4f %.4f\n', j, stn_id(j), lat(j),lon(j));
        end
    end    
    
    G = struct('stn_id',stn_id, ...
               'lat',lat, ...
               'lon',lon, ...
               'elev',elev, ...
                'state',state, ...
                'stn_name',stn_name, ...
                'gsn',gsn, ...
                'hcn',hcn, ...
                'crn', crn, ...
                'wmo_id',wmo_id);        

    
    
     if (do_sort)
        [~,ix] = sort(stn_id);
        fields=fieldnames(G);
        for i=1:length(fields)
            f=fields{i};
            G.(f) = G.(f)(ix);
        end

%         [stn_id,ix] = sort(stn_id);
%         stn_name = stn_name(ix);
%         lat = lat(ix);
%         lon = lon(ix);
%         state = state(ix);
%         gsn = gsn(ix);
%         hcn = hcn(ix);
%         crn = crn(ix);
%         wmo_id = wmo_id(ix);
     end  
    
     if (do_table)
         G = struct2table(G);
     end
     
     missing = G.elev==-999.9;
     G.elev(missing) = nan;
    
end
