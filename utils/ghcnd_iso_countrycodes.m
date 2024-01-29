function ghcn_tbl=ghcnd_iso_countrycodes(ghcndir, isodir)        % countries in GHCN with slightly different names in ISO list which need manual matching.

% creates Country-Code/Region/Continent/Country-Name table from ghcnd_countries.txt file and iso country code files.
%
%   Inputs:  
%       ghcndir         directory where to find ghcn files, usually where you just downloaded latest ghcn data
%                           icsf computers:     [/Volumes/jcsf_data/ghcn/...]
%       isodir          directory where to find ISO-3166 country-code files 
%                           icsf computers:     [/Volumes/jcsf_data/ghcn/iso]
%                           CSC servers:        [/data/ghcn/iso]
%
%   latest GHCN data, files and documentation are available at https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/
%
%   iso 3166 country codes and continent info downloaded from:
%       https://dev.maxmind.com/static/csv/codes/iso3166.csv
%       https://dev.maxmind.com/geoip/legacy/codes/iso3166/iso3166.csv -> download as iso_country_codes.csv
%       https://dev.maxmind.com/geoip/legacy/codes/country_continent/country_continent.csv
%           note:  country/continent file needed a little manual editing to add some missing countries.
%
%       alternate, similar info/format:  https://github.com/lukes/ISO-3166-Countries-with-Regional-Codes
%           This git archive has more info, including UN-defined subregions.
%           This archive apparently combines iso country code info from Wikipedia and region info from UN.
%           **** I should implement the sub_regions and intermediate regions via this file, Ian!
%
%   This program creates 2 output tables:
%       regions.csv         With columns:
%                           region_code
%                               2- and 3-letter regions
%                                   2-letter are the iso continents
%                                   3-letter are TTU CSC regions, such as ncamerica, caribbean, etc.
%                                       these might be fully defined in github archive
%                           name    name of region (North America, Caribbean, etc.)
%
%       ghcn_country_region.csv
%                       table with the following columns 
%                               ghcn_cc         ghcn country code
%                               iso_cc          iso-3166 country code
%                               regions         comma-separated list of regions the country belongs in
%                               ghcn_country    ghcn country name
%                               iso_country     iso country name
%
%                       country names are slightly different between the lists, so update manual_matches as needed to
%                           join lists.
%                       Top part of list is GHCN countries and codes
%                       Bottom part of list is iso countries and codes not found in GHCN list.
%
%   Last updated:  2/26/2019 icsf
%____________________________________

    if (~exist('ghcndir','var') || isempty(ghcndir)), ghcndir="/Volumes/jcsf_data/ghcn/2019-02-13"; end
    if (~exist('isodir', 'var') || isempty(isodir)),  isodir ="/Volumes/jcsf_data/ghcn/iso"; end
    
    delim = extractBefore(ghcndir,2);
    if (delim ~= '/' && delim ~= '.'), ghcndir=fullfile("/Volumes/jcsf_data/ghcn",ghcndir); end
    ghcn_fname = fullfile(ghcndir,'ghcnd-countries.txt');
    if (~exist(ghcn_fname,'file'))
        error("error:  cannot open %s for reading. (use full path, or './...' to identify folder)",ghcn_fname);
    end
    
        % create the region table with region codes and region names.
    region_info= [
                    "AF","Africa"; ...
                    "AN","Antarctica"; ...
                    "AS","Asia"; ...
                    "EU","Europe"; ...
                    "NA","North America"; ...
                    "OC","Oceania"; ...
                    "SA","South America"; ...
                    "AUZ","Australia_NZ"; ...
                    "NCA","NCAmerica"; ...
                    "CAR","Carribean"; ...
                    ];
    region_tbl = table(region_info(:,1), region_info(:,2),'VariableNames',{'region_code','name'});
    writetable(region_tbl,fullfile(ghcndir,"regions.csv"),"QuoteStrings",true);
    
    manual_matches = [
        "Bahamas, The", "Bahamas", ""; ...
        "Brunei","Brunei Darussalam",""; ...
        "Burma", "Myanmar", ""; ...
        "Congo (Brazzaville)", "Congo", ""; ...
        "Congo (Kinshasa)", "Congo, The Democratic Republic of the", ""; ...
        "Europa Island [France]", "Europa Island","AF"
        "Falkland Islands (Islas Malvinas) [United Kingdom]" ,"Falkland Islands (Malvinas)",""; ...
        "Federated States of Micronesia","Micronesia, Federated States of",""; ...
        "French Southern and Antarctic Lands [France]","French Southern Territories",""; ...
        "Gambia, The","Gambia",""; ...
        "Iran","Iran, Islamic Republic of",""; ...
        "Jan Mayen [Norway]","Svalbard and Jan Mayen",""; ...
        "Johnston Atoll [United States]","","OC"; ...
        "Juan De Nova Island [France]","","AF"; ...
        "Korea, North","Korea, Democratic People's Republic of",""; ...
        "Korea, South","Korea, Republic of",""; ...
        "Laos","Lao People's Democratic Republic",""; ...
        "Palmyra Atoll [United States]","","OC"; ...
        "Libya","Libyan Arab Jamahiriya",""; ...
        "Macau S.A.R","Macao",""; ...
        "Moldova","Moldova, Republic of",""; ...
        "Midway Islands [United States}","","OC"; ...
        "Netherlands Antilles [Netherlands]","","NA"; ...
        "Pitcairn Islands [United Kingdom]","Pitcairn",""; ...
        "Russia","Russian Federation",""; ...
        "Svalbard [Norway]","Svalbard and Jan Mayen",""; ...
        "Syria","Syrian Arab Republic",""; ...
        "Tromelin Island [France]","","AF"; ...
        "Tanzania","Tanzania, United Republic of",""; ...
        "Virgin Islands [United States]","Virgin Islands, U.S.",""; ...
        "Wake Island [United States]","","OC";    
        ];

            % these are from Ranjini's code.  I think there should be more in each list...
            % 1st is region code, remainder are country codes to add the region code to.
            % ??_i  indicates iso country not in GHCN list.
    region_list{1} = ["NCA","CA","US","MX","BH","GT","ES","HO","NU","CS","PM"];       % ncamerica, north & central america.  Added El Salvador down to Panama.  icsf
    %region_list{2}= ["CAR","AA","AC",       "AV",       "BB","BD",       "BF",       "CJ","CU","DO","DR",       "GJ","GP","HA",       "JM",       "MB","MH",              "NN","RN","RQ","SC","ST","TB",       "TD","TK","UC", "VI","VC",         "VQ"]; % Caribbean list from Ranjini
    region_list{2} = ["CAR",     "AC","AI_i",     "AW_i","BB",     "BL_i","BF","BQ_i","CJ","CU","DO","DR","GD_i",     "GP",     "HT_i","JM","KN_i","MB",     "MF_i","MS_i","NN",     "RQ",     "ST",     "TC_i","TD",     "UC",      "VC_i","VG_i","VQ"]; % Caribbean, icsf list
    region_list{3} = ["AUZ","AS","KT","CK","HM_i","NZ","NF"];        % Australia & NZ
    region_list{4} = ["C_A","BH","CS","ES","GT","HO","NU","PM"];     % Central America


    % missing = false(218,1); 
    found=false;
    fid=fopen(fullfile(ghcndir,'ghcnd-countries.txt'));
    tline = fgetl(fid);
    i=0;
    gcc2=strings(0);
    while ischar(tline)
        i=i+1;
        ghcn_cc=strtrim(string(extractBefore(tline,' ')));    
        iso_cc="";
        regions="";
        cname=strtrim(string(extractAfter(tline,' ')));
                % check if it needs manual matching, or if there's a [...] to remove for matching.
        mm_ix = find(strcmpi(cname,manual_matches(:,1)));
        if (mm_ix)
            gcc2(i) = manual_matches(mm_ix,2);
            regions = manual_matches(mm_ix,3);
            if (strlength(regions)>0)
                found(i)=true; 
            end
        else
            kx=strfind(cname,'[');
            if (isempty(kx))
                gcc2(i) = "";
            else
                gcc2(i)=strtrim(extractBefore(cname,kx(1)));
            end
        end

        t=table(ghcn_cc,iso_cc,regions,cname,"",'VariableNames',{'ghcn_cc','iso_cc','regions','ghcn_country', 'iso_country'});
        if (i==1)
            ghcn_tbl=t;
        else
            ghcn_tbl=[ghcn_tbl;t]; %#ok<AGROW>
        end
        tline = fgetl(fid);
    end
    fclose(fid);
    nccs=size(ghcn_tbl,1);
    match_2 = zeros(nccs,1);
    found(end+1:nccs) = false;
    gcc2(end+1:nccs) = "";
    no_cont = false(nccs,1);
    multimatch=false(nccs);

    iso_tbl=readtable(fullfile(isodir,"iso_country_codes.csv"));
    iso_tbl.CC=strtrim(string(iso_tbl.CC));
    iso_tbl.Country=strtrim(string(iso_tbl.Country));
    niso = size(iso_tbl,1);
    iso_used= false(niso,1);

    region_tbl=readtable(fullfile(isodir,"iso_country_continent.csv"));
    region_tbl.Properties.VariableNames={'cc','continent'};
    region_tbl.cc=strtrim(string(region_tbl.cc));
    region_tbl.continent=strtrim(string(region_tbl.continent));

    for i=1:nccs
        ix=find(strcmpi(ghcn_tbl.ghcn_country(i),iso_tbl.Country),1);
        if (isempty(ix))
            if (strlength(gcc2(i))>0)
                ix=find(strcmpi(gcc2(i), iso_tbl.Country));
                if (length(ix)==1)
                    match_2(i) = ix;
                    found(i) = true;
                end
            end
        end

        if (isempty(ix))

        elseif (length(ix)>1)
            multimatch(i)=true;
            found(i)=false;
        else
            found(i)=true;
            iso_used(ix) = true;
            ghcn_tbl.iso_cc(i)=iso_tbl.CC(ix);
            ghcn_tbl.iso_country(i)=iso_tbl.Country(ix);
            jx = find(strcmpi(ghcn_tbl.iso_cc(i), region_tbl.cc),1);
            if (~isempty(jx))
                ghcn_tbl.regions(i)=region_tbl.continent(jx);
            else
                no_cont(i) = true;
            end
        end
    end

    for i=1:nccs
        if (no_cont(i))
            fprintf('no_continent_info: %3d %s %s\n', i, ghcn_tbl.cc(i),ghcn_tbl.country(i));
        end
    end
    for i=1:nccs
        if (multimatch(i))
            fprintf('multiple_matches:  %3d %s %s\n', i, ghcn_tbl.cc(i),ghcn_tbl.country(i));
        end
    end
    for i=1:nccs
        if (match_2(i))
            fprintf('match_2:           %3d ghcn:  "%s" "%s" ("%s") \t\tiso: "%s" "%s"\n', i, ghcn_tbl.ghcn_cc(i), ghcn_tbl.ghcn_country(i), gcc2(i), iso_tbl.CC(match_2(i)), iso_tbl.Country(match_2(i)));
        end
    end
    fprintf('\n');
    for i=1:nccs
        if (~found(i))
            fprintf('missing:           %3d %s %s\n', i, ghcn_tbl.ghcn_cc(i),ghcn_tbl.ghcn_country(i));
        end
    end
    fprintf('\n');

    for i=1:niso
        if (any(strcmp(iso_tbl.CC(i),["A1","A2","01","O1","--"]))), continue; end
        if (~iso_used(i))
            fprintf('unused_iso: %3d %s %s\n', i, iso_tbl.CC(i),iso_tbl.Country(i));
            jx = find(strcmpi(iso_tbl.CC(i), region_tbl.cc),1);
            if (~isempty(jx))
                region=region_tbl.continent(jx);
            else
                region="";
            end
            t=table("",iso_tbl.CC(i),region,"",iso_tbl.Country(i),'VariableNames',{'ghcn_cc','iso_cc','regions','ghcn_country', 'iso_country'});
            ghcn_tbl=[ghcn_tbl; t]; %#ok<AGROW>
        end
    end

        % append additional region strings

    for i=1:length(region_list)
        region=region_list{i}(1);
        for j=2:length(region_list{i})
            cc=extractBefore(region_list{i}(j),3);
            if (strlength(region_list{i}(j))==2)
                ix=find(strcmpi(cc,ghcn_tbl.ghcn_cc));
            else
                ix=find(strcmpi(cc,ghcn_tbl.iso_cc));
            end
            if (isempty(ix))
                fprintf("oops:  no match for adding region %s to %s\n", region,region_list{i}(j));
            elseif (length(ix)>1)
                fprintf("oops:  multiple matches for adding region %s to %s\n", region,region_list{i}(j));
            else
                if (isempty(ghcn_tbl.regions(ix)))
                    ghcn_tbl.regions(ix)=region;
                else
                    ghcn_tbl.regions(ix)=sprintf("%s,%s", ghcn_tbl.regions(ix),region);
                end
            end
        end
    end


    fprintf('\n');

    writetable(ghcn_tbl,fullfile(ghcndir,"ghcn_country_region.csv"),"QuoteStrings",true);

end
