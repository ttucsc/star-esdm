function cc_codes = GHCN_country_codes(regions)
%
% returns string array of GHCN country codes for various regions.
%
%   Currently defined regions:
%         "AF","Africa"; ...
%         "AN","Antarctica"; ...
%         "AS","Asia"; ...
%         "EU","Europe"; ...
%         "NA","North America"; ...
%         "OC","Oceania"; ...
%         "SA","South America"; ...
%         "AUZ","Australia_NZ"; ...
%         "NCA","NCAmerica"; ...
%         "CAR","Carribean"; ...
%         "ncamerica" -- NA + C_A + CAR
%----------------------------

    regions = string(regions);
    rcodes=strings(0,0);
    for i=1:length(regions)
        rcodes=[rcodes;region_codes(regions(i))];
    end
    
    cc_codes = strings(0,0);
    
    cc_tbl = readtable("/Volumes/jcsf_data/ghcn/2019-02-13/ghcn_country_region.csv");
    for r=1:length(rcodes)
        ix = find(contains(cc_tbl.regions,rcodes(r)));
        for i=1:length(ix)
            if (strlength(cc_tbl.ghcn_cc(ix(i)))>0)
                cc_codes(end+1,1) =  cc_tbl.ghcn_cc(ix(i)); %#ok<AGROW>
            end
        end
    end
    cc_codes = unique(sort(cc_codes));
end

function rcodes = region_codes(region)

    if (strcmpi(region,'ncamerica'))
        rcodes = ["NA";"NCA";"CAR"];
    else
        rcodes = region;
    end
end
    