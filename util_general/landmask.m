function outmask = landmask(latix, lonix, lsmask, masklats,masklons)
    % returns lsmask values for location specified by (lat,lon) index pairs or
    % (lat,lon) pairs
    % If masklats & masklons provided, then we assume latix and lonix are
    % lat/lon values, not indexes.
    %
    %   For an interpolated land/sea mask, use landsea_mask(...)
    %   For mask values for a set of lat/lon values interpolated from a
    %   different gridding, use is_land(...)
    
    if (exist("masklats","var"))    % if given masklats & masklons, then latix and lonix are actually 
        inlats = latix;             % lat/lon pairs we expect to find exact matches for in
        inlons = lonix;             % masklats and masklons.
        for i=1:length(inlats)                       
            latix(i) = find(masklats == inlats(i));
            lonix(i) = find(masklons == inlons(i));
            if (numel(latix)~=1 || numel(lonix) ~= 1)
                error("problem matching ( %.4f , %.4f ) in masklats & masklons", inlats(i), inlons(i))
            end
        end
    end
    outmask = lsmask(to_column(latix),to_column(lonix));
end
