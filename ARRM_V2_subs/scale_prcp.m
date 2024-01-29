function DA = scale_prcp(DA, fldname, direction, prcp_scale_info)

%   scale data forward or reverse using prcp_scale_info.
%   This function naively assumes that all dry days have been replaced with NAs, so none of the arithmetic here will affect them.

        % to scale precip data:
        %   1.  subtract prcp_min, & threshold at 0  (i.e., any thing less tha 0 --> 0)
        %   2.  raise to prcp_scaling power
        %   3.  add dx/2.
        % to scale back to regular precip values:
        %   1.  subtract dx/2
        %   2.  raise to 1/prcp_scaling power
        %   3.  add prcp_min.

        % scale_info = struct("prcp_scaling",prcp_scaling, "dx",dx,"sigma",prcp_sigma,"edges",edges, "prcp_min", prcp_min, "nbins", nbins, "nsigmas", nsigmas);


%   prcp_scale info is pr_min on forward scaling, or struct with scaling info on reverse scaling.

%   max_scaling = 10;
%     nbins       = 1000;
%     nsigmas     = 10;

    if (strcmp(direction, "forward"))
        
        prcp_min    = prcp_scale_info.prcp_min;
        scaling     = prcp_scale_info.prcp_scaling;
        prcp_offset   = prcp_scale_info.prcp_offset;

        prcp = DA.(fldname);

        scaled_prcp = rescale_data_1(prcp, scaling, prcp_min, "forward", prcp_offset);

        DA.(fldname) = scaled_prcp;
        
        DA.prcp_scaling = prcp_scale_info;       % save the scaling info in the DA so we can reverse the scaling later.
    elseif (strcmp(direction, "reverse"))
        prcp_min    = prcp_scale_info.prcp_min;
        scaling     = prcp_scale_info.prcp_scaling;
        prcp_offset   = prcp_scale_info.prcp_offset;

        scaled_prcp = DA.(fldname);
        
        prcp = rescale_data_1(scaled_prcp, scaling, prcp_min, "reverse", prcp_offset);
        
        DA.(fldname) = prcp;
        
%       DA.(fldname) = ARRM_V2_precip_scaling(DA.(fldname), scaling, pr_min, "reverse", dx);     
    else
        error("bad direction info:  %s\n", direction);
    end

end

