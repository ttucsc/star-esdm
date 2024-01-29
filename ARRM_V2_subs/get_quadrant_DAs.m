function [mdl_DAs_quadrant, hist_DAs_quadrant, DSP_mdl, insufficient_data] = get_quadrant_DAs(mdllat_ix, mdllon_ix, mdllats, mdllons, all_mdldata, all_histdata, RP, DSP_base , DSP_mdls, DSP_hists)
    % could do away with this, just do weighted sum and create a single DA from it.


        % create DAs for the four gridcells we're using at the moment,
        % calculate the binning, and histogram the DAs.
        %   We have to do this each time because binning will be different for each set of four gridpoints.
    mdl_DAs_quadrant  = cell(2,2);
    hist_DAs_quadrant = cell(2,2);
    insufficient_data = false;
    DSP_mdl  = DSP_mdls.update( 'file_lats', mdllats(mdllat_ix:mdllat_ix+1), 'file_lons',mdllons(mdllon_ix:mdllon_ix+1));
    DSP_hist = DSP_hists.update('file_lats', mdllats(mdllat_ix:mdllat_ix+1), 'file_lons',mdllons(mdllon_ix:mdllon_ix+1));    
    
    for i=1:2
        for j=1:2

            ix=mdllat_ix+i-1;
            jx=mdllon_ix+j-1;
            mdl_data = all_mdldata(:,ix,jx);
%           DSP_base.print_log("mdl  ");
            mdl_DAs_quadrant{i,j} = ARRM_V2_DisaggregateSignal(mdl_data, RP, DSP_mdl,sprintf('DA_run_%s_%02d_%02d',DSP_base.llgrid_lbl, ix,jx), false, false, false, false, false, false);  % don't do any disaggregation;  we'll do the other steps later. 
            if (mdl_DAs_quadrant{i,j}.insufficient_data)
                DSP_base.warn_log("Insufficient Model data:  (%8.4f, %8.5f)\n", mdllats(ix),mdllons(jx));
                insufficient_data = true;
                break;
            end

            hist_data = all_histdata(:,ix,jx);
%           DSP_base.print_log("hist ");
            hist_DAs_quadrant{i,j} = ARRM_V2_DisaggregateSignal(hist_data, RP, DSP_hist,sprintf('DA_run_%s_%02d_%02d',DSP_base.llgrid_lbl, ix,jx), false, false, false, false, false, false);  % don't do any disaggregation;  we'll do the other steps later. 
            if (hist_DAs_quadrant{i,j}.insufficient_data)
                DSP_base.warn_log("Insufficient Hist data:  (%8.4f, %8.5f)\n", mdllats(ix),mdllons(jx));
                insufficient_data = true;
                break;
            end
        end
    end
        
end
 
