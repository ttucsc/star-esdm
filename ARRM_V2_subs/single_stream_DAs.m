function [DA_mdl, DA_hist] = single_stream_DAs(mdl_DAs_quadrant, hist_DAs_quadrant, weights, run_lbl, gridpt_loc)

        % create matrix with the mdl_DAs raw data
    [nr,nc] = size(mdl_DAs_quadrant);
    mdl_DAs_quadrant = reshape(mdl_DAs_quadrant,1, nr*nc);
    DSP_mdl  = mdl_DAs_quadrant{1}.DP;
    RP       = mdl_DAs_quadrant{1}.RP;
    DSP_mdl.gridpt_loc = gridpt_loc;
    RP.weights = 1;
    
    mdl_data = weighted_sum(mdl_DAs_quadrant, weights, "raw_data");
    
    if (DSP_mdl.isPrecipRun)
        DA_mdl = ARRM_V2_DisaggregateSignal(mdl_data, RP, DSP_mdl, run_lbl, false, false, false, false, false, false);  % dont do anomalies or calc binning;  we'll do that later. 

        DSP_hist = hist_DAs_quadrant{1}.DP;
        DSP_hist.gridpt_loc = gridpt_loc;

        hist_data = weighted_sum(hist_DAs_quadrant, weights, "raw_data");

        DA_hist = ARRM_V2_DisaggregateSignal(hist_data, RP, DSP_hist, run_lbl, false, false, false, false, false, false);  % dont do anomalies or calc binning;  we'll do that later. . 
    else

        DA_mdl = ARRM_V2_DisaggregateSignal(mdl_data, RP, DSP_mdl, run_lbl, true, false, false, false, false, false);  % do anomalies & calc binning;  we'll do the other steps later. 
        DA_mdl.DP.gridpt_loc = gridpt_loc;
    
        minedge = DA_mdl.RP.edges(1);
        maxedge = DA_mdl.RP.edges(end);
        binstep = DA_mdl.RP.binstep;
        nbins   = length(DA_mdl.RP.bins);
    %   DSP_mdl.print_log("model gridcell binning:  maxedge:  %22.18f   minedge: %22.18f  binstep:  %22.18f   range:  %22.18f  nbins: %d  \n", maxedge, minedge, binstep, maxedge-minedge, nbins);
        DSP_mdl.print_log("model gridcell binning:  maxedge:  %9.5f   minedge: %9.5f  binstep:  %9.5f   range:  %9.5f  nbins: %d  \n", maxedge, minedge, binstep, maxedge-minedge, nbins);

        DSP_hist = hist_DAs_quadrant{1}.DP;
        DSP_hist.gridpt_loc = gridpt_loc;

        hist_data = weighted_sum(hist_DAs_quadrant, weights, "raw_data");

        DA_hist = ARRM_V2_DisaggregateSignal(hist_data, RP, DSP_hist, run_lbl, true, false, false, false, false, false);  % do anomalies & calc binning;  we'll do the other steps later. 

        minedgeh = DA_hist.RP.edges(1);
        maxedgeh = DA_hist.RP.edges(end);
        binsteph = DA_hist.RP.binstep;
        nbinsh   = length(DA_hist.RP.bins);
%       DSP_mdl.print_log("hist  gridcell binning:  maxedge:  %22.18f   minedge: %22.18f  binstep:  %22.18f   range:  %22.18f  nbins: %d  \n", maxedgeh, minedgeh, binsteph, maxedgeh-minedgeh, nbinsh);
        DSP_mdl.print_log("hist  gridcell binning:  maxedge:  %9.5f   minedge: %9.5f  binstep:  %9.5f   range:  %9.5f  nbins: %d  \n", maxedgeh, minedgeh, binsteph, maxedgeh-minedgeh, nbinsh);
    end
end

