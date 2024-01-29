function  fignums = anomaly_distributions(varargin)
%   anomaly_distributions(model_list, stn_list_filename, ['argname',argvalue...])
%
%   Inputs:
%       models          list of models to include.  Will match any string, so "GFDL" matches all models with string GFDL
%
%   Optional parameters
%       stn_fname       csv filename with lat & lon columns of locations to use for anomalies
%                               default:  "test_stations_1000.csv"  
%
%   Optional Name/Value pairs
%       'exclude_sites', exclude_sites  array of sites to specifically exclude []  (use "+" for 1st to include-only.
%                                           sites are stnNames or stnIDs from stn_fname file  (matches if name contains)
%       'base_yrs', base_yrs            years to use for anomalies [1950,2005]
%       'trend_yrs', trend_yrs          years to use for trendline [1935,2035]
%       'probs', probs                  probability values to calculate equations for
%       'ncname', ncname                (will create from varname, model, scenario, ensemble, etc. if not provided)
%       'varname', varname              variable, scenario & ensemble   ["tasmax"]  use "stations" for station probs
%       'scenario', scenario                                            ["rcp85"]
%       'ensemble', ensemble                                            ["r1i1p1"]
%       'sitenums', [...]               list of station numbers to run for model. [all]
% other optional ones
%       'do_calc_binning'
%       'fig_base'          150  (use 'figbase' to pass figure # on to Disaggregation code)
%       'figpos'            []
%       'figname'           [generated]
%       'anomfile'          [generated] name of output file for anomalies info
%       'read_anomfile'     [false] if true, reads anomfile instead of generating it.
%       'do_plots'          [true]
%       'partial_match'     [true] if true, does partial matching on model name 

    rp = initparams(varargin{:});
%   if (strcmp(rp.base_dir,"help") || strcmp(rp.base_dir,"-h")), show_help(); return; end
    
    fignums = [];
    if (isempty(rp.sitenums))
        nsites = rp.nsites;
        msites = size(rp.sites,1);
        ss_step = floor(msites/nsites);
        sitenums = 1:ss_step:msites;
    else
        sitenums = rp.sitenums;
        nsites=length(rp.sitenums);
        msites = nsites;
    end
    anomfile = rp.anomfile;
    if (rp.read_anomfile)
        if (~exist(anomfile,'file'))
            fprintf("anomaly_distributions:  anomfile %s does not exist\n", anomfile);
        else
            load(anomfile,'all_hist2','all_bins','rp','all_anoms');
            npts = sum(~isnan(all_anoms(:))); 
            save_mdl_cdf(all_bins, all_hist2,rp,  npts)
            
        end        
    else

        dname = fullfile(rp.basedir,sprintf('log/anomaly_distributions_%s_%s_%04d_%04d_log.txt',rp.models(1), rp.varname, msites, nsites));
        if (exist(dname,'file')), delete(dname); end
        diary(char(dname));

        fignum = rp.fig_base;
        figname = rp.figname;
        fignums = [];

        fprintf('\n\n------Anomaly Distributions  %s----------\n\n',datestr(now));
    
        nsets = length(rp.ncnames);
        for ifn=1:nsets

            ncnames = fullfile(rp.base_dirs{ifn}, rp.ncnames{ifn});

            missing = false;
            for i=1:length(ncnames)               
                if (~exist(ncnames(i),'file')), fprintf('file %s doesn''t exist\n', ncnames(i)); missing = true; end
            end
            if (missing), continue; end

            nyrs = rp.base_yrs(2)-rp.base_yrs(1)+1;
            mpts = 365*nyrs*nsites;

            cdfs = cell(nsites,1);
            plines = cell(nsites,1);
            delta = .025;
            all_edges = -15:delta:15;
            all_bins = delta/2+all_edges(1:end-1);
            nbins = length(all_bins);
            all_hist1 = zeros(365,nbins);
            all_hist2 = zeros(365,nbins);
            all_anoms = zeros(mpts,1);
            mbad = 0;
            ix1 = 0;

            for sss=1:length(sitenums)
                ss = sitenums(sss);

                [RP, DP, stnName, rp] = setup_run(ifn, ss, rp, ncnames, nsets, msites);

                if (isempty(RP)), continue; end

                DA = ARRM_V2_DA_run('RP',RP,'DP',DP);
                DA.report_prob_fails(DA.fail_info_base,   'extreme outliers, base period'   ); 
                DA.report_prob_fails(DA.fail_info_rolling,'extreme outliers, rolling period'); 

                [bstart,bend] = DA.using_range("base_yrs");
                my_cdf = DA.CDF_base;
                problines = DA.problines_base;
                anoms = DA.anoms(bstart:bend);

                doy = repmat((1:365)',nyrs,1);
                ix2=ix1 + length(anoms);        
                all_anoms((ix1+1):ix2) = anoms;
                ix1 = ix2;

                nbad = sum(isnan(anoms));
                mbad = mbad + nbad;
                cdfs{i} = my_cdf;
                plines{i} = problines;

                pl_minus2 = problines(:,2);     % -2.5 sigma
                pl_minus1 = problines(:,3);     % -1 sigma
                pl_mus   = problines(:,4);    % mean
                pl_plus1  = problines(:,5);    % +1 sigma
                pl_plus2  = problines(:,6);    % +2.5 sigma

%                 msig = mean(pl_plus - pl_minus);
                anoms = reshape(anoms,365,nyrs);
                anoms1 = zeros(size(anoms));
                anoms2 = zeros(size(anoms));

                for j=1:365
                    anoms1(j,:) =  (anoms(j,:)-pl_mus(j)) * 2/(pl_plus1(j)-pl_minus1(j));  % centered & normalized to 1 sigma.
                    anoms2(j,:) =  (anoms(j,:)-pl_mus(j)) * 5/(pl_plus2(j)-pl_minus2(j));  % centered & normalized to 1 sigma.
                end
                anoms1 = anoms1(:);
                anoms2 = anoms2(:);
                my_hcounts1 = histcounts2(doy,anoms1,.5:366,all_edges);
                my_hcounts2 = histcounts2(doy,anoms2,.5:366,all_edges);
                all_hist1 = all_hist1 + my_hcounts1; 
                all_hist2 = all_hist2 + my_hcounts2; 
                fprintf('nan count: %3d %-45s  %6d\n', ss, stnName, nbad);
            end

            if (~isempty(anomfile))
                save(rp.anomfile,'all_anoms','all_hist1','all_hist2','all_bins','rp');
            end
            
            if (rp.do_plots)
                fignums = plot_results(all_hist1, all_hist2, all_bins, fignum, figname,rp.models(1), rp.varname, fignums, rp.basedir);
            end

        end
        
    end

    diary off;
end

function save_mdl_cdf(mdl_bins, all_hist2, rp, npts)
%   saves basic cdf, bins & okflags to mat file in subfolder "cdf".
%   

    if (length(rp.models)>1)
        outname = fullfile(rp.basedir,'cdf',sprintf('cdf_basic.%s.%s.mat', "multimodel", rp.varname));
    else
        outname = fullfile(rp.basedir,'cdf',sprintf('cdf_basic.%s.%s.mat', rp.models, rp.varname));
    end
    
%   nbins = length(mdl_bins);
    tot_hist = sum(all_hist2,'omitnan');
    mdl_cdf = cumsum(tot_hist,'omitnan');
    mdl_cdf = mdl_cdf/mdl_cdf(end);
    
        % flag where to use the CDF.
    mdl_pdf = diff([0,mdl_cdf]);        % create pdf from cdf to avoid computer rounding near 1.  Safer than using (tot_hist~=0).
    okflags = mdl_pdf~=0;               % all places where CDF will be increasing
    ix1 = find(okflags,1) - 1;          % last zero-point
    okflags(max(1,ix1)) = true;         %#ok<NASGU> % exclude the last point so we don't interpolate all the way to 1 later.
%   ix1 = find(okflags,1,'last');       % last 1-point
    
%    okflags(ix1) = false;               %#ok<NASGU> % exclude the last point so we don't interpolate all the way to 1 later.
    
    [mu, sigma, skewness, kurtosis] = pdf_stats(mdl_pdf',mdl_bins');
    mdl_stats = struct('mu',mu,'sigma',sigma,'skewness',skewness,'excess_kursosis',kurtosis,'npts',npts); %#ok<NASGU>
    
    models = rp.models; %#ok<NASGU>
    varname = rp.varname; %#ok<NASGU>
    fprintf("saving cdf file %s\n", outname);
    save(outname,"mdl_cdf","mdl_bins","okflags", "mdl_stats","models", "varname");
end    


function fignums = plot_results(all_hist_1, all_hist_2, all_bins, fignum, figname, model, varname, fignums, basedir)

    % _1 is normalized to 1 sigma, _2 to 2.5.  _2 is generally better, so displaying that mainly.
    % need to calc slopes and save results here. still.
    
    all_counts = cumsum(all_hist_2,2);
    pdf_all_2 = all_hist_2  ./ all_counts(:,end);
    cdf_all_2 = all_counts ./ all_counts(:,end);
    zs = pdf_all_2 < 1e-15;
    cdf_all_2(zs) = nan;
    all_counts = cumsum(all_hist_1,2);
    pdf_all_1 = all_hist_1  ./ all_counts(:,end);
    cdf_all_1 = all_counts ./ all_counts(:,end);
    zs = pdf_all_1 < 1e-15;
    cdf_all_1(zs) = nan;
    
        % --------------------------------- 1st figure
        
    fignums(end+1) = fignum;
    h = figure(fignum);
    clf;
    fignum = fignum+1;
    
    
    subplot(2,2,1);
    surf(all_bins-5, 1:365, cdf_all_1,'edgecolor','none')
    hold on;
    surf(all_bins+5, 1:365, cdf_all_2,'edgecolor','none')
    hold off;
    lights_on;
    ylabel('day of year');
    xlabel('anomaly');
    zlabel('probability');
    xlim([-15,15]);
    title(sprintf('%s %s CDF, 1- & 2.5-sig normalized', model, varname), 'interpreter','none')
    grid on;
    
    
    subplot(2,2,2);
    ah_1 = all_hist_1;
    ah_1(ah_1==0) = nan;
    surf(all_bins+5, 1:365, ah_1,'edgecolor','none')
    hold on
    ah_2 = all_hist_2;
    ah_2(ah_2==0) = nan;
    surf(all_bins-5, 1:365, ah_2,'edgecolor','none')
    hold on
    
    lights_on;
    ylabel('day of year');
    xlabel('anomaly');
    zlabel('count');
    xlim([-15,15]);
    grid on;
    
    dx = all_bins(2)-all_bins(1);
    ndist = normpdf(all_bins) * dx;
    tot_hist_2 = sum(all_hist_2);
    flags_2 = tot_hist_2 ~= 0;
    nn = sum(tot_hist_2, 'omitnan');
    my_pdf_2 = tot_hist_2 / nn;
    my_pdf_2 = ic_kde(my_pdf_2, true);
    
    subplot(2,2,3);
    plot(all_bins(flags_2), ndist(flags_2), 'g','linewidth',3);
    hold on;
    plot(all_bins(flags_2), my_pdf_2(flags_2),'b','linewidth',1.5);
    hold off;
    grid on;
    title(sprintf('%s %s, normalized to 2.5-sig_equiv', model, varname), 'interpreter','none');
    legend('normal',sprintf('%s %s',model,varname));
    
%    my_cdf_2 = cumsum(tot_hist_2,'omitnan') / nansum(tot_hist_2);
    my_cdf_2 = cumsum(my_pdf_2,'omitnan');
    
    dx = all_bins(2)-all_bins(1);
    ndist = normpdf(all_bins) * dx;
    tot_hist_1 = sum(all_hist_1);
    flags_1 = tot_hist_1 ~= 0;
    my_pdf_1 = tot_hist_1 / sum(tot_hist_1,'omitnan');
    my_pdf_1 = ic_kde(my_pdf_1, true);
    mnmx1 = xlim();
    xlabel('sigmas');
    ylabel('pdf');
    
    subplot(2,2,4);
    plot(all_bins(flags_1), ndist(flags_1), 'g','linewidth',3);
    hold on;
    plot(all_bins(flags_1), my_pdf_1(flags_1),'b','linewidth',1.5);
    hold off;
    grid on;
    title('normalized to 1-sig_equiv', 'interpreter','none');
    mnmx2 = xlim();
    xlabel('sigmas');
    ylabel('pdf');
    
    ylims = ylim();
    ylims(1) = -.025*ylims(2);
    ylim(ylims);
    xlim([min([mnmx1,mnmx2,-7.5]),max([mnmx1,mnmx2,7.5])]);
    subplot(2,2,3);
    xlim([min([mnmx1,mnmx2,-7.5]),max([mnmx1,mnmx2,7.5])]);
    ylim(ylims);
    
%     my_cdf_1 = cumsum(tot_hist_1,'omitnan') / nansum(tot_hist_1);
    
    if (~isempty(figname))
        h.Position=[100,300,1000,1000];
        fn = fullfile(basedir,sprintf("%s_%04d.tif", figname, fignum));
        saveas(h,fn,'tif');
    end
    
        % --------------------------------- 2nd figure
        
    
    keepers = tot_hist_2 > 0;
    nkeepers = sum(keepers);
    my_cdf_2 = my_cdf_2(keepers);
    ok_bins   = all_bins(keepers);
    sigs = -3:3;
    psigs = normcdf(sigs);
    siglocs = interp1(my_cdf_2,ok_bins,psigs, 'makima',nan);
    ctrix = find(sigs==0);
    ctr = siglocs(ctrix); 
    negbins = ok_bins < ctr;
    posbins = ok_bins >= ctr;
    siglocs(ctrix)=[];
    sigs(ctrix)=[];
    msigs = length(siglocs)+1;
    nsigs = length(siglocs)/2;
    
    ratios = zeros(nkeepers,msigs-1);
    
    ymin = 10^floor(min(log10(my_cdf_2)));

    fignums(end+1) = fignum;
    h = figure(fignum);
    fignum = fignum+1;
    
    xlims1 = [nan,nan];
    xlims2 = [nan,nan];
    for i=1:nsigs
        j=msigs-i;
        ratios(negbins,i) = ok_bins(negbins) / siglocs(i) * sigs(j);
        ratios(posbins,j) = ok_bins(posbins) / siglocs(j) * sigs(j);
        normprobs = normcdf(ok_bins);      
        
        subplot(2,nsigs,i);
        h0 = plot(ok_bins(negbins), normprobs(negbins),'g','linewidth', 6);
        hold on;
        plot(ok_bins(posbins), 1-normprobs(posbins),'g','linewidth', 6);
        h1 = scatter(-ratios(negbins,i),   my_cdf_2(negbins),15,'b','filled');        
        scatter(-ratios(posbins,j), 1-my_cdf_2(posbins),3,'r','filled');      
        h2 = scatter(ratios(posbins,j), 1-my_cdf_2(posbins),15,'r','filled');
        scatter(ratios(negbins,i),   my_cdf_2(negbins),3,'b','filled');        
        hold off;
        set(gca,'yscale','log');
        if (i==1)
            title(sprintf('%s %s prob vs scld %d sigs', model, varname, sigs(j))); 
        else
            title(sprintf('probability vs scaled %d sigs', sigs(j)));
        end
        xlabel('scaled sigmas');
        ylabel('cum prob (left), 1-cum prob (right)');
        legend([h0, h1,h2], 'normal','left tail','right tail');
        grid on;
        ylim([ymin,inf]);
        xlims1 = [min([xlim(),xlims1]),max([xlim(),xlims1])];
        
        
        subplot(2,nsigs,i+nsigs);
        h0 = plot(ok_bins(negbins), normprobs(negbins),'g','linewidth', 5);
        hold on;
        plot(ok_bins(posbins), 1-normprobs(posbins),'g','linewidth', 5);
        h1 = scatter(ratios(negbins,i),   my_cdf_2(negbins), 10,'b','filled');
        h2 = scatter(ratios(posbins,j), 1-my_cdf_2(posbins), 10,'r','filled');
        hold off;
        set(gca,'yscale','log');
        xlabel('scaled sigmas');
        ylabel('cum prob (left), 1-cum_prob (right)');
        legend([h0,h1,h2], 'normal','left tail','right_tail');   
        title(sprintf('probability vs scaled %d sigs', sigs(j))); 
        grid on;
        ylim([ymin,inf]);
        xlims2 = [min([xlim(),xlims2]),max([xlim(),xlims2])];
    end
    
        % get the plots all on the same scaling
    for i=1:nsigs
        subplot(2,nsigs,i);
        xlim(xlims1);
        subplot(2,nsigs,i+nsigs);
        xlims2(1)=2;
        xlim(xlims2);
    end
        
    if (~isempty(figname))
        h.Position=[1200,300,1000,1000];
        fn = fullfile(basedir,sprintf("%s_%04d.tif", figname, fignum));
        saveas(h,fn,'tif');
    end
end


function [RP, DP, stnName, rp] = setup_run(ifn, ss, rp, ncnames, nsets, msites)

    lat = rp.sites.lat(ss);
    lon = rp.sites.lon(ss);
    stnName = rp.sites.stnName{ss};
    stnID = rp.sites.stnID{ss};
    f_info = ARRM_V2_parse_netcdf_filenames(ncnames(1));

    runLbl = sprintf('%s_%.4f_%.4f', underscore(stnName,13),lat,lon);
    
    if (any(strcmp(rp.models,"Stations")))
        runtype = "station";
    else
        runtype = "model";
    end

    fprintf('\n\n--------------Run %4d of %4d, set %3d of %3d  %s %s (%s %s %s)\n\n', ss, msites, ifn, nsets, stnID, stnName, f_info.varname, f_info.model, f_info.ensemble);

    RP  =ARRM_V2_RunParams("temperature", rp.Unmatched, 'runType',{'ARRM_V2','temp'}, ...
                           'pdf_yrstep',30, 'pdf_yrs',31, ...
                           'do_pdfs',true, 'do_base_probs', true,'do_rolling_probs', true, 'do_problines',true, ...
                           'probs', rp.probs, 'do_calc_binning', true ...
                           ); 
    if (runtype == "station")
        DP =ARRM_V2_DataParams(RP.Unmatched, 'fnames',ncnames, 'lats',lat,'lons',lon, 'stations', stnID, ...
                                'base_yrs',rp.base_yrs,'rolling_yrs', [], 'trend_yrs',  rp.trend_yrs,...
                                'scenario',rp.scenario,'ensemble',rp.ensemble, ...
                                'runLbl',runLbl, ... 
                                'cdf_append_pts', rp.cdf_append_pts, ...
                                'figbase',rp.figbase, ...
                                'runType',{'ARRM_V2','temp', runtype});
    else
        DP =ARRM_V2_DataParams(RP.Unmatched, 'fnames',ncnames, 'lats',lat,'lons',lon, ...
                                'base_yrs',rp.base_yrs,'rolling_yrs', rp.rolling_yrs, 'trend_yrs',  rp.trend_yrs,...
                                'scenario',rp.scenario,'ensemble',rp.ensemble, ...
                                'runLbl',runLbl, ...                                           
                                'figbase',rp.figbase, ...
                                'cdf_append_pts', rp.cdf_append_pts, ...
                                'figflags', rp.figflags, ...
                                'runType',{'ARRM_V2','temp', runtype});
    end
    if (~isempty(fields(DP.Unmatched)))
       fprintf('error:  Unused input parameters: \n');
       fprintf('\t%s\n', fields(DP.Unmatched));
       error('error:  Unused input parameters');
    end
        
end



function rp = initparams(varargin)
%             initparams(obs, mdl, lats,lons,base_yrs, varargin)

                % and parse the remaining arguments.
    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = true;

%   default_probs = normcdf([-2.5,-1,0,1,2.5]);
    default_probs  = [0., .00001, .0001,.001, .00621, .01,.0228,.05,.10,.1587,.50,.8413,.90,.95,.9772, .99, .99379, .999,.9999, .99999, 1.0];
%                                                ^           ^             ^         ^            ^            ^
%                                                |           |           +/- 1 std. dev           |            |
%                                                |   ~ -2 st devs (.0228)               ~ +2 st devs (.9772)   |
%                                              -2.5                                                          +2.5

    addParameter(p,"basedir","anom_distrib");
    addParameter(p,"models","stations");
    addParameter(p,"stn_fname","test_stations_1000.csv");
    addParameter(p,"exclude_sites",[]);

    addParameter(p,"base_yrs",[1950,2005]);
    addParameter(p,"trend_yrs",[1935,2005]);
    addParameter(p,"probs",default_probs);
    addParameter(p,"do_calc_binning",true);

    addParameter(p,"ncname",[]);
    addParameter(p,"varname","tasmax");
    addParameter(p,"scenario","hist");
    addParameter(p,"ensemble","r1i1p1");
    addParameter(p,"fig_base",150);
    addParameter(p,"figbase",[]);
    addParameter(p,"figpos",[]);
    addParameter(p,'figname',"default");
    addParameter(p,'anomfile',"default");
    addParameter(p,'read_anomfile',false);
    addParameter(p,'do_plots',true);
    addParameter(p,'partial_match',true);
    addParameter(p,'nsites',[]);
    addParameter(p,'sitenums',[]);
    addParameter(p,"cdf_append_pts",[-3,3.0; -4.75, 4.75]);
    
    parse(p, varargin{:});

    rp = p.Results;
    rp.Unmatched = p.Unmatched;     % so we can pass on 

    setup_dirs(rp.basedir,["anom","log","tif","cdf"]);   % make sure all needed output folders exist.
        
    [rp.models, unmatched] = all_models('match',rp.models,'name',true,'match_station',true,'partial',rp.partial_match,'standard',false);
    
    if (isempty(rp.models))
        fprintf("anomaly_distributions:  error matching model %s\n", rp.models);
        return;
    end
    
    if (~isempty(unmatched))
        fprintf("error:  unmatched models: \n");
        disp(unmatched);
        error("unmatched models");
    end
    
    rp.trend_yrs = [min(rp.trend_yrs(1),rp.base_yrs(1)),max(rp.trend_yrs(2),rp.base_yrs(2))];
    if (~isempty(rp.ncname))
        rp = ARRM_V2_parse_netcdf_filenames(rp.ncname,[],rp);
    else
        [rp.ncnames, rp.base_dirs] = ARRM_V2_make_netcdf_filenames(rp.varname,rp.models, rp.ensemble, rp.scenario, rp.trend_yrs);
    end
    
    rp.sites = readtable(rp.stn_fname);
    if (~isempty(rp.exclude_sites))
        ikeepers = find_matches(rp.sites.stnID, [rp.exclude_sites], false);
        nkeepers = find_matches(rp.sites.stnName,[rp.exclude_sites], true);
        excluded = rp.sites((ikeepers | nkeepers),:);
        rp.sites = rp.sites(~(ikeepers | nkeepers),:);
        fprintf("excluding %d sites\n", size(excluded,1));
        if (size(excluded,1) <= 20), disp(excluded); end
    end
    
        % copy site table into basedir.
    stn_bname = basename(rp.stn_fname);
    if (~exist(fullfile(rp.basedir,stn_bname),'file'))
        writetable(rp.sites,fullfile(rp.basedir,stn_bname));
    end
    
    if (isempty(rp.nsites)), rp.nsites = size(rp.sites,1); end
    
    if (strcmp(rp.anomfile,"default"))
        rp.anomfile = fullfile(rp.basedir,sprintf("anom/anom_distrib_%s_%s_%04d_%04d.mat",rp.models(1), rp.varname, size(rp.sites,1), rp.nsites));
    end
    
    if (strcmp(rp.figname,"default"))
        rp.figname = sprintf("tif/anom_distrib_%s_%s_%04d_%04d",rp.models(1), rp.varname, size(rp.sites,1), rp.nsites);
    end
   
end
