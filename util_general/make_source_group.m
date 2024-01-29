function srcgrp = make_source_group(obsnames, histnames, mdlnames, saveMetaData)
% grp = make_source_group(nc_orig)

    obsnames = string(obsnames);
    histnames = string(histnames);
    mdlnames = string(mdlnames);
    
    if (contains(obsnames,"|"))
        obsnames = split(obsnames,"|");
    elseif  (contains(obsnames," "))
        obsnames = split(obsnames," ");
    end
    if (contains(histnames,"|"))
        histnames = split(histnames,"|");
    elseif (contains(histnames," "))
        histnames = split(histnames," ");
    end
    if (contains(mdlnames,"|"))
        mdlnames = split(mdlnames,"|");
    elseif (contains(mdlnames," "))
        mdlnames = split(mdlnames," ");
    end
    
    nc_orig = {obsnames; histnames; mdlnames};

             % add metadata from original NC files
    grp_lbls  = ["OBS_src","HIST_src","MODEL_src"];
    grp_mlbls = ["OBS_src_metadata","HIST_src_metadata","MODEL_src_metadata"];
    grp_descr = ["metadata from OBS source file", "metadata from HIST source file", "metadata from MODEL source file"];
    if (saveMetaData)
        srcgrp = Group("Source_metadata","description","header info from input files");
        srcgrp.putatt("brief",uint8(false));
        for i=1:length(nc_orig)
            if (isempty(nc_orig{i}))
                fns="";
            else
                fns=join(string(nc_orig{i}),'|');
            end
            srcgrp.putatt(grp_lbls{i},fns);
        end            
        for i=1:length(nc_orig)
            fns = nc_orig{i};
            if (length(fns)==1)
                src_nc=ncdf(fns{1},"create", false);
    %                 grp = nc.clone(false, false);  % don't copy data or child groups
                grp = make_metadata_group(grp_mlbls(i), grp_descr(i),src_nc);
                grp.Name = grp_lbls(i);
                srcgrp.putgrp(grp);
            else
                igrp=Group(grp_mlbls{i},"description",strcat("metadata from input files for ",grp_lbls{i}));
                igrp.putatt("filenames",join(fns,"|"));
                srcgrp.putgrp(igrp);
                for j=1:length(fns)
                    src_nc=ncdf(fns{j});
    %                     grp=nc.clone(false,false);
                    lbl  = sprintf("%s_%02d",grp_mlbls(i),j);
                    dscr = sprintf("metadata from %s source file %d", grp_lbls(i),j);
                    grp = make_metadata_group(lbl, dscr,src_nc);
                    igrp.putgrp(grp);
                end
            end
        end
    else
        srcgrp = Group("Source_metadata","description","filenames for sources");
        for i=1:length(nc_orig)
            fns=join(string(nc_orig{i}),'|');
            srcgrp.putatt(grp_lbls{i},fns);
        end            
        srcgrp.putatt("brief",uint8(true));
    end

end