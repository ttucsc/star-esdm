function [tbls, nbad] = readghcn(fnames, vname,dir, outdir)
% Reads a raw GHCN data file

    if (~exist('dir','var') || isempty(dir)), dir='.'; end
    if (~exist('outdir','var') || isempty(outdir)), outdir='.'; end
    fnames=string(fnames);
    nsites = length(fnames);
    tbls = cell(nsites,1);
    yrrange=[nan,nan];
    for i=1:nsites
%        [~,stn,ext]=fileparts(fnames(i));
        [~,stn,~]=fileparts(fnames(i));
%        if (isempty(ext)), ext=".dly"; end
        ext=".dly";
        fname=fullfile(dir,sprintf("%s%s",stn,ext));

        tempvars=["TMAX","TMIN","TAVG"];

        fid=fopen(fname);
        if (fid<0), error('error opening file %s',fname); end

        ghcndata=nan(12*300*31,1);
        gdates=nan(size(ghcndata));
        tline=fgetl(fid);
        dnum0 = 0;
        got_start=false;
        if (any(tempvars==vname))
            scaling = 0.1;
        else
            scaling = 1.0;
        end
        nyrs=0;
        while(ischar(tline))
            [stnID, dnums, varname,mdata, nbad] = ghcn_parse(tline);
            if (strcmp(varname,vname) && strcmp(stnID, stn))
                if (~got_start)
                    got_start = true;
                    dnum0 = dnums(1)-1;
                end
                ix = dnums-dnum0;
                ghcndata(ix) = mdata * scaling;

                gdates(ix) = dnums;
                dvec=datevec(dnums(1));
                year=dvec(1);
                month=dvec(2);
                if (month==1)
                    fprintf('%4d ',year);
                    nyrs=nyrs+1;
                    if (mod(year,20)==0), fprintf('\n'); end
                end
                yrrange=[nanmin(yrrange(1),year),nanmax(yrrange(2),year)];
                lastline=tline;
                lastix = ix;
                lastdnums = dnums;
            end
            tline=fgetl(fid);
        end
        fprintf('\n');
        lastix = find(~isnan(ghcndata),1,'last');
        ghcndata((lastix+1):end)=[];
        gdates((lastix+1):end)=[];
        dvecs = datevec(gdates);
        yr = dvecs(:,1);
        mo = dvecs(:,2);
        da = dvecs(:,3);
        tbl = table(yr,mo,da, ghcndata,'VariableNames',{'year','month','day',vname});
        outname=fullfile(outdir,sprintf('%s.csv', stnID));
        writetable(tbl,outname);
        tbls{i} = tbl;
        fprintf('file %d written: %s\n', i, outname);
    end
    fprintf('year range: %d %d\n', yrrange);
    if (nsites==0)
        tbls=tbls{1};
    end
end


function [stnID, dnums, varname,mdata, nbad] = ghcn_parse(tline)

    stnID = tline(1:11);
    yr      = str2double(tline(12:15));
    mo      = str2double(tline(16:17));
    varname = tline(18:21);
    dnums = datenum(yr,mo,0)+(1:31)';
    mdata = nan(31,1);
    mflag = blanks(31);
    qflag = blanks(31);
    sflag = blanks(31);
    nbad=0;
    for i=1:31
        ix1 = 21+(i-1)*8 + 1;
        ix2 = ix1 + 4;
        ixm = ix1 + 5;
        ixq = ix1 + 6;
        ixs = ix1 + 7;
        try
            mdata(i) = str2double(tline(ix1:ix2));
            mflag(i) = tline(ixm);
            qflag(i) = tline(ixq);
            sflag(i) = tline(ixs);
        catch
            nbad=nbad+1;
        end
    end
    mdata(mdata==-9999)=nan;
    ok = mflag == ' ' & qflag == ' ';
    mdata(~ok) = nan;
end
