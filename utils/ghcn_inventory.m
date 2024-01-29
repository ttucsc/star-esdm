function G = ghcn_inventory(varName, basedir)

    if (~exist('basedir','var') || isempty_s(basedir)), basedir = '.'; end
    fname=fullfile(basedir,sprintf('ghcnd-inventory_%s.txt',varName));
    fid = fopen(fname);
    if (fid == -1)
        throw(MException('ARRMV2:NOSUCHFILE',sprintf('error:  no such file: %s\n',fname)));
    end
        
    try
        G = cell(1000000,5);

        fclose(fid);
        fid = fopen(fname);
        nlines=0;
        while (fid)
            l = fgetl(fid);
            if (isempty_s(l))
%                fprintf('end of file: %6d\n', nlines);
                break; 
            end
            if (length(l) < 45)
                if (l(1)==-1), break;
                else, continue; 
                end
            end
            g = ghcn_parse(l);
            if (g{5}-g{4} >= 15 && g{4} < 2000 && g{5} >= 1915)
                nlines=nlines+1;
%                if (mod(nlines,5000)==0), fprintf('%6d\n', nlines); end
                G(nlines,:) = g;
            end
        end
        G(nlines+1:end,:) = [];
%        fprintf('%6d\n', nlines);
%         G = cell2table(G(1:nlines,:),'VariableNames',{'stn_id','lat','lon','start_yr','end_yr'});
        G = cell2struct(G(1:nlines,:),{'stn_id','lat','lon','start_yr','end_yr'},1);
        fclose(fid);
    catch
        if (fid ~= -1), fclose(fid); end
    end
end
function G = ghcn_parse(l)

    G = cell(1,5);
    G{1} = string(l(1:11));

    G{2} = str2double(l(13:20));
    G{3} = str2double(l(22:30));
    G{4} = str2double(l(37:40));
    G{5} = str2double(l(42:45));

end
        
        
        
