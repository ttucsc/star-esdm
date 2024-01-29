function [stns,nsites] = QC_get_site_list_from_file(fname,column,hasHeadings)
%   reads list of sites from fname
%       checks to see if fname is an actual file
%       then reads data from file.
%       if its not a file, simply copies fname to sites and returns
%       this presumes that fname is actually a list of one stnID
%
%   Inputs:  
%       fname       filename
%       column      (optional)  column number or column name to read from
%                               if string or char, then assume file has column headings
%       hasHeadings (optional)  true/false.  If true, treats 1st row as column headings.
%                               Ignored if column is char or string.

    if (~exist('column','var') || isempty(column))
        column=1;
    end
    if (~isnumeric(column))
        hasHeadings=true;        
    elseif (~exist('hasHeadings','var') || isempty(hasHeadings))
        hasHeadings=false;
    end

    if (ischar(fname))
        fnm=fname;
    else
        fnm=char(fname);    % fullfile doesn't like strings...wants input as char vectors...grrr.
    end
          % see if station is a filename.  If so, read station list from file
    [dir,~,~] = fileparts(fnm);
    if (isempty(dir)), fnm=fullfile(pwd,fnm); end   % prepend path so matlab doesn't search searchpath for it
    if (exist(fnm,'file')==2)
        if (hasHeadings)
            tbl=readtable(fnm);
        else            
            tbl=readtable(fnm,'ReadVariableNames',false);
%             if (strcmpi(stns{1},'stnID') || strcmpi(stns{1},'stnName'))
            hasHeadings = looks_like_it_has_labels(tbl);
            if (hasHeadings)   % see if tbl actually has column heading
                tbl=readtable(fnm);
            end
        end
        if (isempty(tbl{end,column}))       % remove last entry if last line was a blank.
            tbl=tbl(1:end-1,:);
        end
        if (isstring(column) || ischar(column))
            stns=tbl.(column);
        else
            stns=tbl{:,column};
        end
    else
        stns=fname;
    end
    nsites=length(stns);
end

function tf = looks_like_it_has_labels(tbl)

    lbls = {'stnID','stnName','station','ID'};

    for i=1:size(tbl,2)
        if (any(strcmpi(tbl{1,i},lbls)))
            tf = true;
            return;
        end
    end
    tf = false;
end

