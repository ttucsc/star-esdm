function f = fullfile_s(varargin)
%   like fullfile, but works with strings.  NOTE:  commented out:  replaces any spaces with underscores.
%
%   This is now handled properly in R2018a and later;  needed only for earlier versions of matlab.

    emptys = isempty_s(varargin{:});       % remove any empty strings.
    varargin(emptys) = [];

    args = cell(size(varargin));
    for i=1:length(varargin)
        if (isstring(varargin{i}))
            args{i} = char(varargin{i});
        else
            args{i} = varargin{i};
        end
    end
%    args=underscore(args);
    f = string(fullfile(args{:}));
end
