function rpath = canonical_path(fname)
% returns canonical path to fname.
%   uses python script canonical_path.py in same folder as this mfile to find the canonical path.
%   canonical_path.py is given below in comments.
%
    p=fileparts(which(mfilename));
    [status, cmdout] = system(sprintf('%s "%s"', fullfile(p,"canonical_path.py"), fname));
    
    if(status==0)
        rpath=strtrim(cmdout);
    else
        error(cmdout);
    end
end
%
%       canonical_path.py:
%
%---------------

% #!/usr/bin/env python
% #
% #   python script to output the canonical path to a given file.
% #
% #   NOTE:  returns a path even if the file doesn't exist...just resolves all valid paths along the way.
% #
% import os
% import sys
% 
% print(os.path.realpath(sys.argv[1]))

%---------------
