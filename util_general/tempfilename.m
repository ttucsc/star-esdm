function tmpname = tempfilename(dir,lbl, ext)
% returns a temporary filename based on the process's PID and a random integer value
% output filename will be in directory dir, w/ name using lbl, pid and a 
% 4-digit random number, and with extension ext.

    if (~exist('dir','var') || isempty_s(dir)), dir='/tmp'; end
    if (~exist('lbl','var') || isempty_s(lbl)), lbl='tmpfile'; end
    if (~exist('ext','var') || isempty_s(ext)), ext='tmp'; end
    
    tmpname = fullfile(dir,sprintf('%s_%d_%d.%s', lbl, feature('getpid'), randi(9999),ext));

end

