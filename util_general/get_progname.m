function [prname, stk] = get_progname(report_all, full_names)
% [progname, stk] = get_progname(report_all, skip_latest)
% returns name of initial program, or list of all programs to
% current depth.  Top entry is 1st program called;  last entry is name of
% file where get_progname is called from.
%   

    if (~exist("report_all","var")), report_all = false; end
    if (~exist("full_names","var")), full_names = false; end

    if (isnumeric(report_all))
        ioff = report_all;
        report_all = false;
    else
        ioff = 0;
    end

    if (full_names)
        stk = dbstack(1,"-completenames");       % 1 skips current call, so reports up to where this function is called.
    else
        stk = dbstack(1);        
    end
    stklen = length(stk)-ioff;
    if (report_all)
        iend = 1;
        prname = strings(stklen,1);
    else
        iend = stklen;
        prname = strings(1,1);
    end

    for ii=stklen:-1:iend
        i = stklen-ii+1;
        prname(i) = stk(ii).file;
    end
end

