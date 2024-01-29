function oops_msg = oops(msg, show_stack)
    % for debugging.  Prints some stack info and returns.
    %   oops() to print current prog name and line number.
    %   oops("some text") to display text along with line info
    %   oops(true)          to display stack
    %   oops(false)         just return (useful to just set break on oops line)
    %   oops("some text", true);    display text and stack
    %   Either msg or show_stack can also be an MException.  This will report the MException error.  Example:
    %       catch me
    %           oops("something's wrong!", me);
    %       end
    
    ST = dbstack;
    ST=ST(2:end);
    oops_msg = sprintf("%s:%d:  oops.\n", ST(1).file, ST(1).line);
    if (~exist("msg", "var"))
        fprintf(2, oops_msg);
        return;
    end
    if (isa(msg, "MException"))
        report_me_error(msg);
        return;
    elseif (islogical(msg))
        if (~msg), return; end
        show_stack = true;
    else        
        oops_msg = sprintf("%s  %s", string(msg), oops_msg);
    end
    if (~exist("show_stack","var") || isempty(show_stack)), show_stack = false; end
    
    if (isa(show_stack, "MException"))
        fprintf(2, oops_msg);
        report_me_error(show_stack);
    else
        fprintf(2, oops_msg);
        if (show_stack)
            disp(ST);
        end
    end
end