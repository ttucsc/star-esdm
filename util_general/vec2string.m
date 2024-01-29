function s = vec2string(v, varargin)
% function s = vec2string(v, "format", fmt, "sep", sep, "row_sep", row_sep, "brackets", brackets)
% 
%   creates a 1-line string representation of the values in vector or array v.
%   Inputs:
%       v           vector of values to output. can be numeric;  output will use matlab's default numeric-to-string conversions.
%   Optional key/value pairs:
%       "format", fmt        optional format string to use to output each element.  Example:  "format", "%6.2f "
%       "sep", sep              optional.  char or string to use as separator between elements.  Default:  ', '
%       "row_sep", row_sep      char or string to use as separator between rows       defalt:   '; '
%       "brackets:, brackets    optional.  2-element vector used to bracket the output.  Default:  none
%                               examples:  "brackets",'[]'  or "brackets", ["-->","<--"]
%

    [fmt, sep, row_sep, brackets] = init_params(varargin{:});
    
    if (isempty(v))
        s = brackets;
    else
%       if (~exist("fmt","var") || isempty(fmt)), fmt="%f "; end
        [nr,nc] = size(v);
        s="";
        for i=1:nr
            if (i>1), s= s + row_sep; end
            if (~isempty(fmt))
                for j=1:nc
                    s = s + sprintf(fmt, v(i,j)); 
                    if (j ~= nc), s = s + sep; end
                end
            else
                for j=1:nc
                    s = s + string(v(i,j)); 
                    if (j ~= nc), s = s + sep; end
                end
            end
        end
        if (~isempty(brackets))
            s = sprintf("%s%s%s",brackets(1), s, brackets(2));
        end
    end
end

function [fmt, sep, row_sep, brackets] = init_params(varargin)

                    % parse input
    p = inputParser;

                    % these are the params we want to handle in ARRM_V2_wrapper
    addParameter(p,"format", [],    @(s) isempty(s) || ischar(s) || isstring(s));
    addParameter(p,"sep", ", ",     @(s) isempty(s) || ischar(s) || isstring(s));
    addParameter(p,"row_sep", "; ", @(s) isempty(s) || ischar(s) || isstring(s));
    addParameter(p,"brackets", [],  @(s) isempty(s) || ((ischar(s) || isstring(s)) && (length(s)==2 || strlength(s)==2)));
        
    parse(p, varargin{:});

    fmt      = string(p.Results.format);
    sep      = string(p.Results.sep);
    row_sep  = string(p.Results.row_sep);
    brackets = char(p.Results.brackets);
    
    sep     = strrep(sep,     "\n", newline);
    row_sep = strrep(row_sep, "\n", newline);
end
