function c = toCharArray(s, maxlen, maxdigs)
% Convert a cell array of char variables, or array of strings, to a fixed-length array of chars.
% Inputs:
%   s       input array
%               cell array of chars, such as {'hello';'worldliness'};
%               array of strings, such as ["goodbye","moon"];
%               array of integer or real numbers
%   maxlen      # of chars for each element
%                   if missing, use longest string in input
%                   maxlen required for numeric data
%   maxdigs     # of significant digits, for real data
%
% Outputs
%   c           2-D array of chars, with all data space-padded
%                   if maxlen is negative (or missing for text input), results are left-justified
%
    if (nargin == 1), maxlen = []; end
    npts=length(s);
    if (iscell(s))
        if (isempty_s(maxlen)), maxlen = -max(cellfun(@length,s)); end
        c=repmat(' ',npts,abs(maxlen));
        for i=1:length(s)
            c(i,:) = sprintf('%*s',maxlen, s{i});
        end
    elseif (isstring(s))
        if (isempty_s(maxlen)), maxlen = -max(strlength(s)); end
        c=repmat(' ',npts,abs(maxlen));
        for i=1:length(s)
            c(i,:) = sprintf('%*s',maxlen, s(i));
        end
    elseif (isnumeric(s))
        if (isinteger(s))
            if isempty_s(maxlen), error('toCharArray: must specify maxlen for integer arrays'); end
            c=repmat(' ',npts,abs(maxlen));
            for i=1:length(s)
                l = length(s(i));
                if (maxlen < 0)
                    c(i,1:l) = s(i);
                else
                    
                    c(i,(maxlen-l+1):maxlen) = s(i);
                end
            end
        elseif isreal(s)
            if (nargin ~= 3), error('toCharArray: must specify maxlen and maxdigs for real data'); end
            c=repmat(' ',npts,abs(maxlen));
            for i=1:length(s)
                ss= sprintf('%*.*g',maxlen, maxdigs, s(i));
                l=lentgh(ss);
                if (maxlen < 0)
                    c(i,1:l) = ss;
                else
                    c(i,(maxlen-l+1):maxlen) = s;
                end
            end
        else
            error('toCharArray:  input cannot be complex');
        end  
    else
        error('toCharArray:  input must be strings, cell array of chars, or numeric');        
    end
end