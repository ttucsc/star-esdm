function [totSize, info]  = getSize(obj, lbl, do_disp)
%   Returns the total size of all the members of obj (public & private)


    if (nargin==1), lbl=""; end
    if (nargin<3), do_disp=false; end
    
    info="";
    if (isstruct(obj))
        props = fieldnames(obj);
        pubprops = props;
    elseif (isobject(obj) && ~isstring(obj))
        warning off MATLAB:structOnObject
        pubprops = properties(obj);
        obj=struct(obj);
        props = fieldnames(obj);
        warning on MATLAB:structOnObject
    else
        s=whos('obj');
        totSize = s.bytes;
        info=sprintf("%s: %8d",lbl,s.bytes);
        return;
    end
    
    totSize = 0;
    len=max(strlength(props));
    for ii=1:length(props)
        prop = obj.(props{ii}); 
        if (~ismember(pubprops,props{ii}))
            ptxt = "  (private)";
        else
            ptxt = "";
        end
        if (isobject(prop) || isstruct(prop))
            [s.bytes,txt]=getSize(prop,props{ii});
            txt=sprintf("\n%s",txt);
        else
            s = whos('prop');
            txt=props{ii};
        end
        totSize = totSize + s.bytes;

        if (strlength(info)==0)
            fmt = sprintf("%%s%%s %%%ds: %%8d %%s", len);
        else
            fmt = sprintf("%%s\\n%%s %%%ds: %%8d %%s", len);
        end
        info = sprintf(fmt, info, lbl, txt, s.bytes, ptxt);
    end  
    
    if (do_disp)
        disp(info);
    end
end
