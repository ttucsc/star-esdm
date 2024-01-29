function [varName, varID, xtype, dimids, natts] = ncfind_varinfo_ic(fname, varlist)
%
%   returns varinfo for 1st match of a variable name in varlist.
%   If varlist is missing, use list below.
%
%   NOTE:  varID is (zero-based) netcdf varID, not matlab (one-based) index
    
    if (~exist('varlist','var') || isempty(varlist))
        varlist = ["Tmax","TMAX","Tmin","TMIN","Prec","PREC",'tasmax','tasmin','pr', "temp_F","rh_F","rhsmax","rhsmin","relhum"];
    end
    
    [varNames, varIDs, xtypes, dimidss, nattss] = ncget_varsinfo_ic(fname);
    
        % in case not found.
    varName="";
    varID=[];
    xtype=[];
    dimids=[];
    natts=[];
    
        % find first var in the list.
    for i=1:length(varlist)
        vix = find(strcmp(varlist(i), varNames),1);
        if (~isempty(vix))
            varName = string(varNames{vix});
            varID   = varIDs(vix);
            xtype   = xtypes(vix);
            dimids  = dimidss{vix};
            natts   = nattss{vix};
            break;
        end
    end   
end