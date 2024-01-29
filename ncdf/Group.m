classdef Group < ncObj
    %ncdf Group class.  Parent class of ncdf objects.
    %   Can be the netcdf4's global (base) group, or any group under the global group. 
%
%   12/2/2021:  changed fix_unlimited_dimensions(...) to clear_unlimited_dimensions(...)
%               and added a "set_unlimited_dimensions(...) function.
    
    properties
        Name;
        Filename;
        Dimensions;
        Variables;
        Attributes;
        Groups;
        Format;
    end
    
    methods
        function obj = Group(nm, varargin)
            if (nargin == 0), return; end
            if (isstruct(nm) )  % assume a netcdf Variable struct
                mygrp = nm;
                flds = fieldnames(mygrp);
                for i=1:length(flds)
                    fld = flds{i};
                    if (isprop(obj,fld))         % copy only fields that are properties
                        if (obj.isnoncatprop(fld))
                            obj.(fld) = mygrp.(fld);
                        elseif (strcmp(fld,'Dimensions'))
                            dims=mygrp.Dimensions;
                            for j=1:length(dims), obj.put('Dimensions',Dimension(mygrp.Dimensions(j))); end
                        elseif (strcmp(fld,'Variables'))
                            vars=mygrp.Variables;
                            for j=1:length(vars), obj.put('Variables',Variable(mygrp.Variables(j))); end
                        elseif (strcmp(fld,'Attributes'))
                            atts=mygrp.Attributes;
                            for j=1:length(atts), obj.put('Attributes',Attribute(mygrp.Attributes(j))); end
                        elseif (strcmp(fld,'Groups'))
                            grps=mygrp.Groups;
                            for j=1:length(grps)
                                obj.put('Groups',Group(mygrp.Groups(j))); 
                            end
                        else
                            error("error:  unknown field %s in %s\n", fld, obj.Name);
                        end
                    end
                end
            elseif (isa(nm,'Group'))       % copy the properties one by one, in case name is a derived class of Group
                props=properties(obj);
                for i=1:length(props)
                    prop=props{i};
                    obj.(prop) = nm.(prop);
                end
            else                
                if (nargin >= 1), obj.Name = char(nm); end
            end
            if (~isempty(varargin))         % now insert any additional info specified in varargin.
                for i=1:2:length(varargin)
                    arg = varargin{i};
                    val = varargin{i+1};
                    if (isprop(obj,arg))
                        if (iscell(val))
                            obj.put(arg, val{:});
                        else
                            obj.put(arg, val);
                        end
                    else
                        obj.put('Attributes',Attribute(arg,val));
                    end
                end
            end
        end
        
        function vars = varlist(obj, i)        % returns a string array containing qualified names of all variables in group, or ith variable name if i is specified.
            if (isempty(obj.Variables))
                vars=strings(0,0);
            else
                vars = string({obj.Variables.Name});                
                if (exist("i","var") && ~isempty(i))
                    vars = vars(i);
                end
            end
        end
        
        function vars = varlist_deep(obj, deep)        % returns a cell array containing qualified names of all variables in group, recurses into sub-groups.
            if (isempty(obj.Variables))
                vars=strings(0,0);
            else
                vars = string({obj.Variables.Name});
            end
            if (exist('deep','var') && ~deep), return; end
            for i=1:length(obj.Groups)
                vl = obj.Groups(i).varlist();
                len = length(vl);
                if (len > 0)
                    gname=obj.Groups(i).Name;
                    vn=strings(1,len);
                    for j=1:len
                        vn(j)=sprintf("%s/%s", gname, vl(j));
                    end
                    vars = cat(2,vars,vn);
                end
            end
        end
        
        function grps = grplist(obj,deep)        % returns a cell array containing qualified names of all Groups in group
                                                 % if deep is present and true, includes qualified names of subgroups.
            if (isempty(obj.Groups)), grps=strings(0,0); return; end
            grps = string({obj.Groups.Name});
            if (~exist('deep','var') || ~deep), return; end
            for i=1:length(obj.Groups)
                gl = obj.Groups(i).grplist();
                len = length(gl);
                if (len > 0)
                    gname=obj.Groups(i).Name;
                    gn=strings(1,len);
                    for j=1:len
                        gn(j)=sprintf("%s/%s", gname, gl(j));
                    end
                    grps = cat(2,grps,gn);
                end
            end
        end
        
        function atts = attlist(obj,deep)        % returns a cell array containing qualified names of all attributes in group
                                                 % if deep is present and true, includes qualified names of attributes
                                                 % in subgroups.
            if (isa(obj,'Attribute')), atts=obj.Name; return; end
            atts = string({obj.Attributes.Name});
            if (~exist("deep","var") ||  ~deep), return; end
            for i=1:length(obj.Groups)
                al = obj.Groups(i).attlist(deep);
                len = length(al);
                if (len > 0)
                    gname=obj.Groups(i).Name;
                    an=strings(1,len);
                    for j=1:len
                        an(j)=sprintf("%s/%s", gname, al(j));
                    end
                    atts = cat(2,atts,an);
                end
            end
        end
        
        function att_struct = attstruct(obj, deep)
            if (~exist('deep','var')), deep = false; end
            atts = obj.attlist();
            for ax=1:length(atts)
                a = obj.Attributes(ax);
                att_struct.(a.Name) = a.Value;
            end
            if (deep)
                for gx=1:length(obj.Groups)
                    g = obj.Groups(gx);
                    gn = g.Name;
                    att_struct.(gn) = g.attlist(deep);
                end
            end
        end
                    
        function dims = dimlist(obj, deep)        % returns a cell array containing qualified names of all dimensions in group
            if (isempty(obj.Dimensions))
                dims=strings(0,0);
            else
                dims = string({obj.Dimensions.Name});
            end
            if (exist('deep','var') && ~deep), return; end
            for i=1:length(obj.Groups)
                dl = obj.Groups(i).dimlist();
                len = length(dl);
                if (len > 0)
                    gname=obj.Groups(i).Name;
                    dn=strings(1,len);
                    for j=1:len
                        dn(j)=sprintf("%s/%s", gname, dl{j});
                    end
                    d2=dims;        % (this gets rid of a matlab warning about preallocating)
                    dims = cat(2,d2,dn);
                end
            end
        end
        
         function vars = nondimvarlist(obj,deep)        % returns a cell array containing qualified names of all variables which are NOT dimensions
            if (~exist('deep','var') || isempty_s(deep)), deep = true; end
            vars = obj.varlist(deep);
            dims = obj.dimlist(deep);
            if (isempty(dims)), return; end
            for i=1:length(dims)
                dimix = find(strcmp(vars,dims(i)),1);
                if (~isempty(dimix))
                    vars(dimix) = [];
                end
            end
         end
         
         function var = getvar(obj, varargin)
             var=obj.get("Variables",varargin{:});
         end
%          function dim = getdim(obj, varargin)
%              dim = obj.get("Dimensions",varargin{:});
%          end
%          function att = getatt(obj, varargin)
%              att = obj.get("Attributes",varargin{:});
%          end
         function grp = getgrp(obj, varargin)
             grp = obj.get("Groups",varargin{:});
         end
             
        
         function myclone = clone(obj, do_copy_data, do_copy_groups)
             if (~exist('do_copy_data',  'var') || isempty(do_copy_data)),   do_copy_data   = true; end
             if (~exist('do_copy_groups','var') || isempty(do_copy_groups)), do_copy_groups = true; end
             
            if (isa(obj,'ncdf'))
                myclone = ncdf('', 'Name',obj.Name, 'Filename',obj.Filename, 'Format',obj.Format);
            else
                myclone = Group('', 'Name',obj.Name);
            end

                        % copy over all the attributes
            for i=1:length(obj.Attributes)
                att = obj.Attributes(i).clone();
                myclone.putatt(att);
            end

                        % copy over all the dimensions
            for i=1:length(obj.Dimensions)
                dim = obj.Dimensions(i).clone();
                myclone.putdim(dim);
            end

                        % copy over all the Variables
            for i=1:length(obj.Variables)
                var = obj.Variables(i).clone(myclone, do_copy_data);
                myclone.putvar(var);
            end

                        % copy over all the Groups
            if (do_copy_groups)
                for i=1:length(obj.Groups)
                    grp = obj.Groups(i).clone(do_copy_data, do_copy_groups);
                    myclone.putgrp(grp);
                end
            end
         end             
         
         function icwriteschema(obj, ncid, nofillmode)
         % writes out the group schema (using low-level netcdf functions)
         
                    % write the global attributes
             NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');
             for i=1:length(obj.Attributes)
                 att=obj.Attributes(i);
                 try
 %                  fprintf("--- name %s  value %s\n", att.Name, att.Value);
                    netcdf.putAtt(ncid, NC_GLOBAL,att.Name, att.Value);
                 catch me
                     fprintf('icwriteschema putAtt oops!\n');
                       report_me_error(me);
                 end
             end
                    % write out the dimensions
%              dimids = zeros(length(obj.Dimensions);
             for i=1:length(obj.Dimensions)
                 dim = obj.Dimensions(i);
                 if (dim.Unlimited)
%                     dimids(i)=netcdf.defDim(ncid, dim.Name, netcdf.getConstant('NC_UNLIMITED'));
                    netcdf.defDim(ncid, dim.Name, netcdf.getConstant('NC_UNLIMITED'));
                 else
%                     dimids(i)=netcdf.defDim(ncid, dim.Name, dim.Length);
                    netcdf.defDim(ncid, dim.Name, dim.Length);
                 end                 
             end
             
                    % write out the variables
             netcdf.inqDimIDs(ncid);
             for i=1:length(obj.Variables)
                 v=obj.Variables(i);
%                  fprintf("writing variable %s\n", v.Name);
                 dimids = v.dimids(ncid);
                 varid = netcdf.defVar(ncid, v.Name, v.xtype(obj.Format), dimids);
                 if (~isempty_s(v.FillValue))
%                      fprintf("defining FillValue %f\n", v.FillValue);
                    myfillvalue = cast(v.FillValue,v.Datatype);               % make sure VarFill's data type matches the actual datatype.  For example, there's a problem with 1e20 ~= single(1e20)
                    try
                        netcdf.defVarFill(ncid,varid, nofillmode, myfillvalue); 
                    catch
                        fprintf("icwriteschema fillmode oops\n");
                    end
                 end

                 if (~isempty_s(v.ChunkSize))           % empty ChunkSize uses netcdf defaults, which may not be very good...
                     if (ischar_s(v.ChunkSize))         % if char, can be set to CONTIGUOUS or DEFAULT
                         if (strcmpi(v.ChunkSize,"CONTIGUOUS"))
                            try
                                v.ChunkSize = [];
                                netcdf.defVarChunking(ncid, varid, 'CONTIGUOUS');
                            catch me
                                fprintf("error encountered setting chunking for variable %s : %s", v.Name, v.ChunkSize);
                                rethrow(me);
                            end
                         elseif (~strcmpi(v.ChunkSize,"DEFAULT"))
                             v.ChunkSize = [];
                             error("error: icwriteschema: bad ChunkSize:  variable %s : %s", v.Name, v.ChunkSize);
                         end
                         
%                      fprintf("defining ChunkSize %f\n",v.ChunkSize);
                     else
                         try
                            netcdf.defVarChunking(ncid, varid, 'CHUNKED',v.ChunkSize); 
                         catch me
                             fprintf(2,"error: icwriteschema: bad chunk size:  variable %s : %s", v.Name, vec2string(v.ChunkSize, "brackets",'[]')); 
                             rethrow(me)
                         end
                     end
                     if (~isempty_s(v.DeflateLevel) && v.DeflateLevel > 0)
                         netcdf.defVarDeflate(ncid,varid, v.Shuffle,true, v.DeflateLevel); 
                     end
                 end
                 for j=1:length(v.Attributes)
                    att=v.Attributes(j);
                    if (strcmpi(att.Name,'_FillValue'))
                        myattfillvalue = cast(att.Value,v.Datatype);    % make sure data type matches variable's datatype.  For example, there's a problem with 1e20 ~= single(1e20)
                        if (isempty(v.FillValue))
                            netcdf.defVarFill(ncid,varid, nofillmode, myattfillvalue); 
%                       elseif (v.FillValue ~= att.Value)
                        elseif (myfillvalue ~= myattfillvalue)
                            fprintf("warning;  Mismatch between v.FillValue and v.attribute's FillValue\n");
                            disp(v);
                            disp(att);
                        else
                            netcdf.defVarFill(ncid,varid, nofillmode, myattfillvalue);                             
                        end
                    else
                        netcdf.putAtt(ncid,varid,att.Name, att.Value);
                    end
                 end
%                  if (~isempty_s(v.FillValue))
%                      netcdf.defVarFill(ncid,varid, false, v.FillValue); 
%                  end
             end
             
             for i=1:length(obj.Groups)
                 g=obj.Groups(i);
                 gpid = netcdf.defGrp(ncid, g.Name);
                 g.icwriteschema(gpid, nofillmode);
%                  try
%                     netcdf.close(gpid);
%                  catch
%                      fprintf('oops!\n');
%                  end
             end
         end
         
         function initvars(obj, ncid, pncid)         % group ncid and parent (base) ncid.
             
             for i=1:length(obj.Variables)
                 v=obj.Variables(i);
                 FillVal = v.FillValue;
                 nd=length(v.Dimensions);
                 if (~isempty_s(FillVal))
                    varid = netcdf.inqVarID(ncid, v.Name);
                    start=zeros(nd,1);                 
                    for j=1:nd
                        start(j)=v.Dimensions(j).Length-1;
                    end
                    if (any(start<0)), continue; end % if nothing to write for this variable, go on to next.
%                    fprintf('nc initializing var %s: ', v.Name);
%                    t1=tic();
                    oldfill = netcdf.setFill(ncid, 'NC_NOFILL');
%                    fprintf('old fill:  %g   ', oldfill);
                    netcdf.putVar(ncid, varid, start, v.FillValue);
                    netcdf.setFill(ncid, oldfill); 
%                    newfill = netcdf.setFill(ncid, oldfill); 
%                    fprintf('new fill:  %g   ', newfill);
%                    elapsed = toc(t1);
%                    fprintf('%8.3f\n', elapsed);
                 end
             end
             for i=1:length(obj.Groups)
                 g = obj.Groups(i);
                 gpid = netcdf.inqNcid(pncid, g.Name);
                 g.initvars(gpid, pncid);
             end
         end
         
         function clear_unlimited_dimensions(obj, dims)
                % clears Unlimited flag on one or more dimensions.
                % Also clears Unlimited flag for the variables which depend on the specified dimensions.

            
            mydimlist = obj.dimlist();
            if (~exist('dims','var') || isempty(dims))
                dims = [];
                mydims = mydimlist;
            else
                mydims = dims;
            end
             
            for i=1:length(mydims)
                if (isnumeric(mydims))
                    if (mydims(i) < 0 || mydims(i) > length(mydimlist)) 
                        continue; 
                    end
                    d = obj.Dimensions(mydims(i));
                else
                    if (~any(strcmp(mydims(i), mydimlist)))
                        continue; 
                    end
                    d = obj.getdim(mydims(i));
                end
                 
                d.Unlimited = false;
                  
                for j=1:length(obj.Variables)
                    jx = find(strcmp(d.Name, obj.Variables(j).dimlist()));
                    if (~isempty(jx))
                        obj.Variables(j).Dimensions(jx).Unlimited = false;
                    end
                end
            end
            
            if (~isempty(obj.Groups))
                for k=1:length(obj.Groups)
                    obj.Groups(k).clear_unlimited_dimensions(dims);
                end
            end
         end
         
         function set_unlimited_dimensions(obj, dims, only)
                % sets Unlimited flag on one or more dimensions.
                % Also sets Unlimited flag for the variables which depend on the specified dimensions.
                % if only is true, then clears all other dimensions which are unlimited.

            
            mydimlist = obj.dimlist();
            if (~exist('dims','var') || isempty(dims))
                dims = [];
                mydims = mydimlist;
            else
                mydims = dims;
            end
            
            if (exist('only','var') && only)
                obj.clear_unlimited_dimensions();
            end
             
            for i=1:length(mydims)
                if (isnumeric(mydims))
                    if (mydims(i) < 0 || mydims(i) > length(mydimlist)) 
                        continue; 
                    end
                    d = obj.Dimensions(mydims(i));
                else
                    if (~any(strcmp(mydims(i), mydimlist)))
                        continue; 
                    end
                    d = obj.getdim(mydims(i));
                end
                 
                d.Unlimited = true;
                  
                for j=1:length(obj.Variables)
                    jx = find(strcmp(d.Name, obj.Variables(j).dimlist()));
                    if (~isempty(jx))
                        obj.Variables(j).Dimensions(jx).Unlimited = true;
                    end
                end
            end
            
            if (~isempty(obj.Groups))
                for k=1:length(obj.Groups)
                    obj.Groups(k).clear_unlimited_dimensions(dims);
                end
            end
         end
         
    end
    
end

