classdef ncObj < handle
    %ncOBJ base class for ncdf objects.
    %   inherited by ncdf
    
    properties(Constant = true, Access = private)
        catnames   = {'Variables','Attributes','Dimensions','Groups','ncdf'};
        classnames = {'Variable', 'Attribute', 'Dimension', 'Group', 'ncdf'};
        % list of all valid properties for each type
        allfields = {'Name','Value','Length','Unlimited','Size','DataType','ChunkSize','FillValue','DeflateLevel','Shuffle','Filename','Format','vdata'};
        fields = {  {'Name','Dimensions','Size','Datatype','Attributes','ChunkSize','FillValue','DeflateLevel','Shuffle','vdata'}, ...
                    {'Name','Value'}, ...
                    {'Name','Length','Unlimited'}...
                    {'Name','Dimensions','Variables','Attributes','Groups','Filename','Format'}, ...
                    {'Name','Dimensions','Variables','Attributes','Groups','Filename','Format'}, ...
                };  
        % list of valid properties which are not categories (arrays of values)
        noncats = { {'Name','Size','Datatype','ChunkSize','FillValue','DeflateLevel','Shuffle','vdata'}, ...
                    {'Name','Value'}, ...
                    {'Name','Length','Unlimited'} ...
                    {'Name','Filename','Format'}, ...
                    {'Name','Filename','Format'}, ...
                  };      
    end
        
    methods
        function obj = ncObj(name, cat)
            if (nargin < 2), return; end
            if     (strcmp(cat,'Variables')),  obj = Variable(name);
            elseif (strcmp(cat,'Attributes')), obj = Attribute(name,[]);
            elseif (strcmp(cat,'Dimensions')), obj = Dimension(name);
            elseif (strcmp(cat,'Groups')),     obj = Group(name);
            else
                error('ncobj():  invalid category: %s', cat);
            end
        end

        function tf = isempty_s(obj)
            if (length(obj) > 1), tf = false; return; end
            flds = fieldnames(obj);
            for i=1:length(flds)
                if (~isempty_s(obj.(flds{i})) && ~(strcmp(flds{i},'Shuffle') && obj.(flds{i})==0))
                    tf = false; 
                    return; 
                end
            end
            tf = true;
        end
        function ix = grpix(obj, name)
%             if (~strcmp('Group',class(obj)))
            if (~isa(obj,'Group'))
                ix = [];
            elseif (isempty(obj.Groups)), ix = [];
            else
                ix = find(strcmp(string({obj.Groups.Name}),name));
            end
        end
        function ix = varix(obj, name)
            if (~isa(obj,'Group')) 
                ix = [];
            elseif (isempty(obj.Variables)), ix = [];
            else
                ix = find(strcmp(string({obj.Variables.Name}),name));
            end
        end
        function ix = dimix(obj, name)
            if (~find(strcmp({'Group','Variable'},class(obj))))
                ix = [];
            elseif (isempty(obj.Dimensions)), ix = [];
            else
                ix = find(strcmp(string({obj.Dimensions.Name}),name));
            end
        end
        function ix = attix(obj, name)
            if (~find(strcmp({'Group','Variable'},class(obj))))
                ix = [];
            elseif (isempty(obj.Attributes)), ix = [];
            else
                ix = find(strcmp(string({obj.Attributes.Name}),name));
            end
        end
        
        function [ix, cat] = findix(obj, name, cat, do_create)
            % returns index and category of name in object.  descends tree to find name if necessary.
            % adds name to object if it doesn't exist yet and do_create is true.
            %
            if (~exist('cat','var')), cat=[]; end    
            if (~exist('do_create','var') || isempty_s(do_create)), do_create = false; end
            if (~isempty_s(cat))
                if (strcmp(cat,'Variables'))
                    ix = obj.varix(name);
                    if (isempty(ix) && do_create), ix = obj.putvar(Variable(name)); end
                elseif (strcmp(cat,'Attributes'))
                    ix = obj.attix(name);
                    if (isempty(ix) && do_create), ix = obj.putatt(Attribute(name,'')); end
                elseif (strcmp(cat,'Dimensions'))
                    ix = obj.dimix(name);
                    if (isempty(ix) && do_create), ix = obj.putdim(Dimension(name)); end
                elseif (strcmp(cat,'Groups'))
                    ix = obj.grpix(name);
                    if (isempty(ix) && do_create), ix = obj.putgrp(Group(name)); end
                else
                    error('ncObj:findix(): unknown category %s', cat);
                end
            else
                for i=1:length(ncObj.catnames)
                    cat = ncObj.catnames{i};
                    if (isprop(obj,cat))
                        ix = obj.findix(name,cat);
                        if (~isempty(ix)), return; end
                    end
                    cat=[];
                end
                if (do_create), error('ncObj:findix(%s), cannot create:  cat is empty', name); end
            end
        end    
        
        function dim = getdim(obj, varargin)
            dim = obj.get("Dimensions",varargin{:});
        end
        function att = getatt(obj, varargin)
            att = obj.get("Attributes",varargin{:});
        end
        function data = getvardata(obj, arg)        % returns vdata for arg (just the data).  arg can be qualified
            data = obj.get(sprintf('%s/vdata',arg));
        end

        function data = getattvalue(obj, arg)        % returns attribute specified by arg.  arg can be qualified
            data = obj.get(sprintf('%s/Value',arg));
        end

        function outobj = get(obj, varargin)         % returns the object specified by arg.  (Variable, Group, Attribute, Dimension)  
            % returns item within obj named arg.
            % arg may be qualified name, in which case it descends the tree to find it.
            % if a variable and a dimension of same name exists, returns the variable, unless arg qualifies it to
            % request the dimension, as in:
            %
            %       myobj.get('/Dimensions/latitude') or myobj.get("Variables/time","units/Value")
            if (length(varargin)==1)
                arg = varargin{1};
            else
                arg = join(string(varargin),"/");
            end
                
            spath=ncObj.split(arg);
            name = spath{1};
            spath = spath(2:end);
            
            if (obj.isnoncatprop(name))
                if (~isempty_s(spath)), error('ncObj:put:  can''t get regular property, ''%s'', except when specified at end-of-path', name); end
                outobj = obj.(name);
                return;
            end
            
            [~, cat] = ncObj.findcat(name);       % see if 1st entry in spath is a category name
            if (~isempty_s(cat))                  % Yes, it is a category.
                if (length(spath)>=1)            % Still more of path to parse.  descend to next level
                    name = spath{1};
                    ix = obj.findix(name, cat);   % find index of name in category (insert if not present)
                    if (isempty(ix))
%                         if (~any(strcmp(name, ["_FillValue","GroupParams"])))
%                             fprintf("name=%s  cat=%s  arg=%s\n", name, cat, arg);
%                         end
                        error('error: ncobj.get(%s): cannot find matching item in %s', name, cat); 
                    end
                    if (length(spath) > 1)
                        outobj = obj.(cat)(ix).get(spath(2:end));
                    else % path is only 1 long.
                        outobj = obj.(cat)(ix);
                    end
                else
                    outobj = obj.(cat);
                end
            else
                if (isprop(obj,name))
                    outobj = obj.(name);
                else
                    [ix,cat] = findix(obj,name);
                    if (isempty_s(cat))
                        error('error:ncObj.get():  cannot locate %s in %s', name, obj.Name);
                    end
                    if (isempty_s(spath))
                        outobj = obj.(cat)(ix);
                    else
                        outobj = obj.(cat)(ix).get(spath);
                    end
                end
            end
        end
        
        
        
        function put(obj, arg, varargin)        % inserts arg into obj
            spath=ncObj.split(arg);
            name = spath{1};
            spath = spath(2:end);
            
            if (obj.isnoncatprop(name))
                if (~isempty_s(spath)), error('ncObj:put:  can''t add regular property except to end-of-path object'); end
                if (length(varargin)>1), error('ncObj:put:  too many parameters: %s', name); end
                obj.(name) = varargin{1};
                if (strcmp(name,'vdata') && isa(obj,'Variable'))
                    obj.Size = ncObj.mysize(obj.vdata); 
                elseif (strcmp(name,'Filename'))
                    for g=1:length(obj.Groups)  % propogate new filename to all the subgroups.
                        obj.Groups(g).put('Filename',varargin{1});
                    end
                elseif  (strcmp(name,'FillValue'))
                    obj.check_FillValue();          % make sure FillValue and any _FillValue attribute agree.
                end
                return;
            end
            
            [~, cat] = ncObj.findcat(name);       % see if 1st entry in spath is a category name
            if (~isempty_s(cat))                  % Yes, it is a category.
                if (length(spath)>1)            % Still more of path to parse.  descend to next level
                    name = spath{1};
                    ix = obj.findix(name, cat);   % find index of name in category (insert if not present)
                    if (isempty(ix))
                        obj.put(sprintf('%s/%s',cat,name));
                        ix = obj.findix(name, cat);   % find index of name in category (insert if not present)
                    end
                    obj.(cat)(ix).put(spath(2:end), varargin{:});
                elseif (length(spath) == 1)     % path is only 1 long.
                    name = spath{1};
                    if     (strcmp('Variables',cat)),   obj.putvar( name, varargin{:});
                    elseif (strcmp('Attributes',cat)),  obj.putatts(name, varargin{:});
                    elseif (strcmp('Dimensions',cat)),  obj.putdim( name, varargin{:});
                    elseif (strcmp('Groups',cat)),      obj.putgrp( name, varargin{:});
%                     else
%                         if (length(varargin) > 1), error('error:  ncObj.put:  too many parameters: %s', arg); end
%                         obj.(name) = varargin{1};       % it's a property name.  varargin
                    end
                else
                    [pargs,cat] = ncObj.ncparse(varargin{:}, cat);
                    if (~iscell(pargs)), pargs = mat2cell(pargs,size(pargs,1)); end
                    if     (strcmp('Variables',cat)),   obj.putvar( pargs{:});
                    elseif (strcmp('Attributes',cat)),  obj.putatts(pargs{:});
                    elseif (strcmp('Dimensions',cat)),  obj.putdim( pargs{:});
                    elseif (strcmp('Groups',cat)),      obj.putgrp( pargs{:});
%                     else
%                         if (length(varargin) > 2), error('error:  ncObj.put:  too many parameters: %s', arg); end
%                         obj.(name) = varargin{2};       % it's a property name.
                    end
                end
            else
                    
                [ix, cat] = obj.findix(name, []);
                if (isempty_s(cat)), error('ncObj:put(%s):  cannot determine where to put data',name); end
%                if (length(spath)>1)            % Still more of path to parse.  descend to next level
                if (~isempty_s(spath))             % Still more of path to parse.  descend to next level
                    if (isempty(ix))
                        obj.put(sprintf('%s/%s',cat,name));
                        [ix,~] = obj.findix(name, cat);   % find index of name in category (insert if not present)
                    end    %?????
                    obj.(cat)(ix).put(spath, varargin{:});
%                 elseif (length(spath) == 1)     % path is only 1 long.  
%                     [~, cat] = obj.findix(name, []);
%                     if     (strcmp('Variables',cat)),   obj = obj.putvar( name, spath{1}, varargin{:});
%                     elseif (strcmp('Attributes',cat)),  obj = obj.putatts(name, spath{1}, varargin{:});
%                     elseif (strcmp('Dimensions',cat)),  obj = obj.putdim( name, spath{1}, varargin{:});
%                     elseif (strcmp('Groups',cat)),      obj = obj.putgrp( name, spath{1}, varargin{:});
% %                     else
% %                         if (length(varargin) > 1), error('error:  ncObj.put:  too many parameters: %s', arg); end
% %                         obj.(name) = varargin{1};       % it's a property name.  varargin{2} has value to put in prop                    
%                     end
                else
                    if (mod(length(varargin),2) ~= 0), error('error: ncObj.put(%s): argument list must be name/value pairs', name); end 
                    for i=1:2:length(varargin)
                        arg=varargin{i};
                        val=varargin{i+1};
                        obj.(cat)(ix).put(arg,val);
                    end
%                     [pargs,cat] = ncObj.ncparse(varargin{:}, cat);
%                     if (~iscell(pargs)), pargs = mat2cell(pargs,size(pargs,1)); end   % (needed to deal out pargs as separate values in putdim(...))
%                     if     (strcmp('Variables',cat)),   obj = obj.putvar( pargs{:});
%                     elseif (strcmp('Attributes',cat)),  obj = obj.putatts(pargs{:});
%                     elseif (strcmp('Dimensions',cat)),  obj = obj.putdim( pargs{:});
%                     elseif (strcmp('Groups',cat)),      obj = obj.putgrp( pargs{:});
% %                     else
% %                         if (length(varargin) > 2), error('error:  ncObj.put:  too many parameters: %s', arg); end
% %                         obj.(name) = varargin{2};       % it's a property name.  varargin{2} has value to put in prop                      
%                     end
                end
            end
        end

        
            % these put functions take either a single object, an array of objects, or a cell array of name/value pairs
            % which can be parsed into the correct type of object.
            
        function obj = update(obj, varargin)        % update an object, either with non-empty fields of a similar object, or from name,value pairs
            if (isempty(varargin)), return; end
            if (strcmp(class(obj),class(varargin{1})))
                upd = varargin{1};
                p=upd.properties();                
                for i=1:length(p)
                    v = upd.(p{i});
                    if (~isempty_s(v)), obj.(p{i}) = v; end
                end
            else
                for i=1:2:length(varargin)
                    arg=varargin{i};
                    val=varargin{i+1};
                    obj.(arg) = val;
                end
            end
        end
                        
                        
            
        function ix = putdim(obj, varargin)
            if (~isprop(obj,'Dimensions')), error('putdim:  obj %s does not have Dimensions field', class(obj)); end
            if (isa(varargin{1},'Dimension'))       % array of 1 or more Dimensions
                if (length(varargin) > 1)
                    ix = zeros(size(varargin));
                    for i=1:length(varargin)
                        ix(i) = obj.putdim(varargin{i});
                    end
                    return;
                else                    
                    mydim = varargin{1};
                end
            else
                mydim = Dimension(varargin{:});
            end
            if (isempty(obj.Dimensions))
                obj.Dimensions = mydim;
                ix = 1;
            else
                ix = obj.dimix(mydim.Name);
                if (isempty(ix)), ix = length(obj.Dimensions)+1; end
                obj.Dimensions(ix) = mydim;
            end
        end
        
        function ix = putvar(obj, varargin)
            if (isa(varargin{1},'Variable'))
                if (length(varargin)>1)
                    ix = zeros(size(varargin));
                    for i=1:length(varargin)
                        ix(i) = obj.putvar(varargin{i});
                    end
                    return;
                else
                    myvar = varargin{1};
                end
            else
                myvar = Variable(varargin{:});      % doesn't exist yet.  create it.
            end
            if (isempty(obj.Variables))
                obj.Variables = myvar;
                ix = 1;
            else
                ix = obj.varix(myvar.Name);
                if (isempty(ix)), ix = length(obj.Variables)+1; end
                obj.Variables(ix) = myvar;
            end
        end
        
        function ix  = putgrp(obj, varargin)
            obj.Format='netcdf4';       % if we're adding a group, then make the object a netcdf4 if it isn't already.
            if (isa(varargin{1},'Group'))       % we've been passed a group, along with things to put into.
                if (length(varargin)>1)
                    ix=obj.putgrp(varargin{1});
                    jx=obj.putgrp(varargin{2:end});
                    ix=cat(ix,jx);
                    return;
                else
                    mygrp = varargin{1};
                end
            else
                mygrp = Group(varargin{:});
            end
                        % if the groupname is empty, create a name for it, and make sure it isn't in the groups already.
            if (isempty_s(mygrp.Name) || strcmp(mygrp.Name,"/"))
                ng=length(obj.Groups);                    
                myname=sprintf("group%d",ng+1);
                while (~isempty(obj.grpix(myname)))
                    ng=ng+1;
                    myname=sprintf("group%d",ng+1);
                end
                mygrp.Name = myname;
            end                
            mygrp.Filename = obj.Filename;      % copy the filename and format in.
            mygrp.Format='netcdf4'; % force it to be a netcdf4.  don't just copy obj.Format;
            if (isempty(obj.Groups))
                obj.Groups = mygrp;
                ix = 1;
            else
                ix = obj.grpix(mygrp.Name);
                if (isempty(ix)), ix = length(obj.Groups)+1; end
                obj.Groups(ix) = mygrp;
            end
        end
        
        function putatt(obj, varargin)
            obj.putatts(varargin{:});
        end
        
        function putatts(obj, varargin)
            i=1;
            while (i <= length(varargin))
                if (isa(varargin{i},'Attribute'))
                    myatt = varargin{i};
                    i=i+1;
                else
                    myatt = Attribute(varargin{i:(i+1)});
                    i=i+2;
                end
                if (isempty(obj.Attributes))
                    obj.Attributes = myatt;
                else
                    ix = obj.attix(myatt.Name);
                    if (isempty(ix)), ix = length(obj.Attributes)+1; end
                    obj.Attributes(ix) = myatt;
                end
                if (isa(obj,'Variable') && strcmp(myatt.Name,'_FillValue'))
                    obj.FillValue = myatt.Value;
                end
            end
        end
        
        function tf = isnoncatprop(obj, name)
            fldix = find(strcmp(ncObj.classnames,class(obj)),1);
            flds = ncObj.noncats{fldix};
            tf = any(strcmp(name,flds));
        end
                    
        function data = loadvar(obj,src, vname, start, count, stride)
            % NOTE:  CAVEAT EMPTOR:  FillValues NOT replaced by nans.
            %           why?  because netcdf3 files use _FillValue (or some approx.), so can't guarantee that I
            %           replaced them.  
            
            
            if (~exist('stride','var')), stride=[]; end
            if (~exist('count', 'var')), count=[];  end
            if (~exist('start','var')), start=[]; end
            
            if (isa(obj,"Group"))
                datatype = obj.get(vname+"/Datatype");
            else
                datatype = [];
            end
                

            if (nargout > 0)
                data = cast(ncObj.readvar(src,vname,start,count,stride), datatype);
                obj.put(sprintf('%s/vdata',vname),data);
            else
                        % avoid duplicating data if not needed.  Helpful when reading a lot of data.
                obj.put(sprintf('%s/vdata',vname),cast(ncObj.readvar(src,vname,start,count,stride), datatype));
            end
        end

        function nc = NCstruct(obj)
            nc = struct;
%             fprintf('class: %s\n', class(obj));
%            disp(obj);
            flds = properties(obj);
            for i=1:length(flds)
                fld = flds{i};
%                 fprintf('\tfld: %s\n', fld);
                if (obj.isnoncatprop(fld))
                    if (~strcmp(fld,'vdata'))
                        if (isstring(obj.(fld)))
                            nc.(fld) = char(obj.(fld));
                        else
                            nc.(fld) = obj.(fld);                            
                        end
                    end
                else                    
                    if (isempty(obj.(fld)))
                        nc.(fld)=[];
                    else
                        for j=1:length(obj.(fld))
                            o=obj.(fld)(j);
%                             fprintf('fld %s j %d\n', fld,j);
                            nc.(fld)(j) = o.NCstruct();
                        end
                    end
                end
            end            
        end
        
        function ncDisp(obj)
            fprintf('class: %s\n', class(obj));
            disp(obj);
            flds = properties(obj);
            for i=1:length(flds)
                fld = flds{i};
%                 fprintf('\tfld: %s\n', fld);
                if (obj.isnoncatprop(fld))
%                    fprintf('\tnon-cat fld: %s\n', fld);
                else                    
                    if (isempty(obj.(fld)))
%                        fprintf('field %s is empty\n', fld);
                    else
                        for j=1:length(obj.(fld))
                            o=obj.(fld)(j);
                            fprintf('%s(%d):\n', fld,j);
                            if isa(obj.(fld)(j),'ncObj')
                                obj.(fld)(j).ncDisp();
                            else
                                if (isstring(o))
                                    fprintf('**********string*********: %s\n', o);
                                else
                                    disp(o);
                                end
                            end
                        end
                    end
                end
            end            
        end
    end

    methods(Static)
    
        function spath = split(spath)
            if (isempty_s(spath))
                return;
            end
            if (ischar(spath) || isstring(spath))
                spath = split(spath,'/');
            end
            spath = spath(~(cellfun(@isempty,spath)));
            spath = spath(~(cellfun(@strcmp,spath,repmat({'/'},size(spath)))));
        end
        
        
        function [ix, cat] = findcat(name)
            ix = find(strncmpi(ncObj.catnames,name,3),1);
            if (isempty_s(ix)), cat=[]; return; end
            cat=ncObj.catnames{ix};
        end
        
        function sz = mysize(v)
            if (isvector(v) || isrow(v))
                sz = length(v);
            else
                sz = size(v);
            end
        end

        function [argout, cat] = ncparse(varargin)
            % varargin:  array of ncObj's, or cell array.  Last element can be category name.
            % If category provided (or can be determined from class of varargin{1}), cat is set, and list is parsed into
            % array of objects of type cat.
            % If cat not know, returns a struct with fields set.
            
            if (isempty(varargin)), cat=[]; argout = []; return; end
            cix = ncObj.findcat(varargin{end});
            if (~isempty(cix))
                cat = ncObj.catnames{cix};
                varargin(end) = [];
                trycats=cix;
            else
                trycats=1:length(ncObj.catnames);
            end
            for i=trycats
                if (isa(varargin{1},ncObj.classnames{i}))
                    cat = ncObj.catnames{i};
                    argout = varargin;      % we were given an array of Variables, Attributes, Dimensions or Groups.
                    return;
                end
            end
            
                % OK.  now parse the arguments into specified category.
            nargs = length(varargin);
            if (~isempty(cix))
                flds = ncObj.fields{cix};
                if (cix == 1)   % Variable
                    argout = Variable;
                    if (mod(length(varargin),2)==1)         % first item is vdata
                        argout.vdata = varargin{1};
                        argout.Size = ncObj.mysize(argout.vdata);       % update the size field
                        istart = 2;
                    else
                        istart = 1;
                    end
                elseif (cix == 2)    % Attribute
                    argout = Attribute;
%                    if (~contains(varargin(1),flds))        % must be attribute name,length, pairs, instead of labeled fields
                    if (~any(strcmp(varargin(1),flds)))        % must be attribute name,length, pairs, instead of labeled fields
                        for i=1:2:nargs
                            ix = ceil(i/2);
                            argout(ix) = Attribute(varargin{i},varargin{i+1});
                        end
                        return;
                    end
                    istart = 1;
                elseif (cix == 3)   % Dimension
                    argout = Dimension;
%                     if (~contains(varargin(1),flds))        % not tagged name,length, unlim triplet
                    if (~any(strcmp(varargin(1),flds)))        % not tagged name,length, unlim triplet
                        argout = Dimension(varargin{1},varargin{2},false);
                        if (nargs > 2), argout.Unlimited = varargin{3}; end
                        return;
                    end
                    istart = 1;
                else                % group.
                    argout = Group;
                    istart = 1;
                end
                
                for i=istart:2:length(varargin)
                    arg = varargin{i};
                    val = varargin{i+1};
%                     if (~contains(arg,flds)), error('error parsing arguments:  %s not valid field for %s',arg, cat);  end
                    if (~any(strcmp(arg,flds))), error('error parsing arguments:  %s not valid field for %s',arg, cat);  end
                    argout.(arg) = val;
                end
                return;
            else                
                argout = struct;
                flds = ncObj.allfields;
                for i=1:2:length(varargin)
                    arg = varargin{i};
                    val = varargin{i+1};
%                     if (~contains(arg,flds)), error('error parsing arguments:  %s not valid field',arg);  end
                    if (~any(strcmp(arg,flds))), error('error parsing arguments:  %s not valid field',arg);  end
                    argout.(arg) = val;
                end
                return;
            
            end
                
        end
                   
        function s = trimcats(src)
            src = ncObj.split(src);
            len = length(src);
            s=[];
            for i=1:len
%                 if (~contains(src{i},ncObj.catnames))
                if (~any(strcmp(src{i},ncObj.catnames)))
                    s = sprintf('%s/%s',s,src{i});
                end
            end
        end
        
        function data = readvar(src, varName, start, count, stride)
                        % NOTE:  CAVEAT EMPTOR:  FillValues may NOT be replaced by nans.
            %           why?  because netcdf3 files use _FillValue (or some approx.), so can't guarantee that I
            %           replaced them.  

            varName=ncObj.trimcats(varName);
            if (~exist('start','var') || isempty_s(start))
                
                data = ncread(src, varName);
            elseif (~exist('count', 'var') || isempty_s(count))
                data = ncread(src, varName, start);                
            elseif (~exist('stride','var') || isempty_s(stride))               
                data = ncread(src, varName, start, count);
            else
                data = ncread(src, varName, start, count, stride);
            end
        end
        
        function mxmem = maxmem(mx)
            persistent mymaxmem;
            
            if (isempty(mymaxmem))
                mymaxmem = 16*1024*1024*1024;   % default to 16 GB
            end
            
            mxmem = mymaxmem;
            if (nargin)
                mymaxmem = mx;
            end
        end
        
        function tf = verbose(vbose)
            persistent myvbose;
            
            if (isempty(myvbose))
                myvbose = false; 
            end
            tf = myvbose;
            if (nargin)
                if (~islogical(vbose) && (~isnumeric(vbose)))
                    error("ncObj::verbose: error:  input must be logical or 0/1"); 
                else
                    myvbose = logical(vbose); 
                end
            end
        end            
    end
    
end

