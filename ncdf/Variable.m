classdef Variable < ncObj
    % netcdf Variable class.
    %Constructor:  variable_name, [values], [name,value...]
    %   where name,value is fieldname and value.
    %       Dimensions and Attributes can be cell arrays of Dimensions or Attributes.  
    %       Other fields must be a single value.
    %   Detailed explanation goes here
    
    properties
        Name;
        Dimensions;
        Size;
        Datatype;
        Attributes;
        ChunkSize;
        FillValue;
        DeflateLevel;
        Shuffle = false;
        vdata;
    end
    
    methods
        function obj = Variable(nm, varargin)
            % Variable(nm)          nm is a Variable to copy
            % Variable(nm, Variable) 1st varargin is a variable to be copied, and given name nm.
            % Variable(nm, name/value pairs)  nm is a variable to copy, name/value pairs to set/override properties.
            % Variable(nm, vals)    name, data values
            % Variable(nm, vals, name/value pairs)
            % Variable(nm, name/value pairs)
            if (nargin == 0), return; end
            if (isstruct(nm) )  % assume a netcdf Variable struct
                myvar = nm;               
                flds = fieldnames(myvar);
                for i=1:length(flds)
                    fld = flds{i};
                    if (isprop(obj,fld))         % copy only fields that are properties
                        if (obj.isnoncatprop(fld))
                            obj.(fld) = myvar.(fld);
                        elseif (strcmp(fld,'Dimensions'))
                            dims=myvar.Dimensions;
                            for j=1:length(dims), obj.put('Dimensions',Dimension(myvar.Dimensions(j))); end
                        elseif (strcmp(fld,'Attributes'))
                            atts=myvar.Attributes;
                            for j=1:length(atts)
                                obj.put('Attributes',Attribute(myvar.Attributes(j)));     % copy _FillValue attribute in to variable's FillValue field
                            end
                        end
                    end
                end
                obj.check_FillValue();      % make sure variable's FillValue property and _FillValue attribute agree.
                return;
            elseif (isa(nm, 'Variable'))
                props=properties(obj);
                for i=1:length(props)       % copy the properties one by one, in case name is a derived class of Variable
                    prop=props{i};
                    obj.(prop) = nm.(prop);
                end
            else
                obj.Name = char(nm); 
            end
            if (~isempty(varargin))
                if (mod(length(varargin),2)==1)
                    if (isa(varargin{1},'Variable'))
                        myvar = varargin{1};
                        props = properties(obj);
                        for i=1:length(props)
                            prop = props{i};
                            obj.(prop) = myvar.(prop);
                        end
                        obj.Name = nm;
                    else
                        obj.vdata = varargin{1};
                    end
                    istart=2;
                else
                    istart=1;
                end
                for i=istart:2:length(varargin)
                    arg = varargin{i};
                    val = varargin{i+1};
                    if (isprop(obj, arg))
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
            if (~isempty_s(obj.vdata))
                if (isempty(obj.Size))
                    obj.Size = ncObj.mysize(obj.vdata);
                elseif (obj.Size ~= ncObj.mysize(obj.vdata))
                    error('error:  Variable %s:  Size  doesn''t match vdata''s size', obj.Name);
                end
                if (isempty_s(obj.Datatype))
                    obj.Datatype = class(obj.vdata);
                end
            end
            obj.check_FillValue();      % make sure variable's FillValue property and _FillValue attribute agree.
        end
        
        function check_FillValue(obj)   
        % make sure Variable's FillValue field and _FillValue attribute agree (if _FillValue attribute exists).
        % If they both are non-empty and differ, FillValue takes precedence.  This should only happen if parent is a
        % netcdf3.
        % Also, use attribute "missing_value" if FillValue has been set to "disable", which happens when you create a
        % netcdf file using nofillmode==true 
            
                % Kludge here to fix a problem with netcdf files written with nofillmode=true.
            if (isprop(obj,"FillValue"))
                if (strcmp(obj.FillValue,'disable'))
                    fillv=[];
                    try
                        fillv = obj.gettattvalue("missing_value");
                    catch
                    end
                    if (isempty(fillv))
                        try
                            fillv = obj.getattvalue("_FillValue");
                        catch
                        end
                    end
                    if (isempty(fillv) && any(strcmp(obj.Datatype,["single","double"])))
                        fillv = cast(1e20,obj.Datatype);
                    end
                    obj.FillValue = fillv;
                end


                    % and make sure obj.FillValue matches attribute FillValue if present.
                try
                    fv = obj.getattvalue('Attributes/_FillValue');

                    if (isempty(obj.FillValue))             % FillValue is empty.  use _FillValue
                        obj.FillValue = fv.Value;
                    elseif (isempty(fv.Value) || fv.Value ~= obj.FillValue)      % _FillValue exists, and is different from FillValue.
                        fv.Value = obj.FillValue;
                    end
                catch
                end
            
            end
        end
        
        function ids = dimids(obj, grp)
            % returns dimids for variable.
            %  if grp is a netcdf ncid, then calls netcdf function. (0-based)
            %  otherwise if grp is a ncdf Group, then gets dimension indexes (1-based)
            ids=nan(length(obj.Dimensions),1);
            if (isa(grp, 'Group'))
                for i=1:length(ids)
                    ids(i)=grp.dimix(obj.Dimensions(i).Name);
                end
            else
                for i=1:length(ids)
                    ids(i)=netcdf.inqDimID(grp,obj.Dimensions(i).Name);
                end
            end
        end
        
        function ty = xtype(obj,format)
            % returns string for netcdf class for variable.
            if (strcmpi(format,'netcdf4'))
                if (strcmp(obj.Datatype,'logical'))
                    ty = 'NC_BOOL';
                elseif (strcmp(obj.Datatype,'int8'))
                    ty = 'NC_BYTE';
                elseif (strcmp(obj.Datatype,'uint8'))
                    ty = 'NC_UBYTE';
                elseif (strcmp(obj.Datatype,'char'))
                    ty = 'NC_CHAR';
                elseif (strcmp(obj.Datatype,'int16'))
                    ty = 'NC_SHORT';
                elseif (strcmp(obj.Datatype,'uint16'))
                    ty = 'NC_USHORT';
                elseif (strcmp(obj.Datatype,'int32'))
                    ty = 'NC_INT';
                elseif (strcmp(obj.Datatype,'uint32'))
                    ty = 'NC_UINT';
                elseif (strcmp(obj.Datatype,'int64'))
                    ty = 'NC_INT64';
                elseif (strcmp(obj.Datatype,'uint64'))
                    ty = 'NC_UINT64';
                elseif (strcmp(obj.Datatype,'single'))
                    ty = 'NC_FLOAT';
                elseif (strcmp(obj.Datatype,'double'))
                    ty = 'NC_DOUBLE';
                elseif (strcmp(obj.Datatype,'string'))      % store strings as char
                    ty = 'NC_CHAR';
                else
                    error('Variable:  Datatype not defined for variable %s', obj.Name);
                end
            else
                if (strcmp(obj.Datatype,'int8'))
                    ty = 'NC_BYTE';
                elseif (strcmp(obj.Datatype,'uint8'))
                    ty = 'NC_UBYTE';
                elseif (strcmp(obj.Datatype,'char'))
                    ty = 'NC_CHAR';
                elseif (strcmp(obj.Datatype,'int16'))
                    ty = 'NC_SHORT';
                elseif (strcmp(obj.Datatype,'uint16'))
                    ty = 'NC_SHORT';
                elseif (strcmp(obj.Datatype,'int32'))
                    ty = 'NC_INT';
                elseif (strcmp(obj.Datatype,'uint32'))
                    ty = 'NC_INT';
                elseif (strcmp(obj.Datatype,'int64'))
                    ty = 'NC_DOUBLE';
                elseif (strcmp(obj.Datatype,'uint64'))
                    ty = 'NC_DOUBLE';
                elseif (strcmp(obj.Datatype,'single'))
                    ty = 'NC_FLOAT';
                elseif (strcmp(obj.Datatype,'double'))
                    ty = 'NC_DOUBLE';
                else
                    error('Variable:  Datatype not defined for variable %s', obj.Name);
                end
            end
        end
        
        function dims = dimlist(obj)        % returns a cell array containing qualified names of all dimensions in group
            if (isempty(obj.Dimensions))
                dims=strings(0,0);
            else
                dims = string({obj.Dimensions.Name});
            end
        end
        
        function atts = attlist(obj)        % returns a cell array containing qualified names of all dimensions in group
            if (isempty(obj.Attributes))
                atts={};
            else
                atts = {obj.Attributes.Name};
            end
        end

        function [dimix, dimlen, unlim, dimname] = diminfo(obj, dname)
        % diminfo:  returns dimension index, length, unlimited setting and actual name for dname in ncobj
            % obj       is ncdf Variable.
            % dname     dimension name.
            %               For time, lat, lon or day_of_year, matches against likely alternatives, and returns actual
            %               dimension name in dimname.

            latnames=["lat","latitude","lats","latitudes"];
            lonnames=["lon","longitude","lons","longitudes"];
            timenames=["time","times"];
            doynames=["doy","day_of_year","dayofyear"];

            dlist = obj.dimlist();
            dname = string(dname);

            dimlen = [];
            dimname = strings(0,0);
            unlim = [];

            dimix=find(strcmpi(dname,dlist),1);        % see if variable exists.  
            if (~isempty(dimix))
                if (nargout==1), return; end
                dimname = obj.Dimensions(dimix).Name;
                dimlen  = obj.Dimensions(dimix).Length;
                unlim   = obj.Dimensions(dimix).Unlimited;
            else
                % no exact match.  look for alternate name for standard dimensions.
                if     (any(strcmpi(dname,latnames))), names = latnames;
                elseif (any(strcmpi(dname,lonnames))), names = lonnames;
                elseif (any(strcmpi(dname,timenames))), names = timenames;
                elseif (any(strcmpi(dname,doynames))), names = doynames;
                else 

                    return
                end
                for i=1:length(names)
                    dimix = find(strcmpi(names(i), dlist),1);
                    if (~isempty(dimix))
                        if (nargout==1), return; end
                        dimname = obj.Dimensions(dimix).Name;
                        dimlen  = obj.Dimensions(dimix).Length;
                        unlim   = obj.Dimensions(dimix).Unlimited;
                        return;
                    end
                end
            end       
        end
            
        function rotate_dims(obj, new_order, skip_vdata)
        % rotates dimensions on a variable (also swaps around size, chunksize, etc.).
        % if skip_vdata is present and true, then actual data is NOT rotated.  This is useful for variables that are too
        % large to rotate in memory.
        %   NOTE:  if skip_vdata, user must remember that vdata isn't rotated...
            if (~exist('skip_vdata','var') || isempty(skip_vdata)), skip_vdata = false; end
            if (isnumeric(new_order))
                permute_order = new_order;
            else
                nd = length(new_order);
                permute_order = nan(1,nd);
                dimnames = string({obj.Dimensions.Name});
                for i=1:nd
                    dimix = obj.diminfo(new_order(i));  % diminfo is smart enough tofind variable "longitude" if new_order(i) is "lon";
                    if (isempty(dimix)), error("rotate_dims:  Cannot find dimension %s in dimlist %s", new_order(i), vec2string(dimnames,'brackets','[]')); end
                    permute_order(i) = dimix;
                end
                if (any(isnan(permute_order))), error("error: can't match new_order to variable's dimensions for %s", obj.Name); end
            end
            obj.Dimensions(permute_order)=obj.Dimensions;
            obj.Size(permute_order) = obj.Size;
            if (~skip_vdata && ~isempty(obj.vdata))
                obj.vdata = permute(obj.vdata,permute_order);
            end
                % swap ChunkSize as well, if set.  (if string, could be "DEFAULT" or "CONTIGUOUS")
            if (~isempty(obj.ChunkSize) && isnumeric(obj.ChunkSize))
                obj.ChunkSize(permute_order) = obj.ChunkSize;
            end
        end
        
        
        function data = loadvar(obj, src, vname, start, count, stride, verbose)
            if (~exist('vname','var') || isempty_s(vname))
                vname=obj.Name;         % warning:  ncObj's read uses matlab's ncread function to read the data.
                                        % ncread uses a fully qualified name (i.e. path down to group) to read data from
                                        % a variable in a subgroup.  Therefore, this will read will only work here 
                                        % if variable is in ncdf, not in a subgroup!  
                                        % otherwise it must be fully qualified for ncread to work.
            end
            if (~exist('stride','var')), stride=[]; end
            if (~exist('count', 'var')), count=[];  end
            if (~exist('start','var')), start=[]; end
            if (~exist('verbose','var')), verbose=obj.verbose; end
            
            if (~isempty(count))
                mysize = prod(count,'all') * 8;
            else
                mysize = obj.total_nbytes();
            end
            if (verbose || mysize > ncObj.maxmem)
                                                            % if bigger than 1 Gbyte, read it in pieces.  Otherwise
                                                            % the underlying library uses a huge amount of space to read
                                                            % it in.
                
                obj.vdata = obj.read_by_parts(src, start, count, stride, verbose);
            else
                obj.vdata = ncObj.readvar(src, vname, start, count, stride);
                obj.Size = size(obj.vdata);
            end
            if (nargout > 0), data = obj.vdata; end
        end
        
        function data = read_by_parts(obj, src, start1, count1, stride1, verbose)
                % reads a variable's data in small parts to avoid matlab's excessive memory consumption when reading
                % data via ncread(...).  For numeric data, ncread reads the data as double, and can use 2x the amount of
                % memory temporarily.  This can cause problems for very large arrays.
                %
                %   NOTE:  this probably won't work right for char data, or huge integer arrays...
                %
            if (isempty(obj.Size)), error("error:  Variable.read_by_parts:  cannot read variable %s by parts because size is unknown (file: %s)", v.Name, src); end
            
            if (verbose), vbose = obj.verbose(verbose); end
            
                    % read 1 element so we know what kind of data we'll be getting.
            nd = length(obj.Size);
            
            if (exist('start1','var') && ~isempty(start1))
                start = start1;
            else
                start = ones(1,nd);
            end
%           istart = start(end)-1;      % one before where we want to start reading.
            if (exist('count1','var') && ~isempty(count1))
                sz = count1;
                nsteps = count1(end);
                count = count1;
%               count(end) = 1;
            else
                nsteps = obj.Size(end) - start(end)+1;
                count = obj.Size;                
%               count(end) = 1;
                sz = count;
                sz(end)=nsteps;
            end
            if (exist('stride1','var') && ~isempty(stride1))
                stride = stride1;
            else
                stride = ones(1,nd);
            end
                    % fix this, Ian!
            if (any(stride ~= 1)), error("sorry...can't do read_by_parts with a stride different from one.  Contact Ian to get this fixed!"); end
            
            if (obj.numeric())      % if data is numeric, see if data is of type single.  If so, allocate an array of singles.
                if (obj.element_size() == 4)    % can't do this for int32's though, because they might be packed, which would unpack to doubles.
                    data = zeros(sz, 'single'); % should improve this to handle int16 & int32's, Ian!
                else
                    data = zeros(sz); % otherwise, allocate an array of doubles
                end
            else
                s1=ones(1,nd);
                c1=ones(1,nd);
                d = ncread(src, obj.Name, s1, c1);    % 
                outtype = class(d);
                data = zeros(obj.Size, outtype);
            end
                
            if (nsteps < 100)
                ix = 1:(nsteps+1);
                nreads = nsteps;
            else
                ix = round(linspace(1,nsteps+1,101));
                nreads = 100;
            end
%           rdsize = prod(count);
            for i=1:nreads
                start(end)=start1(end) + ix(i)-1;
                count(end)=ix(i+1)-ix(i);
                ix1=ix(i);
                ix2=ix(i+1)-1;
                if (nd == 2)
                    data(:,ix1:ix2)     = ncread(src, obj.Name, start, count, stride);
                elseif (nd == 3)
                    data(:,:,ix1:ix2)   = ncread(src, obj.Name, start, count, stride);
                elseif (nd == 4)
                    data(:,:,:,ix1:ix2) = ncread(src, obj.Name, start, count, stride);
                else
                    error("read_by_parts for more than 4 dimensions not tested yet!");
%                     ix1 = (i-1)*rdsize + 1;
%                     ix2 =     i*rdsize;
%                     d = ncread(src, obj.Name, start, count);
%                     data(ix1:ix2) = d(:);
                end
                if (obj.verbose()), show_progress(i,nreads); end
            end
            
            if (verbose), obj.verbose(vbose); end
        end
        
        function write(obj, src, data, start, stride, permute_order)
                % if start or stride specified (i.e., not empty), use matlab's ncwrite to write a section of the variable's data.
                % Otherwise, writing entire variable's data.
                    % If > 1 GB of data, or if permute_order specified, write data slice by slice using write_by_parts.
                    % if < 1GB && no permute_order given, use ncwrite.
            if (~exist('data','var')), data = obj.vdata; end
            if (~exist('permute_order','var'))
                permute_order = []; 
            elseif (~isempty(permute_order) && (~isempty(start) || ~isempty(stride)))
                error("Variable.write():  %s:  error:  both permute_order and start/stride specified.  Can't handle that.", obj.Name);
            end
            
            w = whos('data');   % get info on variable data, so we can figure out how large it is.
            if (~isempty(permute_order) || w.bytes > ncObj.maxmem)      % if over maxmem, use write_by_parts(...), which is much slower than ncwrite(...)
                obj.write_by_parts(src, data, permute_order)                % but doesn't use as much memory.  Matlab's ncwrite(...)  uses 3x the size of the variable to write it out.
            elseif (~exist('start','var') || isempty(start))
                ncwrite(src, obj.Name, data);
            elseif (~exist('stride','var') || isempty(stride))
                ncwrite(src, obj.Name, data, start);
            else
                ncwrite(src, obj.Name, data, start, stride);
            end
        end
                
            
        function write_by_parts(obj, src, data, permute_order)                        

            if (~exist('data','var') || isempty(data)), data=obj.vdata; end
            
            mysize = size(data);
                % fix the size if we are missing singleton dimensions at the end.
                %   Matlab truncates trailing singleton dimensions when extracting slices of a multi-dimension array.
                %   but we need mysize to be the same size as permute_order.
            if (length(mysize) < length(obj.Size))
                nmissing = length(obj.Size) - length(mysize);
                mysize = [mysize, ones(1,nmissing)];
            end
            nd = ndims(data);
            start = ones(1,nd);
            if (exist('permute_order','var') && ~isempty(permute_order))
                mysize = mysize(permute_order);
                nwrites = mysize(end);
                if (nd > 6), error("error:  can't permute write_by_parts for %d dimensions", nd); end
                fetchdim = find(permute_order==nd,1);
                for i=1:nwrites 
                    start(end)=i;
    %                 if (nd == 2)
    %                     ncwrite(src, obj.Name, data(:,i),     start);
    %                 elseif (nd == 3)
    %                     ncwrite(src, obj.Name, data(:,:,i),   start, mysize);
    %                 elseif (nd == 4)
    %                     ncwrite(src, obj.Name, data(:,:,:,i), start, mysize);
    %                 else
                    dataslice = permute(extract_slice(data,fetchdim, i), permute_order);
                    ncwrite(src, obj.Name, dataslice, start);
                    show_progress(i,nwrites);
                end            
            else
                nwrites = mysize(end);
                for i=1:nwrites
                    try
                        start(end)=i;
%                         if (nd == 2)
%                             ncwrite(src, obj.Name, data(:,i),     start);
%                         elseif (nd == 3)
%                             ncwrite(src, obj.Name, data(:,:,i),   start);
%                         elseif (nd == 4)
%                             ncwrite(src, obj.Name, data(:,:,:,i), start);
%                         elseif (nd == 5)
%                             ncwrite(src, obj.Name, data(:,:,:,:,i), start);
%                         elseif (nd == 6)
%                             ncwrite(src, obj.Name, data(:,:,:,:,:,i), start);
%                         elseif (nd == 7)
%                             ncwrite(src, obj.Name, data(:,:,:,:,:,:,i), start);
%                         elseif (nd == 8)
%                             ncwrite(src, obj.Name, data(:,:,:,:,:,:,:,i), start);
%                         else
%                             error("error:  write_by_parts: %s:  can't write out more than 8 dimensions. Size:  %s", obj.Name, vec2string(obj.Size,'brackets','[]'));
%                         end
                        dataslice = extract_slice(data,nd, i);
                        ncwrite(src, obj.Name, dataslice, start);
                        if (obj.verbose()), show_progress(i,nwrites); end
                    catch me
                        fprintf("oops on loop %d of %d\n", i, nwrites);
                        rethrow(me);
                    end
                end
            end
        end
        
        function nbytes = element_size(obj)
            % returns # of bytes for a single element.
%             if (~isempty(obj.Attributes))
%                 has_offset = any(strcmp("add_offset",{obj.Attributes.Name}));
%                 has_scale  = any(strcmp("scale_factor",{obj.Attributes.Name}));
%             else
%                 has_offset = false;
%                 has_scale = false;
%             end
%           if     (any(strcmp(obj.Datatype,["int64","uint64","double"])) && ~(has_offset && has_scale)), nbytes = 8;
            if     (any(strcmp(obj.Datatype,["int64","uint64","double"]))),  nbytes = 8;
            elseif (any(strcmp(obj.Datatype,["single", "int32","uint32"]))), nbytes = 4;
            elseif (any(strcmp(obj.Datatype,["int16","uint16"]))), nbytes = 2;
            else
                nbytes = 1;
            end
        end
        
        function totbytes = total_nbytes(obj)
            % returns total # of bytes for object (whether or not its vdata is present)
            if (isempty(obj.Size))
                totbytes = [];
            else
                totbytes = prod(obj.Size) * obj.element_size;
            end
        end
        
        function tf = numeric(obj)
            % returns false for character or string data, true for everything else.
            tf = ~any(strcmp(obj.Datatype,["char","string","uchar"]));
        end           
        
        function myclone = clone(obj, parent, do_copy_data)
            
            myclone = Variable(obj.Name, 'Size', obj.Size, 'Datatype', obj.Datatype, 'ChunkSize', obj.ChunkSize, ...
                               'FillValue', obj.FillValue, 'DeflateLevel', obj.DeflateLevel, 'Shuffle', obj.Shuffle);
            if (~exist('parent','var')), parent = []; end
            if (~exist('do_copy_data','var')), do_copy_data = true; end
            if (do_copy_data)
                myclone.vdata = obj.vdata;
            end
                        % copy over all the attributes
            for i=1:length(obj.Attributes)
                att = obj.Attributes(i).clone();
                myclone.putatt(att);
            end
            
                        % copy over all the Dimensions.  Look to see if it's already been added to the parent.
            if (~isempty(parent))     % look for each dimension in parent group and use it i present,
                                            % otherwise, clone obj's Dimension, and add to parent group.
                for i=1:length(obj.Dimensions)                    
                    try
                        dim = parent.get(obj.Dimensions(i).Name);
                        myclone.putdim(dim);
                    catch
                        dim = obj.Dimensions(i).clone();
                        myclone.putdim(dim);
                        parent.putdim(dim);
                    end
                end
            else
                for i=1:length(obj.Dimensions)
                    dim = obj.Dimensions(i).clone();
                    myclone.putdim(dim);
                end            
            end
        end
        
         function cast(obj, typ)
             if (~isa(obj.vdata,typ) || ~strcmpi(obj.Datatype,typ))
                 obj.Datatype=typ;
                 obj.vdata=cast(obj.vdata,typ);
             end
         end
        
    end
    
end

