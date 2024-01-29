classdef ncdf < Group
    %ncdf   netcdf class
    %   class to encapsulate the contents of a netcdf file.
    %
    %   should make a child class ncdf_ic which inherits from this, and has the CSC-specific code in it.
    
    properties
%         Filename;
%         Format;
    end
    
    methods
        function obj = ncdf(nc,varargin)
            % ncdf(nc, name, value, name, value,...)
            % creates a netcdf object from nc.  nc can be a filename or a ncinfo struct.  
            % varargin are a series of name,value pairs defining contents of ncdf object, such as Filename, Format, etc.
            % leave nc empty to create a new nc and initialize from input params.
            % pass in "create_ok",true  to create a new file.  Otherwise Will throw exception if file doesn't exist.
            obj.Name='/';
            if (nargin == 0), return;  end
            create_ok = false;
                % if 1st param is a filename, read info from file into an ncinfo struct
            if ((ischar_s(nc) && ~isempty_s(nc)))
                if (~isfile(nc))
                    obj.Filename = nc;
                else
                    nc = ncinfo(nc);
                end
            end
                % if nc is a ncinfo struct or is a ncdf Group, copy info into obj
            if (isstruct(nc) || isa(nc,'Group'))
                if (isstruct(nc))
                    flds = fieldnames(nc);
                else
                    flds = properties(obj);
                end
                for i=1:length(flds)
                    fld = flds{i};
                    if (~isprop(obj,fld)), continue; end
                    if (strcmp(fld,'Dimensions'))
                        dims=nc.Dimensions;
                        for j=1:length(dims)
                            obj.putdim(Dimension(dims(j)));
                        end
                    elseif (strcmp(fld,'Variables'))
                        vars=nc.Variables;
                        for j=1:length(vars)
                            obj.putvar(Variable(vars(j)));
                        end
                    elseif (strcmp(fld,'Attributes'))
                        atts=nc.Attributes;
                        for j=1:length(atts)
                            obj.putatts(Attribute(atts(j)));
                        end
                    elseif (strcmp(fld,'Groups'))
                        grps=nc.Groups;
                        for j=1:length(grps)
                            obj.putgrp(Group(grps(j)));
                        end
                    else
                        obj.(fld) = nc.(fld);
                    end
                end                
            end
            for i=1:2:length(varargin)
                arg = varargin{i};
                val = varargin{i+1};
                if (strcmp(arg,'Name'))
                    obj.Name = char(val);
                elseif (strcmp(arg,'Filename'))
                    obj.Filename = char(val);
                elseif (strcmpi(arg,'Dimensions'))
                    obj.putdim(val);
                elseif (strcmpi(arg,'Variables'))
                    obj.putvar(val);
                elseif (strcmpi(arg,'Attributes'))
                    obj.putatts(val);
                elseif (strcmpi(arg,'Groups'))
                    obj.putgrp(val);
                elseif (strcmpi(arg,'Format'))
                    obj.Format = char(val);
                elseif (strncmpi(arg,'create',6) || strcmpi(arg,"create_ok"))
                    create_ok = logical(val);
                else
                    obj.putatts(arg,val);
                end
            end
            
            if (~create_ok && ~isfile(obj.Filename)), error("ncdf:error:  file does not exist: %s", obj.Filename); end  
            
            if (isempty_s(obj.Name)), obj.Name='/'; end
%            if (~isempty(obj.Groups))   % make sure the format is set to netcdf4 if there are subgroups.
            if (isempty(obj.Format))
                obj.Format='netcdf4';   % make them all netcdf4's for now.
            end
%            elseif (isempty_s(obj.Format))
%                obj.Format='64bit';
%            end

        end 
        
        function data = loadvar(obj, varName, start, count, stride, src, verbose)
            % data = readnc(obj, vname, start, count, stride, src)
            % reads a variable's data from obj.Filename, and updates the object's variable field.
            
            if (~exist('start','var')), start=[]; end
            if (~exist('count', 'var')), count=[];  end
            if (~exist('stride','var')), stride=[]; end
            if (~exist('src','var') || isempty_s(src)), src=obj.Filename; end
            if (~exist('verbose','var') || isempty(verbose)), verbose = false; end
            
            v=obj.get(varName);
            try
%               fprintf("loading %s\n", varName);
                data = v.loadvar(src, varName, start, count, stride, verbose); 
            catch me
                oops();
                rethrow(me);
            end
            
%             data = loadvar@ncObj(obj,src,varName,start,count,stride);
        end

        function data = readvar(obj, varName, start, count, stride, verbose)
            % data = readvar(obj, varName, start, count, stride)
            % reads a single variable from the file specified in obj.Filename
            % returns the data directly.  does not update the object itself.
            % NOTE:  CAVEAT EMPTOR:  FillValues may NOT be replaced by nans.
            %           why?  because netcdf3 files use _FillValue (or some approx.), so can't guarantee that they get replaced.
            
            if (~exist('stride','var')), stride=[]; end
            if (~exist('count', 'var')), count=[];  end
            if (~exist('start','var')), start=[]; end
            if (~exist('verbose','var')), verbose = false; end
            
                    % superceded.  Variable.readvar( ) now can handle everything.
%             if (verbose || (isempty(stride) && isempty(count) && isempty(start)))
%                 v = obj.get(varName);
%                 if (verbose || v.total_nbytes() > 1024*1024*1024)         % if bigger than 1 Gbyte, read it in pieces.  Otherwise
%                                                             % the underlying library uses a huge amount of space to read
%                                                             % it in.
%                     vbose = obj.verbose(verbose);
%                     data = v.read_by_parts(obj.Filename, start, count, stride);
%                     obj.verbose(vbose);
%                     return;
%                 end
%             end
% verbose=false;
            if (verbose)
                v=obj.getvar(varName);
                data = v.read_by_parts(obj.Filename, start, count, stride, verbose);
%               data = readvar@ncObj(obj.Filename, varName, start, count, stride);
            else                
                data = readvar@ncObj(obj.Filename, varName, start, count, stride);
            end
            
        end    
        
        function loadvars(obj, varNames,dims_only)
            % reads a select list of variables, or all variables, or dimensions only from file obj.Filename
            % does not return any variable data;  instead, returns the updated object
            % if varNames empty or not present, read all variables (or all dimensions only)
            % if varNames is logical, it is used as the dims_only flag 
            if (~exist('varNames','var') || isempty_s(varNames))
                varNames = obj.varlist();
            elseif (ischar(varNames))
                varNames = {varNames};
            elseif (isstring(varNames))
                varNames = cellstr(varNames);
            elseif (islogical(varNames))
                dims_only = varNames;
                varNames = obj.varlist();
            end
            if (~exist('dims_only','var')), dims_only = false; end

            if (dims_only)
                mydims = obj.dimlist();
            end
            
            for i=1:length(varNames)
                vname = varNames{i};
                if (~dims_only || any(strcmp(mydims, vname)))
%                     fprintf('reading %s\n', vname);
                    obj.loadvar(vname);
                end
            end
        end    
        
        function writeschema(obj, do_overwrite, varargin)
            % writes netcdf schema to file.
            %   do_overwrite    boolean.  If false or missing, will not overwrite existing file.
            %   varargin        optional string parameter keyword/value pairs
            %                       "Filename"  fname       specify output filename.  [default:  use obj.Filename]
            %                       "Dimensions" true/false; write Dimension variables to output file as well as schema
            %                       "Variables"  true/false; write all variables to output file as well as schema
            %                       "nofillmode" true/false; set nofillmode to false (initialize variables w/ NA-value) 
            %                                       or set nofillmode to true (do not initialize with fillValues)
            %
            %       This function only works when creating a new netcdf file.
            %       If adding to an existing netcdf file, add variable to ncdf object, along with dimension(s).
            %
            
            do_dims = false;
            do_vars = false;
            nofillmode = false;
            do_share = false;
            if (~exist('do_overwrite','var')), do_overwrite = false; end
            for i=1:2:length(varargin)
                arg = varargin{i};
                val = varargin{i+1};
                if (    strcmpi(arg,'filename')),       obj.Filename = char(val);
                elseif (strncmpi(arg,'dimension',9)),   do_dims = val;
                elseif (strncmpi(arg,'variable',8)),    do_vars = val;
                elseif (strcmpi(arg,'nofillmode')),     nofillmode=val;
                else, error('NCDF:UNKNOWN_ARG','error:  unknown argument: %s',arg);
                end
            end
            if (ischar_s(nofillmode))
                if     (strcmpi(nofillmode,'fill')),   nofillmode = false;
                elseif (strcmpi(nofillmode,'nofill')), nofillmode = true;
                else, error(NCDF:BADFILLMODE, 'error:  bad nofillmode, %s',nofillmode); 
                end
            end
            if (exist(obj.Filename,'dir')), error('error:  ncdf.writeschema(%s) exists and is a directory', obj.Filename); end
            if (exist(obj.Filename,'file'))
                if (~do_overwrite)
                    error('NCDF:FILE_EXISTS','error: ncdf.writeschema(%s): output file already exists', obj.Filename); 
                else
                    delete(obj.Filename);
                end
            end
            
%             if (nofillmode==false)
%                 ncinf = obj.NCstruct();
%                 ncwriteschema(obj.Filename, ncinf);
%             else
                icwriteschema(obj, do_overwrite, nofillmode, do_share);
%             end
            
            if (do_vars)
                obj.writevars();
            elseif (do_dims)
                obj.writevars([], true);
            end
        end
        
        function writeatt(obj, attName, attValue)
            % writes attribute and value to obj.Filename
            % KLUDGE FOR NOW.  THIS NEEDS TO PARSE attName for group/path info.  
            % writeatt(obj, attName, attValue):
            s=strsplit(string(attName),"/");
            
            if (strlength(s(1))==0)
                s=s(2:end);
            end
            if (length(s)==1)
                location = "/";
            else
                location = "/"+join(s(1:end-1),"/");
            end
            attName = s(end);
            if (isstring(attValue))           
                ncwriteatt(obj.Filename, char(location), char(attName), char(attValue));
            else
                ncwriteatt(obj.Filename, char(location), char(attName), attValue);
            end
        end
        
        function writevar(obj, varName, vardata, varargin)
            % writes variable data to obj.Filename .  If vardata not present or empty, takes data from the obj's vardata.
            % writevar(obj, varName, varargin):  possible versions
            %   varargin:   vardata), or 
            %               vardata, start), or 
            %               vardata, start, stride) or
            %               vardata, start, stride, permute_order)

            v = obj.get(varName);
    
                
            if (~exist('vardata','var') || isempty_s(vardata))
                v.write(obj.Filename, v.vdata, varargin{:});    % will write variable in parts if needed, either to permute array before writing or if > 1 GB of data.
            else
                v.write(obj.Filename, vardata, varargin{:});    % will write variable in parts if needed, either to permute array before writing or if > 1 GB of data.
            end
            
%             try
%                 fillval = v.FillValue;
%                 if (~isempty_s(fillval))        % I *think* ncwrite will take care of this automatically...
%                     vardata(isnan(vardata)) = fillval;
%                 end
%             catch
%                 fprintf('ncdf.writevar oops!\n');
%             end
%                 fprintf('%s.writevar(%s, %d %d)\n', obj.Name, varName, size(vardata,1),size(vardata,2));

%                ncwrite(obj.Filename, varName, vardata);
%             else
%                 if (length(varargin) == 1)
%                     start=varargin{1};
%  %                   try
% %                         fprintf("writing var %s, start is: ", varName);
% %                         disp(start);
%                         ncwrite(obj.Filename, varName, vardata, start);
%   %                  catch
%    %                     fprintf('writevar oops!\n');
%    %                 end
%                 elseif (length(varargin) == 2)
%                     start=varargin{1};
%                     stride = varargin{2};
% %                     fprintf("writing var %s, start, stride are: ", varName);
% %                     disp(start);
% %                     disp(stride);
%                     ncwrite(obj.Filename, varName, vardata, start, stride);
%                 else
%                     error('ncdf:writevar():  too many arguments');
%                 end
%             end            
        end
        
        function writevars(obj, varNames, do_dims, do_vars)
            % writes multiple variables to obj.Filename, stores data in obj.
            %   varNames            cell array, list of variables, or empty/not present
            %                           if empty or not present, get full list from obj
            %   do_dims, do_vars    boolean flags. if missing, do all vars in variables list
            %                           do_dims:  write dimension variables
            %                           do_vars:  write non-dimension variables
            %   
            %   Note:  to write only parts of multiple variables (using start and stride), use obj.writevar(...) to
            %           write each variable separately
            
            if (~exist('varNames','var') || isempty_s(varNames))
                varNames = obj.varlist();
            elseif (ischar(varNames))
                varNames = {varNames};
            elseif (isstring(varNames))
                varNames = cellstr(varNames);
            end
            if (~exist('do_dims','var')), do_dims = false; end
            if (~exist('do_vars','var')), do_vars = false; end
            if (~do_dims && ~do_vars)
                 do_both = true;  
            else
                do_both = false; 
            end

            myvars = obj.nondimvarlist();
            mydims = obj.dimlist();
            
            for i=1:length(varNames)
                vname = varNames{i};
                if (do_both || (do_dims && any(strcmp(mydims, vname))) || (do_vars && any(strcmp(myvars, vname))))
%                   vdata = obj.getvardata(vname);
                    v = obj.getvar(vname);
                    if (~isempty_s(v.vdata))
                        obj.writevar(vname);
                    end
                end
            end
        end
                
        function icwriteschema(obj, do_overwrite, nofillmode, do_share)
            if (~exist('do_overwrite','var') || isempty_s(do_overwrite)), do_overwrite = false; end
            if (~exist('nofillmode','var')   || isempty_s(nofillmode)),   nofillmode   = false; end
            if (~exist('do_share','var')     || isempty_s(do_share)),     do_share     = false; end
            
            if (do_overwrite)
                mode = netcdf.getConstant('CLOBBER');
            else
                mode = netcdf.getConstant('NOCLOBBER');
            end
            if (do_share), mode = bitor(mode, netcdf.getConstant('SHARE')); end
            do_enddef = true;
            if (strcmpi(obj.Format,'classic'))
                if (isempty(obj.Groups))
                    mode = bitor(mode, netcdf.getConstant('CLASSIC_MODEL'));
                else
                    mode = bitor(mode, netcdf.getConstant('netcdf4'));
                    obj.Format = "netcdf4";
                    do_enddef = false;
                end
            elseif (strncmpi(obj.Format,'64bit',5))
                mode = bitor(mode, netcdf.getConstant('64BIT_OFFSET'));
                do_enddef = true;
            elseif (strcmpi(obj.Format,'netcdf4_classic'))
                mode = bitor(mode, netcdf.getConstant('FORMAT_NETCDF4_CLASSIC'));                    
                do_enddef = false;
            else   
                mode = bitor(mode, netcdf.getConstant(obj.Format));
                do_enddef = false;
            end
                
            ncid = netcdf.create(obj.Filename,mode);
            
            icwriteschema@Group(obj,ncid, nofillmode);        % use Group's icwriteschema(...) function to recursively write out the schema 

%             if (nofillmode)     % this is supposed to keep it from initializing the file with nans, but it doesn't seem to work.  file always gets initialized...
% %                fprintf('calling initvars\n');
% %                t1=tic;
%                 oldfill = netcdf.setFill(ncid, 'NC_NOFILL');
% %                fprintf('in ncdf, oldfill was %g\n', oldfill);
%                 obj.initvars(ncid, ncid);           % write a single FillValue to all variables with a FillValue
%                                                     % this should allocate space for all the variables, but not
%                                                     % write FillValues into the entire space.
% %                newfill = netcdf.setFill(ncid, oldfill);
%                 netcdf.setFill(ncid, oldfill);
% %                fprintf('in ncdf, newfill was %g\n', newfill);
% %                elapsed = toc(t1);
% %                fprintf('%8.3f secs to initvars\n',elapsed);
%             else
            if (~nofillmode)
                netcdf.setFill(ncid, 'NC_FILL');
            end
            if (do_enddef)
                netcdf.endDef(ncid);
            end
            netcdf.sync(ncid);
            netcdf.close(ncid);
        end
        
%         function addvar(obj, varName, vardata, varargin)
%         % adds variable to netcdf file and ncdf object.
%         % varargin should specify Dimensions, Attributes, 
%             obj.putvar(varName, vardata, varargin)
%             obj.writevar(varName, vardata);
%             
%         end
        
    end
end


