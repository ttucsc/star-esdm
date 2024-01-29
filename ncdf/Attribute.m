classdef Attribute < ncObj
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name;
        Value;
    end
    
    methods
        function obj = Attribute(name,val)
            if (exist('val','var'))
                if (isstring(val))
                    val=char(val);              % underlying matlab netcdf stuff doesn't like strings...
                elseif(islogical(val))
                    val=uint8(val);             % underlying matlab netcdf stuff doesn't like logical either 
                end     
            end
            if (nargin == 0),  return; end
            if (isstruct(name))     % assume a netcdf Attribute struct.
                flds=fieldnames(name);
                for i=1:length(flds)
                    fld=flds{i};
                    if (isprop(obj,fld))        % copy only Attribute fields
                        obj.(fld) = name.(fld);
                    end
                end
                % obj.Name=name.Name; 
                % obj.Value = name.Value;
            elseif (isa(name,'Attribute'))
                obj.Name = name.Name;
                obj.Value = name.Value;
            elseif (nargin == 1), obj.Name=char(name); obj.Value = [];  
            elseif (nargin == 2), obj.Name=char(name); obj.Value = val; 
            else
                error('Attribute():  too many arguments');
            end
        end
%         function tf = isempty_s(obj)
%             tf = length(obj) == 1 && isempty_s(obj.Name) && isempty_s(obj.Value);
%         end
        function myclone = clone(obj)
            myclone = Attribute(obj.Name, obj.Value);
        end
        
    end
   
end

