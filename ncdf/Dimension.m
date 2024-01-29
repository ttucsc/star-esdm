classdef Dimension < ncObj
    % Object representing a netcdf dimension
    % 
    
    properties
        Name
        Length
        Unlimited
    end
    
    methods
        function obj = Dimension(nm, len, unlim)
            % Dimension(nm)         nm (name of dimension)
            % Dimension(nm, len)    nm is new name, len is dimension
            % Dimension(nm, len, unlim)  name, length, unlimited.  (length & unlimited optional, default to 0 & false)
            if (nargin == 0), return; end
            if (isstruct(nm))       % assume nma netcdf Dimension struct.               
                flds = fieldnames(nm);
                for i=1:length(flds)
                    fld = flds{i};
                    if (isprop(obj,fld))         % copy only Attribute fields
                        obj.(fld) = nm.(fld);
                    end
                end
            elseif (isa(nm,'Dimension') || isa(len,'Dimension'))
                props=properties(obj);
                if (isa(nm, 'Dimension'))
                    mydim=nm; 
                    myname = mydim.Name; 
                else
                    mydim = len; 
                    myname = nm; 
                end
                for i=1:length(props)       % copy the properties one by one, in case name is a derived class of Dimension
                    prop=props{i};
                    obj.(prop) = mydim.(prop);
                end
                obj.Name = char(myname);
                return;
            else
                obj.Name=char(nm);
                obj.Length=0;
                obj.Unlimited=false;
            end
            if (nargin == 1), return; end
            if (~isinf(len) && (~isnumeric(len) || mod(len,1)~= 0)), error('error:  Dimension:  Length must be integer'); end
            obj.Length=len;
            obj.Unlimited=isinf(len);
            if (nargin == 2), return; end
            if (~islogical(unlim) && unlim ~= 0 && unlim ~= 1), error('error:  Dimension:  Unlimited must be logical'); end
            obj.Unlimited = logical(unlim);
        end
        
        function myclone = clone(obj)
            myclone = Dimension(obj.Name, obj.Length, obj.Unlimited);
        end

    end
end

