function atts = ncdf_get_attributes(nc,group)
% atts = ncdf_get_attributes(nc,group)
%
%   returns struct with attributes from nc (global) or group attributes from nc/group
%       if group is not present in nc.Groups, returns and empty object.
%   
    atts = [];
    if (exist('group','var'))
        try
            g=nc.get(sprintf("Groups/%s",group));
        catch
            return;
        end
    else
        g = nc;
    end
    natts = length(g.Attributes);
    for i=1:natts
        att=g.Attributes(i);
        atts.(att.Name)=att.Value;
    end
end
