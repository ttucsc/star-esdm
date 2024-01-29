function grp = make_metadata_group(grpname, grp_description, nc, brief)
% grp = make_metadata_group(grpname, grp_description, nc)
%
%   returns a ncdf Group with metadata from source file nc_src (nc_src can be a ncdf object or a filename)
%   the src filename and directory are added as the first parameters.
%   if brief is true, group contains only the filename and directory.
%   if brief is false or missing, group contains the metadata (attributes) from nc_src as attributes of the group.
    
    if (isstring(nc) || ischar(nc))
        nc = ncdf(nc, "do_create", false);
    end
    if (~exist("brief","var")), brief = false; end
    
    descr_name = sprintf("Group_Description", grpname);
    try
%       descr = nc.get("/Attributes/Group_Description");
        descr = nc.getatt("Group_Desciption");
        grp_description = sprintf("%s\n%s", grp_description, descr);
    catch
    end

    grp = Group(grpname,descr_name, grp_description, "source_filename", basename(nc.Filename), "source_folder", dirname(nc.Filename),"brief",brief);
    if (~exist("brief","var") || ~brief)
        for i=1:length(nc.Attributes)
            grp.putatt(nc.Attributes(i));
        end
    end
end
