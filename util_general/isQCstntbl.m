function yn = isQCstntbl(tbl)
% function yn = isQCstntbl(tbl)
% returns true if tbl is output from reading a netcdf station table with QC_get_site_table(...)
    yn=false;
    if (istable(tbl))
        if (isfield(tbl.Properties.UserData,'isQCstntbl' )), yn=tbl.Properties.UserData.isQCstntbl;  end
        if (isfield(tbl.Properties.UserData,'isQC_stntbl')), yn=tbl.Properties.UserData.isQC_stntbl; end    % used to have an underscore in it, so need to check just in case it's an old file. 
    end
end

