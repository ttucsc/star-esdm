function yn = isNCtbl(tbl)
% function yn = isNCtbl(tbl)
% returns true if tbl is output from reading a gridded netcdf  table with ARRM_V2_get_grid_table(...)
    yn=false;
    if (istable(tbl))
        if (isfield(tbl.Properties.UserData,'isNCtbl' )), yn=tbl.Properties.UserData.isNCtbl;  end
    end
end

