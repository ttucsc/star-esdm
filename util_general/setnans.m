function data = setnans(data,vals)
    data(isnan(data))=vals;
end

