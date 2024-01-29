function nc = NCstruct(obj)
    if (~isa(obj, 'ncObj'))
        error('error:  NCstruct:  obj %s is not an ncObj',class(obj));
    else
        nc = obj.NCstruct();
    end
end
