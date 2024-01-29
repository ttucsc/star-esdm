function data = extract_keepers(indata, keepers, dimix)
    % returns a subset of indata, taken along dimension dimix.  keepers can be either a logical vector or a vector of
    % index values to keep.
    % This is a brute force function for data of up to 8 dimensions, since matlab doesn't have a way of extracting from
    % an arbitrarily-dimensioned object.
    %
    %   This function will work to extract data from matrices of up to 8 dimensions.
    %
    %   code relies on matlab not objecting to indata having fewer dimensions than specified in the statements below.
    %   If that changes in the future, this will need to switch to a set of nested switch statements similar to the
    %   comment code below.
    
    dimsize = size(indata);
    nd = ndims(indata);
    if (dimix > nd || dimix < 1)
        error("error:  dimix %d out of range.  indata size size: %s)", dimix, vec2string(dimsize, "brackets",'[]'))
    end

    switch dimix
        case 1
            data = indata(keepers,:,:,:,:,:,:,:);
        case 2
            data = indata(:,keepers,:,:,:,:,:,:);
        case 3
            data = indata(:,:,keepers,:,:,:,:,:);
        case 4
            data = indata(:,:,:,keepers,:,:,:,:);
        case 5
            data = indata(:,:,:,:,keepers,:,:,:);
        case 6
            data = indata(:,:,:,:,:,keepers,:,:);
        case 7
            data = indata(:,:,:,:,:,:,keepers,:);
        case 8
            data = indata(:,:,:,:,:,:,:,keepers);
    otherwise
        error("extract_keepers:  indata has too many (%d) dimensions:  %s (8 max)", nd, vec2string(dimsize, "brackets",'[]'));
    end                 
end

%   Use something like this if matlab changes how things work and the compact code above breaks for dimensions less than 8.
    
%     switch nd
%         case 2
%             switch fetchdim
%                 case 1
%                     data = indata(i,:);
%                 case 2
%                     data = indata(:,i);
%             end
%         case 3
%             switch fetchdim
%                 case 1
%                     data = indata(i,:,:);
%                 case 2
%                     data = indata(:,i,:);
%                 case 3
%                     data = indata(:,:,i,:);
%             end
%         case 4
%             switch fetchdim
%                 case 1
%                     data = indata(i,:,:,:);
%                 case 2
%                     data = indata(:,i,:,:);
%                 case 3
%                     data = indata(:,:,i,:);
%                 case 4
%                     data = indata(:,:,:,i);
%             end
%         case 5
%             switch fetchdim
%                 case 1
%                     data = indata(i,:,:,:,:);
%                 case 2
%                     data = indata(:,i,:,:,:);
%                 case 3
%                     data = indata(:,:,i,:,:);
%                 case 4
%                     data = indata(:,:,:,i,:);
%                 case 5
%                     data = indata(:,:,:,:,i);
%             end
%             
%         case 6
%             switch fetchdim
%                 case 1
%                     data = indata(i,:,:,:,:,:);
%                 case 2
%                     data = indata(:,i,:,:,:,:);
%                 case 3
%                     data = indata(:,:,i,:,:,:);
%                 case 4
%                     data = indata(:,:,:,i,:,:);
%                 case 5
%                     data = indata(:,:,:,:,i,:);
%                 case 6
%                     data = indata(:,:,:,:,:,i);
%             end          
%         otherwise
%             error("extract_dims:  indata has too many (%d) dimensions:  %s (6 max)", nd, vec2string(dimsize, "brackets",'[]'));
%     end
            

