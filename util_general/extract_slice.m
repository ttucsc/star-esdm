function data = extract_slice(indata, fetchdim, ix)
    % returns a hyperslice ix of indata, where the slice is taken along the dimension fetchdim
    % This is a brute force function for data of up to 8 dimensions, since matlab doesn't have a way of extracting from
    % an arbitrarily-dimensioned object.
    
    dimsize = size(indata);
    nd = ndims(indata);
    if (fetchdim > nd || fetchdim < 1)
        error("error:  extract_slice:  fetchdim %d not a valid dimension.  indata size: %s)", fetchdim, vec2string(dimsize, "brackets",'[]'))
    elseif (ix > dimsize(fetchdim) || ix < 1)
        error("error:  slice index %d out of range for dimension %d.  indata size size: %s)", ix, fetchdim, vec2string(dimsize, "brackets",'[]'))
    end

    switch fetchdim
        case 1
            data = indata(ix,:,:,:,:,:,:,:);
        case 2
            data = indata(:,ix,:,:,:,:,:,:);
        case 3
            data = indata(:,:,ix,:,:,:,:,:);
        case 4
            data = indata(:,:,:,ix,:,:,:,:);
        case 5
            data = indata(:,:,:,:,ix,:,:,:);
        case 6
            data = indata(:,:,:,:,:,ix,:,:);
        case 7
            data = indata(:,:,:,:,:,:,ix,:);
        case 8
            data = indata(:,:,:,:,:,:,:,ix);
    otherwise
        error("extract_dims:  indata has too many (%d) dimensions:  %s (8 max)", nd, vec2string(dimsize, "brackets",'[]'));
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
            

