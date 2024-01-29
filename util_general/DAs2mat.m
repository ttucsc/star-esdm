function m = DAs2mat(s, fieldname)
%converts a cell array of data to a 2-D column matrix.
% if s is a cell array of structs, extract field fieldname and concatenate into matrix.
% otherwise, s must be a cell array of vectors of matrices, all the same size.
% each vector from s is extended with NAs to the maximum length of any of the vectors.

    n=numel(s);
    if (isstruct(s{1}) || isa(s{1},'ARRM_V2_DisaggregateSignal') )                        % cell array of structs
        sz = size(s{1}.(fieldname));      % get max length of all the data
        sz=sz(sz~=1);
        if (isempty(sz)), sz=1; end
        nd = length(sz);
                            % note to self:  can handle any # of dims with a little creativity.
                            % reshape to column vector, join, then reshape aftewards.
                            % only question is what happens with singleton dimensions?
        if (nd > 3), error('DAs2mat: can''t handle a struct of matrics with %d dimensions',nd); end
        sz_out = [sz,n];
        m=nan(sz_out);
        for i=1:n                               % put each data set into the output matrix            
            if (nd==1)
                if (isrow(s{i}.(fieldname)))
                    m(:,i) = s{i}.(fieldname)'; 
                else
                    m(:,i) = s{i}.(fieldname);
                end
            elseif (nd==2)
                m(:,:,i) = s{i}.(fieldname);
            else
                m(:,:,:,i) = s{i}.(fieldname);
            end
        end
    else                                        % cell array of vectors or matrices
        sz = size(s{1});                  % get max length of all the data
        sz=sz(sz~=1);
        if (isempty(sz)), sz=1; end
        nd = length(sz);
        if (nd > 3), error('DAs2mat: can''t handle a struct of matrics with %d dimensions',nd); end
        sz_out = [sz,n];
        m=nan(sz_out);
        for i=1:n                               % put each data set into the output matrix
            if (nd==1)
                if (isrow(s{i}))
                    m(:,i) = s{i}'; 
                else
                    m(:,i) = s{i};
                end
            elseif (nd==2)
                m(:,:,i) = s{i};
            else
                m(:,:,i) = s{i};
            end
        end
    end
end

