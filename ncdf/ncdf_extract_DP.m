function DP = ncdf_extract_DP(ncobj, varargin)

    DP = [];
   if (any(strcmp(ncobj.Name,["DataParams", "DownscalingParams"])))
        DP = extract_DP(ncobj, varargin{:});
        if (isempty(DP.data_yrs))   % kludge.  data_yrs is empty.  should make sure data_yrs gets saved, Ian!
            DP.data_yrs = DP.data_final_yrs;
        end
   else
       for i=1:length(ncobj.Groups)
           DP = ncdf_extract_DP(ncobj.Groups(i), varargin{:});
           if (~isempty(DP))
               return;
           end
       end
   end               
end

function DP = extract_DP(ncobj, varargin)

    if (strcmp(ncobj.Name, "DownscalingParams"))
        DP = ARRM_V2_DownscalingParams(); % create a DownscalingParams object with default values.
    else
        DP = ARRM_V2_DataParams(); % create a DataParams object with default values.
    end

        % copy all DataParam fields from the group's attributes
    for i=1:length(ncobj.Attributes)
        nam = ncobj.Attributes(i).Name;
        val = ncobj.Attributes(i).Value;
        if (isprop(DP,nam))
            try         % in try/catch block in case nam is a dependent property.  
                        % there is a way to test if something is a dependent property, but it only works for
                        % non-hidden properties, so we just use a try/catch block here instead...
                        % Ref: https://www.mathworks.com/help/matlab/matlab_oop/getting-information-about-properties.html
                        
                            % special cases when extracting:
                if (strcmp(nam,"cdf_append_pts"))    % 2-D matrix for cdf_append_pts get stores as a 1-D vector as an attribute, and not ordered as we expect.
                    DP.(nam) = [val(1), val(3); val(2), val(4)];
                else        % all the rest should be ok.
                    DP.(nam) = val;
                end
            catch
            end
        end
    end
    
        % then apply any needed updates specified in varargin.
    DP = DP.update(varargin{:});
    
end
