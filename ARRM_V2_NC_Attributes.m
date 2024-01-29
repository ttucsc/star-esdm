classdef ARRM_V2_NC_Attributes
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        title               = [];
        long_title          = [];
        institution         = [];
        comments            = [];
        references          = [];
        history             = [];
        hostname            = [];
        CreatedBy           = [];
        CreationDate        = [];
        Conventions         = "CF-1.6";
        Conventions_help    = "http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html";
        
        Unmatched           = [];
        runType             = 'unknown';
    end
    
    methods
        function obj = ARRM_V2_NC_Attributes(varargin)

                    % parse input for runType
            p = inputParser;
            p.StructExpand = true;
            p.KeepUnmatched = true;    
            addParameter(p,'runType',obj.runType);
            parse(p, varargin{:});
            runType = p.Results.runType;
            Unmatched = p.Unmatched;        % save rest of input params for later
            
                % apply any localized parameters if the function ARRM_V2_NC_Attribures_localized( ) exists.
            if (exist('ARRM_V2_NC_Attributes_TTUCSC.m','file'))
                obj = update(obj,ARRM_V2_NC_Attributes_TTUCSC(runType));
            end
            if (exist('ARRM_V2_NC_Attributes_localized.m','file'))
                obj = update(obj,ARRM_V2_NC_Attributes_localized(runType));
            end
                % now apply any command-line params (which are now in Unmatched);
            obj = update(obj, 'runType',runType, Unmatched);
            
        end
        
        function obj = update(obj, varargin)
            
            % Updates RunParams object.  Uses current values as defaults, and parses any key/value pairs
            % Call this when you want to set multiple parameters with a single function call.
            % To update a single parameter, just set the parameter directly.
    
            % I need to add verification to these still, Ian!

            p = inputParser;
            p.StructExpand = true;
            p.KeepUnmatched = true;
    
            addParameter(p,'title',obj.title);
            addParameter(p,'long_title',obj.long_title);
            addParameter(p,'institution',obj.institution);
            addParameter(p,'comments',obj.comments);
            addParameter(p,'references', obj.references);
            addParameter(p,'history',obj.history);
            addParameter(p, 'hostname',obj.hostname);
            addParameter(p, 'CreatedBy',obj.CreatedBy);
            addParameter(p, 'CreationDate',obj.CreationDate);
            addParameter(p, 'Conventions',obj.Conventions);
            addParameter(p, 'Conventions_help',obj.Conventions_help);      
            addParameter(p, 'Unmatched', []);
            addParameter(p, 'runType', obj.runType);
            
            parse(p, varargin{:});
            
            f=fields(p.Results);
            for i=1:length(f)
                if (~strcmpi(f{i}, 'Unmatched'))
                    obj.(f{i}) = p.Results.(f{i});
                end
            end
            obj.Unmatched = p.Unmatched;
        end
        
        function out = replace_keywords(obj, dataParams, kwds, instring)
            %   replaces keywords either in instring or in all the properties of obj.
            %   if instring is provided, then all keywords in instring are replaced, and new string returned.
            %   Otherwise, all instances of the keywords are replaced in all the objects' properties.
            %
            %   NOTE:  if replacing LAT1 & LATEND, or LON1 and LONEND, be sure to specify LATEND or LONEND before LAT or LON 
            %   
            if (~exist('kwds','var')), kwds=[]; end
            if (exist('instring','var'))
                out = dataParams.do_replace(kwds, instring);
            else
                props=obj.properties();
                for j=1:length(props)
                    prop=props{j};
                    if (ischar_s(prop))
                        obj.(prop) = dataParams.do_replace(kwds,obj.(prop));
                    end
                end
                out=obj;
            end            
            
        end
    end
end

