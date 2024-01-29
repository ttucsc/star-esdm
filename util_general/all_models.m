function [models, unmatched] = all_models(varargin)
% all_models(...)  ...:  optional "model_set", "cmip[56]", "match","somestring" + "keepers", "excludes", "standard", "names", "match_station", "partial", "insensitive"
%       key/values:  
%           'excludes', [array or cell array] of models to exclude  (exact match only)
%           'standard', [true],false            if true, returns standard models (ones with good years).  if false, returns all
%           'match',    [list of models to match.] This does partial matching. 
%                                                   "gfdl" matches all models with "GFDL" or "gfdl" in the name.
%                                                   put "-" as 1st element to provide a list of models to exclude matches for.
%                                                   ["-","miroc"]  excludes all miroc models.
%           'names', 	true/[false].       if true, returns names only;  else returns table. 
%           'match_station' true/[false].   If true, will also match the work station, stations, stn or stns.
%           'partial', [true]/false         if true, matching is partial.  if false, matching must be exact
%           'insensitive',[true]/false      if true, matching is case-insensitive.
%           'keepers', keepers:             int or bool array of models to keep
%           
%   returns list of all models as array of strings.
%   if present, keepers can be numeric list or boolean list of which ones to keep.
%
%   If the file "model_list.csv" is present, it reads the table of models & info from the file.
%       ARRM_V2's csv file is a table including calender & #lats & #lons.
%   Otherwise, it uses the hard-coded list.
%
%   extra comment added.  

    p = inputParser;
    p.KeepUnmatched = false;
    addOptional(p,'keepers',[]);
    addParameter(p,'excludes', []);
    addParameter(p,'standard', true);
    addParameter(p,'match', []);
    addParameter(p,'names',false);
    addParameter(p,'match_station',false);
    addParameter(p,'partial',true);
    addParameter(p,'insensitive',true);
    addParameter(p,'model_set',"cmip5");
    
    parse(p,varargin{:});
    
    rp = p.Results;
    
    model_set   = rp.model_set;
    keepers     = rp.keepers;    
    excludes    = string(rp.excludes);
    matchlist   = string(rp.match);
    if (isrow(excludes)), excludes=excludes'; end
    
    if (rp.standard && strcmp(model_set,"cmip5"))
        excludes=[excludes;"CMCC-CMS";"CanCM4";"EC-EARTH";"HadCM3";"MIROC4h"];
    end

    if (strcmp(model_set, "cmip5"))
        models=readtable("cmip5_models.csv");
    else
        models=readtable("cmip6_models.csv");
    end
    stn = {"Stations","standard",1,1};
            
    if (rp.match_station)
        models = [models; stn];
    end
    
        % add an index row so user knows what rows the retained model info is from in the original file.
    models.ix = (1:size(models,1))';
    
    nmdls = size(models,1);
    if (islogical(keepers))
        keepers(nmdls+1:end)=[];
    elseif (isnumeric(keepers))
        keepers(keepers>nmdls)=[];
    elseif (ischars(keepers))
        keepnames=string(keepers);
        keepers=false(nmdls,1);
        for i=1:length(keepnames)
            ix=find(strcmpi(keepnames(i),models.model),1);
            if (~isempty(ix)), keepers(ix)=true; end
        end
    else
        error("invalid input: keepers");
    end
                
    if (~isempty(keepers))
        models=models(keepers,:);
    end
    if (~isempty(matchlist))
        if (rp.partial), matchlist = ["~",matchlist]; end
        [matches,unmatched] = find_matches(models.model, matchlist, rp.insensitive);
        models = models(matches,:);
    else
        unmatched = [];
    end
    
    
    if (~isempty(excludes))
        keepers = true(size(models,1),1);
        for i=1:length(excludes)
            ix=find(strcmpi(excludes(i), models.model));
            if (~isempty(ix))
                keepers(ix) = false;
            end
        end
        models = models(keepers,:);
    end
    
    if (rp.names)
        models = models.model;  % return names only.
    else
        models.ix=(1:size(models,1))';
        models=movevars(models,'ix','after','model');
    end
    
    
    
end
