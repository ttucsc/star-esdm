function vals = jc_units_conversion(vals, fromUnits, toUnits)
% function to convert between units
%
%   ***This should be replaced with a full implementation of R's udunits2 routines.*** 
%       At present, there is no implementation in Matlab for the full set of possible units conversions.
%
%   For now, supports conversions between 
%       * precip and other measures (mm, cm, in, feet, meters and kg_m-2_s-1)
%       * temperature conversions (DegF, DegC, Kelvin)
%
%   arrays (xxx)strings give possible spellings for each set of units.
%   string comparisons are done as case-insensitive.
%
%   precip conversion from "kg_m-2_s-2" return units in whatevers-per-day, and vice versa. 
%       1 kg/m-2/s = 86400 mm/day (60 sec/min * 60 min/hour * 24 hrs/day) = 86.4 m/day
%
%   length/distance/height measures, which are simply scaling factors, are done by the product of 
%       lookups into to_meters and from_meters arrays.  No conversion is done if that product is 1.0 
%   temperature conversions need an offset and scaling, so are handled individually.
%
%   no conversion is done (and no error thrown) if the fromUnits and toUnits match, even if they are not valid unit
%   names.
%
%   12/4/2021 icsf:  added singular centimeter and millimeter to lists.  
%                    added code so when fromUnits and toUnits are from the same group, returns immediately. 
%____________________________________________________________________________


    if (strcmpi(fromUnits, toUnits)), return; end       % if strings are the same, do nothing.
    
%           to/from     m     cm     mm   inches       ft  kg_m-2_s-1  .1mm
    to_meters   = [     1.,  .01,  .001,   .0254,   .3048,   86.4,      .0001];  % to meters
    from_meters = [     1., 100., 1000., 1/.0254, 1/.3048, 1/86.4,     10000.];  % from meters to ... 
        
    
    DegCstrings   = ["C","DegC","Deg_C","Deg C","Degrees_C","Degrees C"];
    DegFstrings   = ["F","DegF","Deg_F","Deg F","Degrees_F","Degrees F"];
    Kelvinstrings = ["K","Kelvin","Kelvins"];
    
    Mstrings      = ["m","meter","meters"];                         % 1
    CMstrings     = ["cm","centimeter","centimeters"];              % 2
    MMstrings     = ["mm","millimeter","millimeters"];              % 3
    Inchstrings   = ["in","inch","inches"];                         % 4
    Footstrings   = ["ft","foot","feet"];                           % 5
    KGMSstrings   = ["kg m-2 s-1","kg_m-2_s-1","kgms"];             % 6
    MM1strings    = [".1mm","mm.1",".1millimeter",".1millimeters"]; % 7
    
    if     (any(strcmpi(fromUnits,DegCstrings)))
        if     (any(strcmpi(toUnits,  DegCstrings))),                             return;
        elseif (any(strcmpi(toUnits,  DegFstrings))), vals = vals * 9/5 + 32;     return;
        elseif (any(strcmpi(toUnits,Kelvinstrings))), vals = vals       + 273.15; return;
        else,  error('unknown output conversion string:  %s', toUnits);
        end
    
    elseif (any(strcmpi(fromUnits,DegFstrings)))
        if     (any(strcmpi(toUnits,  DegFstrings))),                                  return;
        elseif (any(strcmpi(toUnits,  DegCstrings))), vals = (vals-32) * 5/9;          return;
        elseif (any(strcmpi(toUnits,  DegFstrings))),                                  return;
        elseif (any(strcmpi(toUnits,Kelvinstrings))), vals = (vals-32) * 5/9 + 273.15; return;
        else,  error('unknown ouput conversion string:  %s', toUnits);
        end
        
    elseif (any(strcmpi(fromUnits,Kelvinstrings)))
        if     (any(strcmpi(toUnits,  Kelvinstrings))),                                  return;
        elseif (any(strcmpi(toUnits,  DegCstrings))), vals =  vals - 273.15;             return;
        elseif (any(strcmpi(toUnits,  DegFstrings))), vals = (vals - 273.15) * 9/5 + 32; return;
        elseif (any(strcmpi(toUnits,Kelvinstrings))),                                    return;
        else,  error('unknown ouput conversion string:  %s', toUnits);
        end
        
    else
        if     (any(strcmpi(fromUnits,   Mstrings))), ixfrom = 1;
        elseif (any(strcmpi(fromUnits,  CMstrings))), ixfrom = 2;
        elseif (any(strcmpi(fromUnits,  MMstrings))), ixfrom = 3;
        elseif (any(strcmpi(fromUnits,Inchstrings))), ixfrom = 4;
        elseif (any(strcmpi(fromUnits,Footstrings))), ixfrom = 5;
        elseif (any(strcmpi(fromUnits,KGMSstrings))), ixfrom = 6;
        elseif (any(strcmpi(fromUnits, MM1strings))), ixfrom = 7;
        else,  error('unknown input conversion string:  %s', fromUnits);
        end
            
        if     (any(strcmpi(toUnits,   Mstrings))), ixto = 1;
        elseif (any(strcmpi(toUnits,  CMstrings))), ixto = 2;
        elseif (any(strcmpi(toUnits,  MMstrings))), ixto = 3;
        elseif (any(strcmpi(toUnits,Inchstrings))), ixto = 4;
        elseif (any(strcmpi(toUnits,Footstrings))), ixto = 5;
        elseif (any(strcmpi(toUnits,KGMSstrings))), ixto = 6;
        elseif (any(strcmpi(toUnits, MM1strings))), ixto = 7;
        else,  error('unknown output conversion string:  %s', toUnits);                        
        end
        
        fctr = to_meters(ixfrom) * from_meters(ixto);
        if (fctr ~= 1.0)
            vals = vals * fctr; 
        end
        return;
    end  
end

