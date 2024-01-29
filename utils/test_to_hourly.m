
vars=["temp_F";"rh_F"];
vars2=["tasmin","tasmax"; "rhsmin","rhsmax"];
outvars=["tas"; "rhs"];
sources=["aashto","ncdc"];
nout = 0;
outnames = strings(8,1);
for iv=1:length(vars)
    for is=1:length(sources)
        ncHourly=sprintf("transportation/%s.stalist60.hourly.1979.2015.nc",vars(iv));
        ncMinvals=sprintf("transportation/%s.GFDL-ESM2G.r1i1p1.rcp85.%s.daily.1950.2100.nc", sources(is), vars2(iv,1));
        ncMaxvals=sprintf("transportation/%s.GFDL-ESM2G.r1i1p1.rcp85.%s.daily.1950.2100.nc", sources(is), vars2(iv,2));

        ncOutName=sprintf("transportation/%s.GFDL-ESM2G.r1i1p1.rcp85.%s.generated_hourly.1950.2100.nc", sources(is), outvars(iv));

        toDegF=true;
        useLocalTime = true;
        stns=[];
        yrRange_input=[];
        yrRange_output=[];
        outMinutes=[];
        inMinutes=[];

        ncComments=["Generating hourly data from observation data and downscaled model data",...
                    sprintf("using hourly %s from %s, ",vars(iv), ncHourly), ...
                    sprintf("and %s %s daily from %s and %s", sources(is), vars2(iv), ncMinvals, ncMaxvals)];
                
        ncAtts=[];
                        
        fnout = QC_daily_to_hourly(ncHourly, ncMinvals, ncMaxvals, ncOutName, ncComments, ncAtts,  toDegF, useLocalTime, stns, yrRange_input, yrRange_output, outMinutes, inMinutes);
        nout = nout+1;
        outnames(nout) = string(fnout);
    end
end
fprintf("output files: \n");
for i=1:nout
    fprintf("\t%2d %s\n", i, outnames(i));
end


        % ncAtts=["test_attribute1","this is a test attribute";  "test_attr2","second test attr"];
