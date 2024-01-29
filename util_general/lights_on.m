function lights_on(n)
    lightangle(120,30);
    lightangle(240,60);
    hl=light('Color',[.5,.5,.5]);
    lightangle(hl, 15,30);
   if (exist('n','var') && n)
        hl2=light('Color',[.5,.5,.5]);
        lightangle(hl2,-15,85);
        hl3=light('Color',[.75,.75,.75]);
        lightangle(hl3,-150,85);
   end
%        hl2=light('Color',[.75,.75,.75]);
end
   