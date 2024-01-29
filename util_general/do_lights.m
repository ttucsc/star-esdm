function do_lights(yes)
    if (~exist('yes','var') || yes)
        lightangle(120,30);
        lightangle(240,60);
        hl=light('Color',[.5,.5,.5]);
        lightangle(hl, 15,30);
        lightangle(190,75);
    else
        delete(findall(gcf, 'Type','light'))
    end
end