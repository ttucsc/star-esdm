function [h, final_pos, ok_flag] = set_figure_position(fignum, pos, vis, pauselen)
% pos = set_figure_position(fignum, pos, vis)
%
% function to reposition a figure (and possibly change its
% visibility).  Needed because when running interactively on a
% Linux system via X-windows, a pause is needed between actions.
% 
%   inputs:
%       fignum      figure number or handle
%       pos         new figure position. [lower_x, lower_y, size_x, size_y]
%                       if any value is pos is inf, then current position value is retained.
%                       For example, to resize but leave bottom left corner at current
%                       position, use [inf, inf, new_width, new_height]
%       vis         (optional) final visibility.  (true = default, visible, false = hidden)
%                       Note:  vis can be an array of 2 elements;  if length>1, then 
%                       the first sets the visibility while repositioning the window, the 2nd the final
%                       visibility.  use [true,true] on windows or macos to
%                       avoid blinking the figure off and back on.  Do NOT
%                       use [true,true] on Linux, or the window probably won't be
%                       repositioned or resized properly.
%       pauselen    (optional) pause length to use.  default: .25 for linux systems, 0 for macos & windows. 
%   returns:
%       h           handle to figure
%       pos         ending position.  
%       ok_flag     true if final position matches requested position,
%                   false if not
%                       If network lagtime is too long, then window may not
%                       be positioned or resized.  If that happens, try
%                       using a longer pauselen.

    if (~exist("vis","var") || isempty(vis))
        vis = [false, true]; 
    elseif (numel(vis) < 2)
        vis = [false, vis];
    end
    if (~exist("pauselen","var") || isempty(pauselen))
        on_linux = strncmpi(computer,"GLNX",4);
        if (on_linux)
            pauselen = .25;
        else
            pauselen = 0.1;
        end
    end
    
    h = figure(fignum);
    if (any(isinf(pos)))
        curpos = h.Position;
        pos(isinf(pos)) = curpos(isinf(pos));
    end
    
    h.Visible=vis(1);  
    pause(pauselen); 
    h.Position=pos; 
    pause(pauselen); 
    h.Visible=vis(end);
    
    if (nargout > 1)
        if (nargout > 2), pause(pauselen); end
        final_pos = h.Position;
        ok_flag = all(pos == final_pos);
    end
end
