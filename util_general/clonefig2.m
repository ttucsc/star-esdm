function hf2 = clonefig2(inFig,OutFig, position, show_fig, close_fig, clear_fig, pause_len)
% hf2 = clonefig2(inFig,OutFig, position, show_fig, close_fig, clear_fig, pause_len)
%
% this program copies a figure to another figure
% example: CloneFig(1,4) would copy Fig. 1 to Fig. 4
% Matt Fetterman, 2009
% pretty much taken from Matlab Technical solutions:
% http://www.mathworks.com/support/solutions/en/data/1-1UTBOL/?solution=1-1UTBOL
%
%  Modified  Ian Scott-Fleming 2023
%   using handle(...) instead of figure(...)
%   to allow cloning figures w/ Visible=false
%   This is helpful for environments where drawing to visible figures takes
%   time and can cause problems when working with multiple figures.
%   
%   Also added return value of handle to cloned figure.
%
%   Also added: position        position & size of figure on screen.
%               show_fig        if true, displays cloned figure on screen
%                                   default:  true.  False useful for drawing figs when running matlab over X
%               close_fig       if true, closes inFig after cloning it
%                                   default:  false. True useful when running either in batch mode or over X, 
%                                   because no graphical UI available for batch mode, and over X there are long delays.
%               clear_fig       if true, clears inFig when finished
%                                   default:  false.  True useful if reusing a hidden figure when running either
%                                   over X or in batch mode.
%               pause_len       time to pause after graphical actions.  Default:  0.1 secs.  
%                                   Necessary because graphics stuff seems to happen in a child thread, and 
%                                   especially when running over X, there can be long delays before a figure is actually
%                                   drawn on the screen
%                                   if pause_len is 2 elements, then (1) is after creating outFig, and (2) after copying, 
%                                       and before closing or clearing inFig.
%
% hf1=figure(inFigNum);
% hf2=figure(OutFigNum);

    hf1=handle(inFig);
    if (ishandle(OutFig))
        hf2=handle(OutFig);
    else
        hf2=figure(OutFig);
    end
    hf2.Visible = false;
    clf(hf2);

    if (~exist("position", "var") || isempty(position)),  position = hf1.Position; end
    if (~exist("show_fig", "var") || isempty(show_fig)),  show_fig = true;   end
    if (~exist("close_fig","var") || isempty(close_fig)), close_fig = false; end
    if (~exist("clear_fig","var") || isempty(clear_fig)), clear_fig = false; end
    if (~exist("pause_len","var") || isempty(pause_len)), pause_len = 0.1;   end

    hf2.Position = position;
    pause(pause_len(1));

    compCopy(hf1,hf2);

    pause(pause_len(end));

    if (show_fig), hf2.Visible = true; end
    if (close_fig)
        close(hf1); 
    elseif (clear_fig)
        clf(hf1);
    end
end

function compCopy(op, np)
%COMPCOPY copies a figure object represented by "op" and its % descendants to another figure "np" preserving the same hierarchy.

    ch = get(op, 'children');
    if ~isempty(ch)
        nh = copyobj(ch,np);
        for k = 1:length(ch)
            compCopy(ch(k),nh(k));
        end
    end
    return;
end