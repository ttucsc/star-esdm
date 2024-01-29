function pdfname = safe_pdf(h, orientation, pdfname, fig_title)
% saves figure in h as a PDF.  
%       Orientation:  "landscape" or "portrait"  (default if empty: portrait) 
%       figname:  if  empty or missing, figname is set to (main_progname)_fig_(fignum)
%       fig_title if present, is added at top of figure
%           useful for figures with multiple plots.

    if (~exist("h","var") || isempty(h))
        h=gcf;
    end
    v=h.Visible;
    figure(h);
    h.Visible = v;

    if (~exist("orientation", "var") || isempty(orientation))
        orientation = "portrait";
    end
    if (~exist("pdfname","var") || isempty(pdfname) || strlength(pdfname)==0)
        pdfname = sprintf("%s_fig_%04d.pdf", get_progname(), h.Number);
    end

    if (strcmpi(orientation,"portrait"))
        set(h,'PaperOrientation','portrait');
        h.PaperPosition = [0,0,11,8.5];
    else        
        set(h,'PaperOrientation','landscape');
        h.PaperPosition = [0,0,8.5,11];
    end

    if (exist("fig_title", "var") && ~isempty(fig_title))
        ha = axes('Position',[0 0 1 1],'Visible','off');
        text(ha,.5,.975,fig_title, "fontsize", 16,"HorizontalAlignment","center");
    end

    saveas(h, pdfname, "pdf");
end