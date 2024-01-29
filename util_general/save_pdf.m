function [pdfname, h] = save_pdf(h, orientation, pdfname, do_matfig)
%   pdfname = save_pdf(h, orientation, pdfname, fig_title)  saves figure in h as a PDF. 
%       h               figure handle.  If missing or empty, h == gcf.
%       Orientation     "landscape" or "portrait"  (default if empty: portrait) 
%       pdfname         if  empty or missing or fignum, figname is set to (main_progname)_fig_(fignum)

%       do_matfig       if present & not empty & true, saves figure as matlab .fig file as well, using same name w/ .fig extension
%
%   should make these kwd/value pairs, Ian!

    if (~exist("h","var") || isempty(h)),  h=gcf; end
    if (~exist("orientation", "var") || isempty(orientation))
        orientation = "portrait";
    end
    if (~exist("do_matfig","var") || isempty(do_matfig)), do_matfig = false; end
    
    if (strcmpi(orientation,"portrait"))
        set(h,'PaperOrientation','portrait');
        h.PaperPosition = [0,0,8.5,11];
    else        
        set(h,'PaperOrientation','landscape');
        h.PaperPosition = [0,0,11,8.5];
    end        

    if (~exist("pdfname","var") || isempty(pdfname))
        pdfname = sprintf("%s_fig_%04d.pdf", basename(get_progname(1), true), h.Number);
    elseif (isnumeric(pdfname))
        pdfname = sprintf("%s_fig_%04d.pdf", basename(get_progname(1), true), pdfname);
    else        
                % make sure extension is pdf
        [~,~,fext] = fileparts(pdfname);
        if (isempty(fext) || ~strcmpi(fext,".pdf"))
            pdfname = sprintf("%s.pdf", pdfname);
        end
    end

    pause(.25);
    saveas(h, pdfname, "pdf");

    if (do_matfig)
        [fpath,fname,~] = fileparts(pdfname);
        figname = fullfile(fpath,fname+".fig");
        savefig(h, figname, "compact");
    end
    
end