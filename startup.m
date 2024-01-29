 setenv('PATH',[getenv('PATH'),':/opt/local/bin:/opt/local/sbin:/Users/iscottfl/bin']);

    % restore missing figure controls.
try %#ok
    if ~verLessThan('matlab','9.5')
        set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig));
        set(groot,'defaultAxesCreateFcn',  @(ax,~)set(ax.Toolbar,'Visible','off'));
    end
end
ARRM_V2_setpath();
