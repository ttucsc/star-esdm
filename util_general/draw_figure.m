function  draw_figure(nsides, fignum)

    if (~exist('fignum','var')), fignum=99; end
    figure(fignum);
    clf;
    
    A_len=1;
    B_len=A_len/sqrt(nsides);

    [pos_A, theta] = draw_fig([0,0], 0, A_len, nsides, 'r-', 2);
    
    for i=1:nsides
        pos_b1 = pos_A(i,:);
        theta_b1 = theta(i);
        draw_fig(pos_b1, theta_b1, B_len, nsides,'b-',1);
    end
end

function [pos, theta] = draw_fig(pos0, theta0, len, nsides, how, width)

    hold on;
    theta_corner = 2*pi/nsides;
    pos=repmat(pos0,nsides,1);
    theta=theta0*ones(nsides+1,1);
    for i=1:nsides
        pos(i+1,:) = pos(i,:) + len*[cos(theta(i)),sin(theta(i))];
        theta(i+1) = theta(i) + theta_corner;
    end
    
    for i=1:nsides
        x = [pos(i,1),pos(i+1,1)];
        y = [pos(i,2),pos(i+1,2)];
        plot(x,y,how, 'linewidth',width);
    end
    
    hold off;
    daspect([1,1,1]);

end

