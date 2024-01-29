function  plot_kdelines(fignum, myplines, mykplines, probs)

%     xvals = norminv(probs);
%     xvals(1)=-9;
%     xvals(end)=-9;
%     
    figure(fignum);
    
    ratios = mykplines./myplines;
    
    [yrlen,nprobs] = size(ratios);
    
    subplot(3,1,1);
    cla;
    plot(1:yrlen, myplines);
    subplot(3,1,2);
    cla;
    plot(1:yrlen, mykplines);
    subplot(3,1,3);
    cla;
    keepers=true(1,nprobs);
    keepers(ceil(nprobs/2))=false;  
    plot(1:yrlen, ratios(:,keepers));
    figure(fignum+1);
    clf;
    ii=0;
    for i=1:nprobs
        if (keepers(i))
            ii=ii+1;
            subplot(10,4,ii);
            plot(1:yrlen,ratios(:,i));
            title(sprintf('%.5f', probs(i)*100));
            ylim([1,1.3]);
            grid on;
            ii=ii+1;
            subplot(10,4,ii);
            plot(1:yrlen,myplines(:,i),1:yrlen, mykplines(:,i));
            title(sprintf('%.5f', probs(i)*100));
            grid on;
        end
    end
end

