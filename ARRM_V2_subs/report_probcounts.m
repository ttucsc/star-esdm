function report_probcounts(probs, probcounts_base, probcounts_rolling, totcount_base, totcounts_rolling, lbl, lbl2, rolling_yrs)
    tcount_base = totcount_base;            % non-rounded
    tcounts_rolling = totcounts_rolling;    % non-rounded
    
    nsets = length(totcounts_rolling);
    totcount_base = round(totcount_base);
    totcounts_rolling = round(totcounts_rolling);
    totcounts_tot = round(sum(tcounts_rolling));
    fprintf('\n%s Probability Counts %s\n', lbl, lbl2);
    nprobs = length(probs);
    pcountsb = round(probcounts_base);
    pcountsr = round(probcounts_rolling);
    mid = ceil(nprobs/2);
    pcountsb(mid+1:end) = round(tcount_base-probcounts_base(mid+1:end));
    probFracsb = probcounts_base/tcount_base;
    probpctsb = probFracsb*100.0;
    if (nsets > 0)
        probFracsr=nan(nprobs,nsets);
        for i=1:nsets
            pcountsr(mid+1:end,i) = round(tcounts_rolling(i)-probcounts_rolling(mid+1:end,i));
            probFracsr(:,i) = probcounts_rolling(:,i)/tcounts_rolling(i);
        end
        probpctsr = probFracsr*100.0;
    end
    
    fprintf('          ');
    fprintf('%8.4f ', probs*100.0);
    fprintf('%8.4f ', 100);
    fprintf('\n');
    fprintf('          ');
    fprintf('%8d ', pcountsb);
    fprintf('%8d ', totcount_base);
    fprintf('\n');
    fprintf('          ');
    fprintf(' %7.3f%%', probpctsb);
    fprintf('  100.00');
    fprintf('\n\n')
    
    if (nsets > 0)
        fprintf('          ');
        fprintf('%8.4f ', probs*100.0);
        fprintf('%8.4f ', 100);
        fprintf('\n');
        for i=1:nsets
            fprintf('   %4d   ', rolling_yrs(i));
            fprintf('%8d ', pcountsr(:,i));
            fprintf('%8d ', totcounts_rolling(i));
            fprintf('\n');
            fprintf('          ');
            fprintf('%7.3f%% ', probpctsr(:,i));
            fprintf('  100.00');
            fprintf('\n')
        end

        pcntst     = sum(probcounts_rolling,2);
        pcountstot = round(pcntst);
        pcountstot(mid+1:end) = max(0, round(sum(tcounts_rolling)-pcntst(mid+1:end)));
        probFracstot = sum(probcounts_rolling,2)/sum(tcounts_rolling);
        probpctstot = probFracstot*100.0;

        fprintf('\n   total  ');
        fprintf('%8.4f ', probs*100.0);
        fprintf('%8.4f ', 100);
        fprintf('\n');
        fprintf('          ');
        fprintf('%8d ', pcountstot);
        fprintf('%8d ', totcounts_tot);
        fprintf('\n');
        fprintf('          ');
        fprintf('%7.3f%% ', probpctstot);
        fprintf('  100.00');
        fprintf('\n')
    end    
        
end
