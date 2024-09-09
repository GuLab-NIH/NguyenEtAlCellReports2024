
% load data
% contact Lead Author for data



%% FOR avr3 - vvr1/vvr2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot % overlap compared to others
count = 1;
figure; hold on
clear xx;
for vr = 1:3
    switch vr
        case 1
            load('extra_31_4.mat');
            clr = 'm';
        case 2
            load('extra_31d_3.mat');
            clr = 'k';
        case 3
            load('extra_31b_6.mat');
            clr = 'g';
    end
    if vr == 2
        xx{count} = cat(1, prcOverlap{2});
        violin2( xx{count}, 'x', count, 'facecolor', clr, 'facealpha', 0.3);
        count = count+1;
        
    else
        
        xx{count} = cat(1, prcOverlap{1,2});
        violin2( xx{count}, 'x', count, 'facecolor', clr, 'facealpha', 0.3);
        count = count+1;
        
        xx{count} = cat(1, prcOverlap{2,3});
        violin2( xx{count}, 'x', count, 'facecolor', clr, 'facealpha', 0.3);
        count = count+1;
    end
end
xlim([.5 5.5]);
ylabel('% overlap');
set(gca, 'xtick', 1:5, 'xticklabel', {'avr1-2', 'avr2-3', 'avr3-vvr1', 'vvr1-2', 'vvr2-3'}, 'xticklabelrotation', 30);
set(gcf, 'position', [222 222 345 258]);


[~,p] = ttest2( xx{1}, xx{2} ); disp(['1-2: ' num2str(p)]);
[~,p] = ttest2( xx{4}, xx{5} ); disp(['4-5: ' num2str(p)]);
[~,p] = ttest2( xx{1}, xx{3} ); disp(['1-3: ' num2str(p)]);
[~,p] = ttest2( xx{2}, xx{3} ); disp(['2-3: ' num2str(p)]);
[~,p] = ttest2( xx{4}, xx{3} ); disp(['4-3: ' num2str(p)]);
[~,p] = ttest2( xx{5}, xx{3} ); disp(['5-3: ' num2str(p)]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot acitivty correlation

count = 1;
figure; hold on
clear xx;
for vr = 1:3
    switch vr
        case 1
            load('extra_31_4.mat');
            clr = 'm';
        case 2
            load('extra_31d_3.mat');
            clr = 'k';
        case 3
            load('extra_31b_6.mat');
            clr = 'g';
    end
    
    
    if vr == 2
        xx{count} = cat(1, activityCorr{2});
        yy = cat(1, activityCorr_shuffled{2});
        
        violin2( xx{count}, 'x', count, 'facecolor', clr, 'facealpha', 0.3);

        [~,p] = kstest2( xx{count}, yy); disp(p);
        count = count+1;
        
    else
        
        xx{count} = cat(1, activityCorr{1,2});
        yy = cat(1, activityCorr_shuffled{1,2});
        
        violin2( xx{count}, 'x', count, 'facecolor', clr, 'facealpha', 0.3);
        count = count+1;
        
        xx{count} = cat(1, activityCorr{2,3});
        yy = cat(1, activityCorr_shuffled{2,3});
        
        violin2( xx{count}, 'x', count, 'facecolor', clr, 'facealpha', 0.3);
        count = count+1;
    end
end
xlim([.5 5.5]);
ylabel('activity correlation');
set(gca, 'xtick', 1:5, 'xticklabel', {'avr1-2', 'avr2-3', 'avr3-vvr1', 'vvr1-2', 'vvr2-3'}, 'xticklabelrotation', 30);
set(gcf, 'position', [444 222 444 258]);

[~,p] = ttest2( xx{1}, xx{2} ); disp(['1-2: ' num2str(p)]);
[~,p] = ttest2( xx{4}, xx{5} ); disp(['4-5: ' num2str(p)]);
[~,p] = ttest2( xx{1}, xx{3} ); disp(['1-3: ' num2str(p)]);
[~,p] = ttest2( xx{2}, xx{3} ); disp(['2-3: ' num2str(p)]);
[~,p] = ttest2( xx{4}, xx{3} ); disp(['4-3: ' num2str(p)]);
[~,p] = ttest2( xx{5}, xx{3} ); disp(['5-3: ' num2str(p)]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot cue correlation


count = 1;
figure; hold on
for vr = 1:3
    switch vr
        case 1
            load('extra_31_4.mat');
            clr = 'm';
        case 2
            load('extra_31d_3.mat');
            clr = 'k';
        case 3
            load('extra_31b_6.mat');
            clr = 'g';
    end
    
    
    if vr == 2
        cueCorrs = cat(1, cueCorr{2});
        dt1 = cueCorrs(:,1);
        dt2 = cueCorrs(:,2);

        edges = prctile(dt1, 0:10:100);
        [~, ~, binIdx] = histcounts( dt1, edges );
        binnedDt1 = nan(2, length(edges)-1);
        binnedDt2 = nan(2, length(edges)-1);

        for ii = 1:(length(edges)-1)
            kk = find(binIdx == ii);
            binnedDt1(1, ii) = nanmean( dt1(kk) );
            binnedDt1(2, ii) = nanstd( dt1(kk) ) / sqrt(length(kk));

            binnedDt2(1, ii) = nanmean( dt2(kk) );
            binnedDt2(2, ii) = nanstd( dt2(kk) ) / sqrt(length(kk));
        end
        subplot(2,3,count); hold on
        errorbar(binnedDt1(1,:), binnedDt2(1,:), binnedDt2(2,:), 'k');
        iidx = intersect( find(~isnan(dt1)), find(~isnan(dt2)) );
        [r,p] = corr(dt1(iidx), dt2(iidx));
        xlabel(['r = ' num2str(r) '; p = ' num2str(p)]);
        count = count+1;
        
    else
        
        cueCorrs = cat(1, cueCorr{1,2});
        dt1 = cueCorrs(:,1);
        dt2 = cueCorrs(:,2);

        edges = prctile(dt1, 0:10:100);
        [~, ~, binIdx] = histcounts( dt1, edges );
        binnedDt1 = nan(2, length(edges)-1);
        binnedDt2 = nan(2, length(edges)-1);

        for ii = 1:(length(edges)-1)
            kk = find(binIdx == ii);
            binnedDt1(1, ii) = nanmean( dt1(kk) );
            binnedDt1(2, ii) = nanstd( dt1(kk) ) / sqrt(length(kk));

            binnedDt2(1, ii) = nanmean( dt2(kk) );
            binnedDt2(2, ii) = nanstd( dt2(kk) ) / sqrt(length(kk));
        end
        subplot(2,3,count); hold on
        errorbar(binnedDt1(1,:), binnedDt2(1,:), binnedDt2(2,:), 'k');
        iidx = intersect( find(~isnan(dt1)), find(~isnan(dt2)) );
        [r,p] = corr(dt1(iidx), dt2(iidx));
        xlabel(['r = ' num2str(r) '; p = ' num2str(p)]);
        count = count+1;
        
        cueCorrs = cat(1, cueCorr{2,3});
        dt1 = cueCorrs(:,1);
        dt2 = cueCorrs(:,2);

        edges = prctile(dt1, 0:10:100);
        [~, ~, binIdx] = histcounts( dt1, edges );
        binnedDt1 = nan(2, length(edges)-1);
        binnedDt2 = nan(2, length(edges)-1);

        for ii = 1:(length(edges)-1)
            kk = find(binIdx == ii);
            binnedDt1(1, ii) = nanmean( dt1(kk) );
            binnedDt1(2, ii) = nanstd( dt1(kk) ) / sqrt(length(kk));

            binnedDt2(1, ii) = nanmean( dt2(kk) );
            binnedDt2(2, ii) = nanstd( dt2(kk) ) / sqrt(length(kk));
        end
        subplot(2,3,count); hold on
        errorbar(binnedDt1(1,:), binnedDt2(1,:), binnedDt2(2,:), 'k');
        iidx = intersect( find(~isnan(dt1)), find(~isnan(dt2)) );
        [r,p] = corr(dt1(iidx), dt2(iidx));
        xlabel(['r = ' num2str(r) '; p = ' num2str(p)]);
        count = count+1;
        
    end
end
subplot(2,3,1); title('avr1-2');
subplot(2,3,2); title('avr2-3');
subplot(2,3,3); title('avr3-vvr1');
subplot(2,3,4); title('vvr1-2');
subplot(2,3,5); title('vvr2-3');
set(gcf, 'position', [123 123 777 444]);





%%
