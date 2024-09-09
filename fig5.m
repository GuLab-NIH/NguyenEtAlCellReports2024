
% load data
% contact Lead Author for data




%% run-consistency in AUDIO vs VISUAL
% panel 5N

figure; hold on
useR = 3;

% ao
iidx = [34 48 68 116]; xidx = 1;
xx = cat(1, cueConsistencies{useR,:}); xx = xx(:, iidx);
bar(xidx, nanmean(xx, 'all'), 'm', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(xx, 'all'), nanstd(xx, [], 'all')/sqrt(size(xx,1)), 'Color', 'm', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'm', 'LineWidth', 2, 'LineStyle', 'none');

% vs
iidx = [14 86 100 132]; xidx = 2;
yy = cat(1, cueConsistencies{useR,:}); yy = yy(:, iidx);
bar(xidx, nanmean(yy, 'all'), 'g', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(yy, 'all'), nanstd(yy, [], 'all')/sqrt(size(yy,1)), 'Color', 'g', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'g', 'LineWidth', 2, 'LineStyle', 'none');
    
xlim([0.5 2.5]); title('Runs Consistency');
set(gca, 'xtick', 1:2, 'xticklabel', {'AO', 'VS'});
set(gcf, 'Position', [123 123 222 222]);

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);



%% Mean-DFOF in AUDIO vs VISUAL
% panel 5L

figure; hold on

% ao
iidx = [31:36 45:50 65:70 113:118]; xidx = 1;
xx = cat(1, allDfofs{:}); xx = xx(:, iidx);
bar(xidx, nanmean(xx, 'all'), 'm', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(xx, 'all'), nanstd(xx, [], 'all')/sqrt(size(xx,1)), 'Color', 'm', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'm', 'LineWidth', 2, 'LineStyle', 'none');

% vs
iidx = [11:16 83:88 97:102 129:134]; xidx = 2;
yy = cat(1, allDfofs{:}); yy = yy(:, iidx);
bar(xidx, nanmean(yy, 'all'), 'g', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(yy, 'all'), nanstd(yy, [], 'all')/sqrt(size(yy,1)), 'Color', 'g', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'g', 'LineWidth', 2, 'LineStyle', 'none');
    
xlim([0.5 2.5]); title('Mean \DeltaF/F'); ylabel('\DeltaF/F');
set(gca, 'xtick', 1:2, 'xticklabel', {'AO', 'VS'});
set(gcf, 'Position', [123 123 234 222]);

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);



%% fields in AUDIO vs VISUAL
% panel 5M


fields = cellfun(@(x) mean(x,1), allFields, 'uni', 0);

figure; hold on

% ao
iidx = [31:36 45:50 65:70 113:118]; xidx = 1;
xx = cellfun(@(x) mean(x(iidx)), fields, 'uni', 0);
xx = cat(1, xx{:});
% xx = fields(:, iidx);
bar(xidx, nanmean(xx, 'all'), 'm', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(xx, 'all'), nanstd(xx, [], 'all')/sqrt(size(xx,1)), 'Color', 'm', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'm', 'LineWidth', 2, 'LineStyle', 'none');

% vs
iidx = [11:16 83:88 97:102 129:134]; xidx = 2;
yy = cellfun(@(x) mean(x(iidx)), fields, 'uni', 0);
yy = cat(1, yy{:});
% yy = fields(:, iidx);
bar(xidx, nanmean(yy, 'all'), 'g', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(yy, 'all'), nanstd(yy, [], 'all')/sqrt(size(yy,1)), 'Color', 'g', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'g', 'LineWidth', 2, 'LineStyle', 'none');
    
xlim([0.5 2.5]); title('Mean % cells'); ylabel('% cells');
set(gca, 'xtick', 1:2, 'xticklabel', {'AO', 'VS'});
set(gcf, 'Position', [123 123 234 222]);

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);



%%
