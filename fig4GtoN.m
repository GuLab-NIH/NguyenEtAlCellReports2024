
% load data
% contact Lead Author for data




%% CUE CORRELATION relationship
% calculate dfof correlation with cue template
% 4G-N

load('fig6_avr1avr2_10bin.mat'); ttl = 'avr1-avr2';
% load('fig6_avr1vvr1_10bin.mat'); ttl = 'avr1-vvr1';
% load('fig6_avr2vvr2_10bin.mat'); ttl = 'avr2-vvr2';
% load('fig6_vvr1vvr2_10bin.mat'); ttl = 'vvr1-vvr2';

dt1 = cat(2, cueDfofs{:,1});
dt2 = cat(2, cueDfofs{:,2});

if size(dt1,1) == 42
    temp1 = zeros(42,1);
    temp1([1:10 15:24 29:38]) = 1;
elseif size(dt1, 1) == 66
    temp1 = zeros(66,1);
    temp1([5:14 21:30 39:48 53:62]) = 1;
end

if size(dt2,1) == 42
    temp2 = zeros(42,1);
    temp2([1:10 15:24 29:38]) = 1;
elseif size(dt2, 1) == 66
    temp2 = zeros(66,1);
    temp2([5:14 21:30 39:48 53:62]) = 1;
end


dt1 = num2cell(dt1, 1);
dt1 = cellfun(@(x) corr(x, temp1), dt1, 'uni', 0);

dt2 = num2cell(dt2, 1);
dt2 = cellfun(@(x) corr(x, temp2), dt2, 'uni', 0);

dt1 = cat(1, dt1{:});
dt2 = cat(1, dt2{:});


% remove cells with no fields
f1 = cat(2, cueFields{:,1});
f2 = cat(2, cueFields{:,2});
i1 = find( sum(f1,1) > 0);
i2 = find( sum(f2,1) > 0);
iidx = intersect(i1, i2);
dt1 = dt1(iidx); dt2 = dt2(iidx);

edges = prctile(dt1, 0:10:100);
% edges = linspace(prctile(dt1,0), prctile(dt1,100), 10);
[~, ~, binIdx] = histcounts( dt1, edges );
binnedDt1 = nan(2, length(edges)-1);
binnedDt2 = nan(2, length(edges)-1);
for ii = 1:(length(edges)-1)
    bb = find(binIdx == ii);
    binnedDt1(1, ii) = nanmean( dt1(bb) );
    binnedDt1(2, ii) = nanstd( dt1(bb) ) / sqrt(length(bb));
    
    binnedDt2(1, ii) = nanmean( dt2(bb) );
    binnedDt2(2, ii) = nanstd( dt2(bb) ) / sqrt(length(bb));
end
figure; hold on
errorbar(binnedDt1(1,:), binnedDt2(1,:), binnedDt2(2,:), 'k');
iidx = intersect( find(~isnan(dt1)), find(~isnan(dt2)) );
[r,p] = corr(dt1(iidx), dt2(iidx));
title(ttl);
xlabel(['r = ' num2str(r) '; p = ' num2str(p)]);
set(gcf, 'position', [222 222 234 222]);




%% cell example
% 4G-H

% load('fig6_avr1avr2_10bin.mat'); ttls{1} = 'avr2'; ttls{2} = 'avr1';
% load('fig6_avr1vvr1_10bin.mat'); ttls{1} = 'avr1'; ttls{2} = 'vvr1';
% load('fig6_avr2vvr2_10bin.mat'); ttls{1} = 'avr2'; ttls{2} = 'vvr2';
% load('fig6_vvr1vvr2_10bin.mat'); ttls{1} = 'vvr2'; ttls{2} = 'vvr1';


mice = 1:6;

dfs_ref = cat(2, cueDfofs{mice, 1});
fields_ref = cat(2, cueFields{mice, 1});
runs_ref = cat(2, cueSFruns{mice, 1});
rbr_ref = cat(1, cueRbrs{mice, 1});

dfs_target = cat(2, cueDfofs{mice, 2});
fields_target = cat(2, cueFields{mice, 2});
runs_target = cat(2, cueSFruns{mice, 2});
rbr_target = cat(1, cueRbrs{mice, 2});

if size(dfs_ref,1) == 42
    template_ref = zeros(42,1);
    template_ref([1:10 15:24 29:38]) = 1;
elseif size(dfs_ref, 1) == 66
    template_ref = zeros(66,1);
    template_ref([5:14 21:30 39:48 53:62]) = 1;
end

if size(dfs_target,1) == 42
    template_target = zeros(42,1);
    template_target([1:10 15:24 29:38]) = 1;
elseif size(dfs_target, 1) == 66
    template_target = zeros(66,1);
    template_target([5:14 21:30 39:48 53:62]) = 1;
end



dt1 = num2cell(dfs_ref, 1);
dt1 = cellfun(@(x) naninterp(x), dt1, 'uni', 0);
dt1 = cellfun(@(x) corr(x, template_ref), dt1, 'uni', 0);
dt2 = num2cell(dfs_target, 1);
dt2 = cellfun(@(x) naninterp(x), dt2, 'uni', 0);
dt2 = cellfun(@(x) corr(x, template_target), dt2, 'uni', 0);
dt1 = cat(1, dt1{:});
dt2 = cat(1, dt2{:});

cueScores = cellfun(@(x) x', cueScores, 'uni', 0);
sc1 = cat(2, cueScores{:, 1})'; sc1 = sc1(:, 2);
sc2 = cat(1, cueScores{:, 2}); sc2 = sc2(:, 1);

xx = find(sc1 < 0.1); 
yy = find(sc2 > 0.2);
indices = intersect(xx, yy);    % top corner cue cell

if 1
for ii = 1:numel(indices)
    ccID = indices(ii); disp(ii);
    
    figure;
    subplot(221); hold on
    aa = runs_ref(:, 1:(ccID-1)); aa = ~isnan(aa); aa = sum(aa, 'all')+1;
    bb = runs_ref(:, 1:(ccID)); bb = ~isnan(bb); bb = sum(bb, 'all');
    imagesc(rbr_ref(aa:bb,:)); set(gca, 'xtick', [1 length(template_ref)], 'xticklabel', {'start', 'end'});
    xlim([0 length(template_ref)]+0.5); ylim([0 (bb-aa+1)]+0.5); title(ttls{1}); ylabel('runs');
    
    subplot(222); hold on
    aa = runs_target(:, 1:(ccID-1)); aa = ~isnan(aa); aa = sum(aa, 'all')+1;
    bb = runs_target(:, 1:(ccID)); bb = ~isnan(bb); bb = sum(bb, 'all');
    imagesc(rbr_target(aa:bb,:)); set(gca, 'xtick', [1 length(template_target)], 'xticklabel', {'start', 'end'});
    xlim([0 length(template_target)]+0.5); ylim([0 (bb-aa+1)]+0.5); title(ttls{2});
    
    subplot(223); hold on
    df = dfs_ref(:, ccID);
    plot(df, 'k', 'LineWidth', 2);
    plot(template_ref*1.1*max(df)); clear df;
    yl = get(gca, 'yLim');
    bins = find(fields_ref(:, ccID));
    if ~isempty(bins); bar(bins, ones(1, length(bins))*(0.05*yl(2)), 'FaceColor', 'r', 'EdgeColor', 'r', 'BarWidth', 1); end
    xlim([0 length(template_ref)] + 0.5); set(gca, 'xtick', [1 length(template_ref)], 'xticklabel', {'start', 'end'}); ylabel('\DeltaF/F');
    title(dt1(ccID));
    
    subplot(224); hold on
    df = dfs_target(:, ccID);
    plot(df, 'k', 'LineWidth', 2);
    plot(template_target*1.1*max(df)); clear df;
    yl = get(gca, 'yLim');
    bins = find(fields_target(:, ccID));
    if ~isempty(bins); bar(bins, ones(1, length(bins))*(0.05*yl(2)), 'FaceColor', 'r', 'EdgeColor', 'r', 'BarWidth', 1); end
    xlim([0 length(template_target)] + 0.5); set(gca, 'xtick', [1 length(template_target)], 'xticklabel', {'start', 'end'});
    title(dt2(ccID));
    
    close();
end
end




%%
