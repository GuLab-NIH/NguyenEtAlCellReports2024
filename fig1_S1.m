
% load data

% contact Lead Author for data



%% PERCENTAGE COMMON-UNIQUE CELLS
% panel 1I-J


prcOfCells = nan(numel(folders), 4);
% 1 - unique audioVR cells / total audioVR cells
% 2 - common audioVR cells / total audioVR cells
% 3 - common nosound cells / total nosound cells
% 4 - unique nosound cells / total nosound cells

nCells = nan(numel(folders), 3);

for ff = 1:length(folders)
    
    dd = dir([folders{ff} '\1\1_2\cellRegistered*']);
    load([folders{ff} '\1\1_2\' dd(1).name]);
    tempCommonCells = cell_registered_struct.cell_to_index_map;
    nComCells = 0;
    nAOcells = 0;
    nNScells = 0;
    
    for ii = 1:size(tempCommonCells, 1)
        x = tempCommonCells(ii,:);
        if prod(x)
            nComCells = nComCells + 1;
        elseif x(1)
            nAOcells = nAOcells + 1;
        elseif x(2)
            nNScells = nNScells + 1;
        end
    end
    prcOfCells(ff,:) = [nAOcells/(nComCells+nAOcells)...
        nComCells/(nComCells+nAOcells)...
        nComCells/(nComCells+nNScells)...
        nNScells/(nComCells+nNScells)];
    
    nCells(ff,:) = [nAOcells nComCells nNScells];
end



% PLOTTING PERCENTAGE OVERLAPPED
figure;
TEMP1 = nCells(:,1) + nCells(:, 2);
TEMP2 = nCells(:,3) + nCells(:, 2);
x = TEMP1./TEMP2;
hold on
bar(1, mean(x), 'm', 'BarWidth', 0.8, 'FaceAlpha', 0.3, 'EdgeColor', 'm');
errorbar(1, mean(x), std(x)/sqrt(length(x)), 'm', 'LineWidth', 2)

bar(2, 1, 'k', 'BarWidth', 0.8, 'FaceAlpha', 0.3, 'EdgeColor', 'k');

[h, p] = ttest(x); disp(p);
if h; plot([1 2], [1.9 1.9], 'k', 'LineWidth', 1); text(1.5, 1.95, '*', 'FontSize', 20, 'HorizontalAlignment', 'center'); end

hold off
xlim([0.5 2.5])
set(gca, 'ytick', 0:0.5:2, 'yticklabel', (0:0.5:2), 'xtick', 1:2, 'xticklabel', {'audioVR', 'nosound'});
ylabel('Percentage');
title('Ratio of number of active cells');
set(gcf, 'Position', [200 200 222 222]);


figure; hold on
violin2(TEMP1, 'x', 1, 'facecolor', 'm', 'facealpha', 0.3);
violin2(TEMP2, 'x', 2, 'facecolor', 'k', 'facealpha', 0.3);

xlim([0.5 2.5]);
[~,p] = ttest(TEMP1, TEMP2);
ylabel('number of cells'); title(p);
set(gcf, 'Position', [111 111 222 222]);



% PLOTTING CELL PERCENTAGE
xx = (nCells ./ sum(nCells, 2))*100;
figure; hold on
for ii = 1:3
    switch ii; case 1; clr = 'm'; case 2; clr = 'k'; case 3; clr = 'k'; end
    violin2(xx(:,ii), 'x', ii, 'facecolor', clr, 'facealpha', 0.3);
end
xlim([0.5 3.5]); ylabel('Percentage');
set(gca, 'xtick', 1:3); set(gca, 'xticklabel', {'u-AO', 'common', 'u-NS'});
set(gcf, 'Position', [200 200 258 222]);

[~,p] = ttest(xx(:,1), xx(:,2)); disp(p);
[~,p] = ttest(xx(:,1), xx(:,3)); disp(p);
[~,p] = ttest(xx(:,2), xx(:,3)); disp(p);









%% PLOT CELL SEQUENCES
% 1M



xx1 = cat(1, fieldsDistribution{:,1,2});
yy1 = cat(2, dfofDistribution{:,1,2});
yy1 = normalize(yy1, 1, 'range')';
iidx1 = find(sum(xx1,2));

xx2 = cat(1, fieldsDistribution{:,2,2});
yy2 = cat(2, dfofDistribution{:,2,2});
yy2 = normalize(yy2, 1, 'range')';
iidx2 = find(sum(xx2,2));

iidx = intersect(iidx1, iidx2);
% iidx = 1:1557;
avr = yy1(iidx,:) .* xx1(iidx,:);
na = yy2(iidx,:) .* xx2(iidx,:);
figure; plotSequenceByMax2(avr', na');





% AO
xx = cat(1, fieldsDistribution{:,1,2});
yy = cat(2, dfofDistribution{:,1,2});
yy = normalize(yy, 1, 'range')';
% iidx = find(sum(xx,2));
% 
yy = yy(iidx,:) .* xx(iidx,:);
figure; plotSequenceByMax(yy); colormap((parula)); title('Fields-Amplitude');
set(gcf, 'Position', [222 222 345 333]);

% NS
xx = cat(1, fieldsDistribution{:,2,2});
yy = cat(2, dfofDistribution{:,2,2});
yy = normalize(yy, 1, 'range')';
% iidx = find(sum(xx,2));
% 
yy = yy(iidx,:) .* xx(iidx,:);
figure; plotSequenceByMax(yy); colormap((parula)); title('Fields-Amplitude');
set(gcf, 'Position', [333 333 345 333]);
% colorbar



%% Vector population analysis
% 1O

% calculate vector population analysis
corrMat = nan(42,42);
for iAVR = 1:42
for iNA = 1:42
    corrMat(iAVR, iNA) = corr( avr(:,iAVR), na(:,iNA) );
end
end

figure; hold on
imagesc(corrMat);
xlim([0 42]+.5); ylim([0 42]+.5);
axis square;
colorbar; caxis([-0.1 0.4]);
xlabel('NA'); ylabel('AVR1');
set(gcf, 'position', [222 222 444 345]);

shCorrMat = nan( 42, 42, 1000 );
for ss = 1:1000
    shEnv_a = avr(randperm(size(avr,1)), :);
    for aa = 1:42
    for bb = 1:42
        shCorrMat(aa, bb, ss) = corr( shEnv_a(:,aa), na(:,bb) );
    end
    end
end
temp = prctile(shCorrMat, 99, 3);

figure; hold on
imagesc(corrMat > temp);
xlim([0 42]+.5); ylim([0 42]+.5);
axis square;
colormap(gray);
xlabel('NA'); ylabel('AVR1');
set(gcf, 'position', [222 222 444 345]);


figure;
subplot(1,2,1); hold on
imagesc(corrMat);
xlim([0 42]+.5); ylim([0 42]+.5);
axis square;
colorbar; caxis([-0.1 0.4]);
xlabel('NA'); ylabel('AVR1');

subplot(1,2,2); hold on
imagesc(corrMat > temp);
xlim([0 42]+.5); ylim([0 42]+.5);
axis square;
colormap(gray);
xlabel('NA'); ylabel('AVR1');

set(gcf, 'position', [222 222 777 345]);



%% Vector population analysis


% calculate vector population analysis for avr1-avr1
load('AoldTwice.mat');

flds = cat(1, fieldsInfo{:});
dfofs = cat(1, allDfof{:});

xx3 = reshape(flds(:,1), 42, [])';
yy3 = reshape(dfofs(:,1), 42, [])';
yy3 = normalize(yy3, 1, 'range');
iidx3 = find(sum(xx3,2));

xx4 = reshape(flds(:,2), 42, [])';
yy4 = reshape(dfofs(:,2), 42, [])';
yy4 = normalize(yy4, 1, 'range');
iidx4 = find(sum(xx4,2));

iidx = intersect(iidx3, iidx4);
avr_a = yy3(iidx,:) .* xx3(iidx,:);
avr_b = yy4(iidx,:) .* xx4(iidx,:);

avr_a = naninterp(avr_a')';
avr_b = naninterp(avr_b')';

corrMat = nan(42,42);
for aa = 1:42
for bb = 1:42
    corrMat(aa, bb) = corr( avr_a(:,aa), avr_b(:,bb) );
end
end

figure; hold on
imagesc(corrMat);
xlim([0 42]+.5); ylim([0 42]+.5);
axis square;
colorbar; caxis([-0.1 0.4]);
xlabel('avr1'); ylabel('avr1');
set(gcf, 'position', [222 222 444 345]);


shCorrMat = nan( 42, 42, 1000 );
for ss = 1:1000
    shEnv_a = avr_a(randperm(size(avr_a,1)), :);
    for aa = 1:42
    for bb = 1:42
        shCorrMat(aa, bb, ss) = corr( shEnv_a(:,aa), avr_b(:,bb) );
    end
    end
end
temp = prctile(shCorrMat, 99, 3);

figure; hold on
imagesc(corrMat > temp);
xlim([0 42]+.5); ylim([0 42]+.5);
axis square;
colormap(gray);
xlabel('avr1'); ylabel('avr1');
set(gcf, 'position', [222 222 444 345]);




% calculate vector population analysis for NA-NA
load('NAtwice.mat');

flds = cat(1, fieldsInfo{:});
dfofs = cat(1, allDfof{:});

xx3 = reshape(flds(:,1), 42, [])';
yy3 = reshape(dfofs(:,1), 42, [])';
yy3 = normalize(yy3, 1, 'range');
iidx3 = find(sum(xx3,2));

xx4 = reshape(flds(:,2), 42, [])';
yy4 = reshape(dfofs(:,2), 42, [])';
yy4 = normalize(yy4, 1, 'range');
iidx4 = find(sum(xx4,2));

iidx = intersect(iidx3, iidx4);
avr_a = yy3(iidx,:) .* xx3(iidx,:);
avr_b = yy4(iidx,:) .* xx4(iidx,:);

corrMat = nan(42,42);
for aa = 1:42
for bb = 1:42
    corrMat(aa, bb) = corr( avr_a(:,aa), avr_b(:,bb) );
end
end

figure; hold on
imagesc(corrMat);
xlim([0 42]+.5); ylim([0 42]+.5);
axis square;
colorbar; caxis([-0.1 0.4]);
xlabel('NA'); ylabel('NA');
set(gcf, 'position', [222 222 444 345]);


shCorrMat = nan( 42, 42, 1000 );
for ss = 1:1000
    shEnv_a = avr_a(randperm(size(avr_a,1)), :);
    for aa = 1:42
    for bb = 1:42
        shCorrMat(aa, bb, ss) = corr( shEnv_a(:,aa), avr_b(:,bb) );
    end
    end
end
temp = prctile(shCorrMat, 99, 3);

figure; hold on
imagesc(corrMat > temp);
xlim([0 42]+.5); ylim([0 42]+.5);
axis square;
colormap(gray);
xlabel('NA'); ylabel('NA');
set(gcf, 'position', [222 222 444 345]);





%% FIELDS DISTRIBUTION ALONG TRACK
% 1N

allFields = cellfun(@(x) mean(x,1), fieldsDistribution, 'uni', 0);
xidx = 1:42;

% UNIQUE cells
figure; hold on
xx = cat(1, allFields{:,1,1}); xx = xx *100; % *100 for percentage
[tempmean, ~, seUpper, seLower] = getMeanAndSE(xx, 1);
plot(xidx, tempmean(xidx), 'm', 'LineWidth', 2);
patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'm', 'FaceAlpha', 0.2, 'LineStyle', 'none');
yy = cat(1, allFields{:,2,1}); yy = yy *100; % *100 for percentage
[tempmean, ~, seUpper, seLower] = getMeanAndSE(yy, 1);
plot(xidx, tempmean(xidx), 'k', 'LineWidth', 2);
patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'k', 'FaceAlpha', 0.2, 'LineStyle', 'none');
yl = get(gca, 'ylim'); plotCueAndRewardZone(yl);

slideWindow = 0;
for ii = 1:(42-slideWindow)
    aa = xx(:,ii:(ii+slideWindow)); aa = reshape(aa, 1, []);
    bb = yy(:,ii:(ii+slideWindow)); bb = reshape(bb, 1, []);
    [hh,~] = ttest2( aa, bb );
    if hh
        text(ii+slideWindow/2, 40, '*', 'Color', 'k');
    end
end

xlabel('Distance along track (cm)'); ylabel('% fields');
set(gca, 'xtick', [0 42] + 0.5, 'xticklabel', [0 210]); xlim([0 42]+0.5);
set(gcf, 'Position', [200 200 444 234]);



% COMMON cells
figure; hold on
xx = cat(1, allFields{:,1,2}); xx = xx *100; % *100 for percentage
[tempmean, ~, seUpper, seLower] = getMeanAndSE(xx, 1);
plot(xidx, tempmean(xidx), 'm', 'LineWidth', 2);
patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'm', 'FaceAlpha', 0.2, 'LineStyle', 'none');
yy = cat(1, allFields{:,2,2}); yy = yy *100; % *100 for percentage
[tempmean, ~, seUpper, seLower] = getMeanAndSE(yy, 1);
plot(xidx, tempmean(xidx), 'k', 'LineWidth', 2);
patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'k', 'FaceAlpha', 0.2, 'LineStyle', 'none');
yl = get(gca, 'ylim'); plotCueAndRewardZone(yl);

slideWindow = 0;
for ii = 1:(42-slideWindow)
    aa = xx(:,ii:(ii+slideWindow)); aa = reshape(aa, 1, []);
    bb = yy(:,ii:(ii+slideWindow)); bb = reshape(bb, 1, []);
    [hh,~] = ttest2( aa, bb );
    if hh
        text(ii+slideWindow/2, 40, '*', 'Color', 'k');
    end
end

xlabel('Distance along track (cm)'); ylabel('% fields');
set(gca, 'xtick', [0 42] + 0.5, 'xticklabel', [0 210]); xlim([0 42]+0.5);
set(gcf, 'Position', [200 200 444 234]);



%% average spatial fields bin per cell
% 1L

allFields  = cellfun(@(x) nanmean(x,2), fieldsDistribution, 'uni', 0);


% UNIQUE cells
figure; hold on
xx = cat(1, allFields{:,1,1}); xx = xx *100; % *100 for percentage
violin2(xx, 'x', 1, 'facecolor', 'm', 'facealpha', 0.3);

yy = cat(1, allFields{:,2,1}); yy = yy *100; % *100 for percentage
violin2(yy, 'x', 2, 'facecolor', 'k', 'facealpha', 0.3);

[~,p] = ttest2(xx, yy);
xlim([0.5 2.5]); title(['UNIQUE: p = ' num2str(p)]);
set(gca, 'xtick', 1:2, 'xticklabel', {'avr1', 'ns'});
set(gcf, 'Position', [200 200 234 222]);



% COMMON cells
figure; hold on
xx = cat(1, allFields{:,1,2}); xx = xx *100; % *100 for percentage
violin2(xx, 'x', 1, 'facecolor', 'm', 'facealpha', 0.3);

yy = cat(1, allFields{:,2,2}); yy = yy *100; % *100 for percentage
violin2(yy, 'x', 2, 'facecolor', 'k', 'facealpha', 0.3);

[~,p] = ttest2(xx, yy);
xlim([0.5 2.5]); title(['COMMON: p = ' num2str(p)]);
set(gca, 'xtick', 1:2, 'xticklabel', {'avr1', 'ns'});
set(gcf, 'Position', [200 200 234 222]);



%% RUN-BY-RUN CONSISTENCY
% 1P

useData = rbrConsistencies; 
useR = 3;       % corrToOthers
for ii = 2:3
    switch ii
        case 1  % ALL CELLS
            indices = 1:2;
        case 2  % COMMON
            indices = 2;
        case 3  % UNIQUE
            indices = 1;
    end
    
    xx = cat(1, useData{:, 1, indices});    % AudioVR
    xx = xx(:, useR);
    yy = cat(1, useData{:, 2, indices});    % NoSound
    yy = yy(:, useR);
    
	figure; hold on
    xx = xx( xx < nanmean(xx) + 3*nanstd(xx) & xx > nanmean(xx) - 3*nanstd(xx) );
    violin2(xx, 'x', 1, 'facecolor', 'm', 'facealpha', 0.3);
    
    yy = yy( yy < nanmean(yy) + 3*nanstd(yy) & yy > nanmean(yy) - 3*nanstd(yy) );
    violin2(yy, 'x', 2, 'facecolor', 'k', 'facealpha', 0.3);
    
    try
        [~, p] = ttest(xx, yy); disp(p);
    catch er
        [~, p] = ttest2(xx, yy); disp(p);
    end
    
    xlim([0.5 2.5]); set(gcf, 'Position', [50 50 500 250]);
    set(gca, 'xtick', 1:2, 'xticklabel', {'AO', 'NS'});
    
    ylabel('rbr correlation');
    set(gcf, 'Position', [111*ii 111*ii 234 222]);
end




%% PLOT CUE CONSISTENCIES
% 1Q


% just separate them
cellType = 3;   % all = 1; common = 2; unique = 3;
corrType = 3;   % toMean/toNext/toOthers

figure; hold on;
xidx = 1:3; zones = [5 20 32];
for dayType = 1:2
    switch dayType
        case 1
            clr = 'm';
        case 2
            clr = 'k';
    end
    xx = cuesConsistencies{dayType, corrType, cellType};
    [tempmean, se, ~, ~] = getMeanAndSE(xx, 1);
    tempmean = tempmean(zones);
    se = se(zones);
    bar(xidx + 4*(dayType-1), tempmean, 'FaceColor', clr, 'FaceAlpha', 0.2);
    errorbar(xidx + 4*(dayType-1), tempmean, se, clr, 'LineWidth', 2,'LineStyle', 'none');
    
end
yl = get(gca, 'ylim'); xlim([0.5 7.5]);
plot([2.5 2.5], yl, 'y-', 'LineWidth', 1.5);
plot([5.5 5.5]+1, yl, 'y-', 'LineWidth', 1.5);
set(gca, 'xtick', [1:3 5:7], 'xticklabel', [1:3 1:3]);
set(gcf, 'Position', [200 200 234 222]);





%% RUN & CUE CONSISTENCIES for Correct vs Error runs
% for avr1 only

% field distribution
fields = cellfun(@(x) mean(x,1), allFields, 'uni', 0);
xidx = 1:42;
for cellType = 1:2 % common/unique
    switch cellType
        case 1      % common
            ttl = 'common';
        case 2      % unique
            ttl = 'unique';
    end
    
    figure; hold on
    
    xx = cat(1, fields{:, cellType, 1}); xx = xx *100;
    [mm, ~, seUpper, seLower] = getMeanAndSE(xx, 1);
    plot(xidx, mm(xidx), 'm', 'LineWidth', 2);
    patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'm', 'FaceAlpha', 0.2, 'LineStyle', 'none');
    
    yy = cat(1, fields{:, cellType, 2}); yy = yy *100;
    [mm, ~, seUpper, seLower] = getMeanAndSE(yy, 1);
    plot(xidx, mm(xidx), 'k', 'LineWidth', 2);
    patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'k', 'FaceAlpha', 0.2, 'LineStyle', 'none');
    
    yl = get(gca, 'ylim'); plotCueAndRewardZone(yl);
    
    slideWindow = 0;
    for ii = 1:(42-slideWindow)
        aa = xx(:,ii:(ii+slideWindow)); aa = reshape(aa, 1, []);
        bb = yy(:,ii:(ii+slideWindow)); bb = reshape(bb, 1, []);
        [hh,~] = ttest2( aa, bb );
        if hh
            text(ii+slideWindow/2, 35, '*', 'Color', 'k');
        end
    end
    xlabel('Distance along track (cm)'); ylabel('% fields'); title(ttl);
    set(gca, 'xtick', [0 42] + 0.5, 'xticklabel', [0 210]); xlim([0 42]+0.5);
    set(gcf, 'Position', [111*cellType 111*cellType 444 234]);
end



% run consistency
for cellType = 1:2
    switch cellType
        case 1      % common
            ttl = 'common';
        case 2      % unique
            ttl = 'unique';
    end
    figure; hold on
    
    xx = cat(1, allRbrCorrs{:, cellType, 1});
    xx = xx( xx < nanmean(xx) + 3*nanstd(xx) & xx > nanmean(xx) - 3*nanstd(xx) );
    violin2(xx, 'x', 1, 'facecolor', 'r', 'facealpha', 0.3);
    
    yy = cat(1, allRbrCorrs{:, cellType, 2});
    yy = yy( yy < nanmean(yy) + 3*nanstd(yy) & yy > nanmean(yy) - 3*nanstd(yy) );
    violin2(yy, 'x', 2, 'facecolor', 'b', 'facealpha', 0.3);
    
    [~,p] = ttest2(xx, yy);
    title([ttl ' - p: ' num2str(p)]);
    xlim([0 2]+0.5);
    set(gcf, 'position', [111*cellType 111*cellType 234 222]);
end


rr = 2.5;
% cue consistency
for cellType = 1:2
    switch cellType
        case 1      % common
            ttl = 'common';
        case 2      % unique
            ttl = 'unique';
    end
    figure; hold on
    
    xx = cat(1, allCueCorrs{:, cellType, 1});
    
    for kk = 1:3
        zz = xx(:,kk);
        zz = zz( zz < nanmean(zz) + 3*nanstd(zz) & zz > nanmean(zz) - 3*nanstd(zz) );
        violin2(zz, 'x', kk, 'facecolor', 'r', 'facealpha', 0.3);
    end
    
    yy = cat(1, allCueCorrs{:, cellType, 2});
    
    for kk = 1:3
        zz = yy(:,kk);
        zz = zz( zz < nanmean(zz) + 3*nanstd(zz) & zz > nanmean(zz) - 3*nanstd(zz) );
        violin2(zz, 'x', kk+4, 'facecolor', 'b', 'facealpha', 0.3);
    end
    
    [~,p1] = ttest2(xx(:,1), yy(:,1));
    [~,p2] = ttest2(xx(:,2), yy(:,2));
    [~,p3] = ttest2(xx(:,3), yy(:,3));
    
    [~,p1] = ttest2(yy(:,1), yy(:,2));
    [~,p2] = ttest2(yy(:,1), yy(:,3));
    [~,p3] = ttest2(yy(:,2), yy(:,3));
    
    
    title([ttl ' - p(s): ' num2str(p1) ' - ' num2str(p2) ' - ' num2str(p3)]);
    yl = get(gca, 'ylim'); xlim([0.5 7.5]);
    plot([rr rr], yl, 'y-', 'LineWidth', 1.5);
    plot([rr rr]+4, yl, 'y-', 'LineWidth', 1.5);
    set(gca, 'xtick', [1:3 5:7], 'xticklabel', [1:3 1:3]);
    set(gcf, 'position', [111*cellType 111*cellType 444 222]);
end





%% PLOTTING CELL EXAMPLE
% panel 1K


count = 1;
for ff = randperm(length(folders), length(folders))
    
    % load common cells
    dd = dir([folders{ff} '\1\1_2\cellRegistered*']);
    load([folders{ff} '\1\1_2\' dd(1).name]);
    commonCells = cell_registered_struct.cell_to_index_map;
    
    % load Audio infos
    load([folders{ff} '\audioVR\pcaica\speedAndAccelerationThresholded\dfofMInterpM.mat']);
    dfofM_AO = dfofMInterpM; clear dfofMInterpM;
    load([folders{ff} '\audioVR\pcaica\speedAndAccelerationThresholded\dfofM.mat']);
    dfofM_AO2 = dfofM; clear dfofM;
    load([folders{ff} '\audioVR\pcaica\speedAndAccelerationThresholded\dfofaveragesmooth.mat']);
    dfof_AO = dfofaveragesmooth; clear dfofaveragesmooth;
    load([folders{ff} '\audioVR\pcaica\newFields20\allCells.mat']);
    fields_AO = allCells.inFieldBins; clear allCells;
    
    % load Nosound infos
    load([folders{ff} '\nosound\pcaica\speedAndAccelerationThresholded\dfofMInterpM.mat']);
    dfofM_NS = dfofMInterpM; clear dfofMInterpM;
    load([folders{ff} '\nosound\pcaica\speedAndAccelerationThresholded\dfofM.mat']);
    dfofM_NS2 = dfofM; clear dfofM;
    load([folders{ff} '\nosound\pcaica\speedAndAccelerationThresholded\dfofaveragesmooth.mat']);
    dfof_NS = dfofaveragesmooth; clear dfofaveragesmooth;
    load([folders{ff} '\nosound\pcaica\newFields20\allCells.mat']);
    fields_NS = allCells.inFieldBins; clear allCells;
    
    pairs = find(prod(commonCells,2));
    
    for nn = 1:numel(pairs)
        x1 = commonCells(pairs(nn), 1); x2 = commonCells(pairs(nn), 2);
        
        corr1 = calculateRBR_singleCell(dfofM_AO{x1}); corr1 = corr1.meantoOthers;
        corr2 = calculateRBR_singleCell(dfofM_NS{x2}); corr2 = corr2.meantoOthers;
        
        if corr1 < 0.05
            continue
        end
        
        figure;
        subplot(2,2,1); hold on
%         imagesc(dfofM_AO2{x1}); set(gca, 'xtick', [1 42], 'xticklabel', {'start', 'end'});
        imagesc(naninterp_for_dfofM_noRemove1stLast(dfofM_AO2{x1})); set(gca, 'xtick', [1 42], 'xticklabel', {'start', 'end'});
        caxis([0 max(dfofM_AO{x1}, [], 'all')]);
        xlim([0 42]+0.5); ylim([0 size(dfofM_AO2{x1},1)]+0.5); title('AudioVR');
        
        subplot(2,2,2); hold on
%         imagesc(dfofM_NS2{x2}); set(gca, 'xtick', [1 42], 'xticklabel', {'start', 'end'});
        imagesc(naninterp_for_dfofM_noRemove1stLast(dfofM_NS2{x2})); set(gca, 'xtick', [1 42], 'xticklabel', {'start', 'end'});
        caxis([0 max(dfofM_NS{x2}, [], 'all')]);
        xlim([0 42]+0.5); ylim([0 size(dfofM_NS2{x2},1)]+0.5); title('NoSound');
        
        subplot(2,2,3); hold on
        df = dfof_AO(:, x1);
        plot(df, 'k', 'LineWidth', 2);
        title(corr1);
        plot(template*1.1*max(df)); clear df;
        
        yl = get(gca, 'yLim');
        bins = fields_AO{x1};
        if ~isempty(bins)
            bar(bins, ones(1, length(bins))*(0.05*yl(2)), 'FaceColor', 'r', 'EdgeColor', 'r', 'BarWidth', 1);
        end
        xlim([0 42] + 0.5);
        
        subplot(2,2,4); hold on
        df = dfof_NS(:, x2);
        plot(df, 'k', 'LineWidth', 2);
        title(corr2);
        plot(template*1.1*max(df)); clear df;
        
        yl = get(gca, 'yLim');
        bins = fields_NS{x2};
        if ~isempty(bins)
            bar(bins, ones(1, length(bins))*(0.05*yl(2)), 'FaceColor', 'r', 'EdgeColor', 'r', 'BarWidth', 1);
        end
        xlim([0 42] + 0.5);
        
        set(gcf, 'Position', [200 200 555 333]);
        
        savefig(['Y:\labMembers\DN\_PROJECT_\finalData\figures\revision\v5\NAexamples\ex' num2str(count) '.fig']); count = count + 1;
%         pause;
        close;
    end
end



%%


function plotTemplate(template)
    hold on
    plot(template, 'k', 'LineWidth', 1.2);
    yl = get(gca, 'ylim');
    patch([120 140 140 120], [yl(1) yl(1) yl(2) yl(2)], 'y', 'FaceAlpha', 0.2);
    hold off
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    xlim([0 210] + 0.5);
    ylim([66 85]);
end

function plotCueAndRewardZone(yl)
    patch([24 28 28 24]+0.5, [yl(1) yl(1) yl(2) yl(2)], 'y', 'FaceAlpha', 0.2, 'LineStyle', 'none');

    plot([10 10]+0.5, [yl(1) yl(2)], 'k:');
    plot([14 14]+0.5, [yl(1) yl(2)], 'k:');
    plot([24 24]+0.5, [yl(1) yl(2)], 'k:');
    plot([28 28]+0.5, [yl(1) yl(2)], 'k:');
    plot([38 38]+0.5, [yl(1) yl(2)], 'k:');
end


