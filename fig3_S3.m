
% load data
% contact Lead Author for data



%% CELL SEQUENCE
% 3B-C


%%%% COMMON CELLS BETWEEN AO1 - VS1

df = cat(2, allDfofs_1{1,1,:});
fields = cat(2, allFields_1{1,1,:});
seq1 = df.*fields;

df = cat(2, allDfofs_1{2,1,:});
fields = cat(2, allFields_1{2,1,:});
seq2 = df.*fields;

iidx = find(sum(seq1,1)>0 & sum(seq2,1)>0);
figure; plotSequenceByMax(seq1(:,iidx)); colormap((parula)); title('AO1');
set(gcf, 'Position', [111 111 345 333]);
figure; plotSequenceByMax(seq2(:,iidx)); colormap((parula)); title('VS1');
set(gcf, 'Position', [222 222 345 333]);



%%%% COMMON CELLS BETWEEN AO2 - VS2

df = cat(2, allDfofs_2{1,1,:});
fields = cat(2, allFields_2{1,1,:});
seq1 = df.*fields;

df = cat(2, allDfofs_2{2,1,:});
fields = cat(2, allFields_2{2,1,:});
seq2 = df.*fields;

iidx = find(sum(seq1,1)>0 & sum(seq2,1)>0);
figure; plotSequenceByMax(seq1(:,iidx)); colormap((parula)); title('AO2');
set(gcf, 'Position', [111 111 345 333]);
figure; plotSequenceByMax(seq2(:,iidx)); colormap((parula)); title('VS2');
set(gcf, 'Position', [222 222 345 333]);




%% FIELDS avr1-vvr1
% 3D

allFields = cellfun(@(x) mean(x,2), allFields_1, 'uni', 0);

xidx = 2:42;
template = [1:10 15:24 29:39]; rew = 25:28;
figure; hold on
xx = cat(2, allFields{1,1,:}); xx = xx*100; xx = xx';
[tempmean, ~, seUpper, seLower] = getMeanAndSE(xx, 1);
plot(xidx, tempmean(xidx), 'm', 'LineWidth', 2);
patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'm', 'FaceAlpha', 0.35, 'LineStyle', 'none');

yy = cat(2, allFields{2,1,:}); yy = yy*100; yy = yy';
[tempmean, ~, seUpper, seLower] = getMeanAndSE(yy, 1);
plot(xidx, tempmean(xidx), 'g', 'LineWidth', 2);
patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'g', 'FaceAlpha', 0.35, 'LineStyle', 'none');

% cue and reward
yl = get(gca, 'ylim');
bar(template-0.5, ones(1, length(template))*yl(2), 'BaseValue', yl(1), 'EdgeColor', 'none', 'BarWidth', 1, 'FaceColor', 'k', 'FaceAlpha', 0.1);
bar(rew-0.5, ones(1, length(rew))*yl(2), 'BaseValue', yl(1), 'EdgeColor', 'none', 'BarWidth', 1, 'FaceColor', 'y', 'FaceAlpha', 0.3);

slideWindow = 0;
for ii = 1:(42-slideWindow)
    aa = xx(:,ii:(ii+slideWindow)); aa = reshape(aa, 1, []);
    bb = yy(:,ii:(ii+slideWindow)); bb = reshape(bb, 1, []);
    [hh,~] = ttest2( aa, bb, 'Tail', 'left' );
    if hh
        text(ii+slideWindow/2, 30, '*', 'Color', 'm');
    end
    [hh,~] = ttest2( aa, bb, 'Tail', 'left' );
    if hh
        text(ii+slideWindow/2, 30, '*', 'Color', 'g');
    end
end

xlabel('Distance along track (cm)'); ylabel('% fields');
xlim([0 42]+0.5); ylim(yl);
set(gca, 'xtick', [0 42] + 0.5, 'xticklabel', [0 210]); xlim([0 42]+0.5);
set(gcf, 'Position', [200 200 444 234]);


%% FIELDS avr1-vvr1 - in/out ratio

inArea = [3:8 17:22 31:36];
outArea = setdiff(1:42, inArea);
allFields = cellfun(@(x) mean(x,2), allFields_1, 'uni', 0);
allFields = cellfun(@(x) mean(x(inArea)) / mean(x(outArea)), allFields, 'uni', 0);

figure; hold on
xx = cat(2, allFields{1,1,:});
yy = cat(2, allFields{2,1,:});

bar(1, nanmean(xx, 'all'), 'm', 'FaceAlpha', 0.3);
errorbar(1, nanmean(xx, 'all'), nanstd(xx, [], 'all')/sqrt(size(xx,1)), 'Color', 'm', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'm', 'LineWidth', 2, 'LineStyle', 'none');
    
bar(2, nanmean(yy, 'all'), 'g', 'FaceAlpha', 0.3);
errorbar(2, nanmean(yy, 'all'), nanstd(yy, [], 'all')/sqrt(size(yy,1)), 'Color', 'g', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'g', 'LineWidth', 2, 'LineStyle', 'none');

xlim([0.5 2.5]); %title('In/Out-field ratio');
set(gca, 'xtick', 1:2, 'xticklabel', {'AO', 'VS'});
set(gcf, 'Position', [123 123 234 222]);

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);



%% FIELD avr1-vvr1 - CUE-AROUND-REWARD / CUE-AWAY-FROM-REWARD ratio

zones{1} = 3:8; zones{2} = 17:22; zones{3} = 31:36;
allFields = cellfun(@(x) mean(x,2), allFields_1, 'uni', 0);


%%%%%% BOTH CUES AROUND REWARD

figure; hold on
% audio 
xx = cat(2, allFields{1,1,:})';
xx1 = mean(xx(:, zones{1}), 2); xx2 = mean(xx(:, cat(2, zones{2:3})), 2);
[tempmean, se, ~, ~] = getMeanAndSE(xx1, 1);
bar(1, tempmean, 'FaceColor', 'm', 'FaceAlpha', 0.2);
errorbar(1, tempmean, se, 'm', 'LineWidth', 2,'LineStyle', 'none');

[tempmean, se, ~, ~] = getMeanAndSE(xx2, 1);
bar(2, tempmean, 'FaceColor', 'm', 'FaceAlpha', 0.2);
errorbar(2, tempmean, se, 'm', 'LineWidth', 2,'LineStyle', 'none');

% visual
yy = cat(2, allFields{2,1,:})';
yy1 = mean(yy(:, zones{1}), 2); yy2 = mean(yy(:, cat(2, zones{2:3})), 2);
[tempmean, se, ~, ~] = getMeanAndSE(yy1, 1);
bar(3, tempmean, 'FaceColor', 'g', 'FaceAlpha', 0.2);
errorbar(3, tempmean, se, 'g', 'LineWidth', 2,'LineStyle', 'none');

[tempmean, se, ~, ~] = getMeanAndSE(yy2, 1);
bar(4, tempmean, 'FaceColor', 'g', 'FaceAlpha', 0.2);
errorbar(4, tempmean, se, 'g', 'LineWidth', 2,'LineStyle', 'none');

[~,p] = ttest(xx1, xx2); disp(p);
[~,p] = ttest(yy1, yy2); disp(p);

yl = get(gca, 'ylim'); xlim([0.5 4.5]); ylim(yl);
set(gca, 'xtick', 1:4, 'xticklabel', {'AVR - far', 'AVR - near', 'VVR - far', 'VVR - near'}, 'xticklabelrotation', 45);
set(gcf, 'Position', [111 111 222 222]);


% FOLD
xx3 = xx2 ./ xx1; xx3(isinf(xx3)) = [];
yy3 = yy2 ./ yy1; yy3(isinf(yy3)) = [];

figure; hold on
[tempmean, se, ~, ~] = getMeanAndSE(xx3, 1);
bar(1, tempmean, 'FaceColor', 'm', 'FaceAlpha', 0.2);
errorbar(1, tempmean, se, 'm', 'LineWidth', 2,'LineStyle', 'none');

[tempmean, se, ~, ~] = getMeanAndSE(yy3, 1);
bar(2, tempmean, 'FaceColor', 'g', 'FaceAlpha', 0.2);
errorbar(2, tempmean, se, 'g', 'LineWidth', 2,'LineStyle', 'none');

[~,p] = ttest2(xx3, yy3); disp(p);
[~,p] = ttest(xx3, 1); disp(p);
[~,p] = ttest(yy3, 1); disp(p);

yl = get(gca, 'ylim'); xlim([0.5 2.5]); ylim(yl);
set(gca, 'xtick', 1:5);
set(gcf, 'Position', [444 111 222 222]);



%% FIELDS avr2-vvr2
% 3D

allFields = cellfun(@(x) mean(x,2), allFields_2, 'uni', 0);

xidx = 2:65;
template = [5:14 21:30 39:48 53:62]; rew = 49:52;
figure; hold on
xx = cat(2, allFields{1,1,:}); xx = xx*100; xx = xx';
[tempmean, ~, seUpper, seLower] = getMeanAndSE(xx, 1);
plot(xidx, tempmean(xidx), 'm', 'LineWidth', 2);
patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'm', 'FaceAlpha', 0.35, 'LineStyle', 'none');

yy = cat(2, allFields{2,1,:}); yy = yy*100; yy = yy';
[tempmean, ~, seUpper, seLower] = getMeanAndSE(yy, 1);
plot(xidx, tempmean(xidx), 'g', 'LineWidth', 2);
patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'g', 'FaceAlpha', 0.35, 'LineStyle', 'none');

% cue and reward
yl = get(gca, 'ylim');
bar(template-0.5, ones(1, length(template))*yl(2), 'BaseValue', yl(1), 'EdgeColor', 'none', 'BarWidth', 1, 'FaceColor', 'k', 'FaceAlpha', 0.1);
bar(rew-0.5, ones(1, length(rew))*yl(2), 'BaseValue', yl(1), 'EdgeColor', 'none', 'BarWidth', 1, 'FaceColor', 'y', 'FaceAlpha', 0.3);

slideWindow = 0;
for ii = 1:(66-slideWindow)
    aa = xx(:,ii:(ii+slideWindow)); aa = reshape(aa, 1, []);
    bb = yy(:,ii:(ii+slideWindow)); bb = reshape(bb, 1, []);
    [hh,~] = ttest2( aa, bb, 'Tail', 'left' );
    if hh
        text(ii+slideWindow/2, 30, '*', 'Color', 'm');
    end
    [hh,~] = ttest2( aa, bb, 'Tail', 'left' );
    if hh
        text(ii+slideWindow/2, 30, '*', 'Color', 'g');
    end
end

xlabel('Distance along track (cm)'); ylabel('% fields');
xlim([0 66]+0.5); ylim(yl);
set(gca, 'xtick', [0 66] + 0.5, 'xticklabel', [0 330]);
set(gcf, 'Position', [200 200 444 234]);



%% FIELDS avr2-vvr2 - in/out ratio

inArea = [7:12 23:28 41:46 55:60];
outArea = setdiff(1:66, inArea);
allFields = cellfun(@(x) mean(x,2), allFields_2, 'uni', 0);
allFields = cellfun(@(x) sum(x(inArea)) / sum(x(outArea)), allFields, 'uni', 0);

figure; hold on
xx = cat(2, allFields{1,1,:});
yy = cat(2, allFields{2,1,:});

bar(1, nanmean(xx, 'all'), 'm', 'FaceAlpha', 0.3);
errorbar(1, nanmean(xx, 'all'), nanstd(xx, [], 'all')/sqrt(size(xx,1)), 'Color', 'm', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'm', 'LineWidth', 2, 'LineStyle', 'none');
    
bar(2, nanmean(yy, 'all'), 'g', 'FaceAlpha', 0.3);
errorbar(2, nanmean(yy, 'all'), nanstd(yy, [], 'all')/sqrt(size(yy,1)), 'Color', 'g', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'g', 'LineWidth', 2, 'LineStyle', 'none');

xlim([0.5 2.5]); % title('In/Out-field ratio'); ylabel('% cells');
set(gca, 'xtick', 1:2, 'xticklabel', {'AO', 'VS'});
set(gcf, 'Position', [123 123 234 222]);

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);



%% FIELD avr2-vvr2 - CUE-AROUND-REWARD / CUE-AWAY-FROM-REWARD ratio

zones{1} = 7:12; zones{2} = 23:28; zones{3} = 41:46; zones{4} = 55:60;
allFields = cellfun(@(x) mean(x,2), allFields_2, 'uni', 0);


%%%%%% BOTH CUES AROUND REWARD

figure; hold on
% audio 
xx = cat(2, allFields{1,1,:})';
xx1 = mean(xx(:, cat(2, zones{1:2})), 2); xx2 = mean(xx(:, cat(2, zones{3:4})), 2);
[tempmean, se, ~, ~] = getMeanAndSE(xx1, 1);
bar(1, tempmean, 'FaceColor', 'm', 'FaceAlpha', 0.2);
errorbar(1, tempmean, se, 'm', 'LineWidth', 2,'LineStyle', 'none');

[tempmean, se, ~, ~] = getMeanAndSE(xx2, 1);
bar(2, tempmean, 'FaceColor', 'm', 'FaceAlpha', 0.2);
errorbar(2, tempmean, se, 'm', 'LineWidth', 2,'LineStyle', 'none');

% visual
yy = cat(2, allFields{2,1,:})';
yy1 = mean(yy(:, cat(2, zones{1:2})), 2); yy2 = mean(yy(:, cat(2, zones{3:4})), 2);
[tempmean, se, ~, ~] = getMeanAndSE(yy1, 1);
bar(3, tempmean, 'FaceColor', 'g', 'FaceAlpha', 0.2);
errorbar(3, tempmean, se, 'g', 'LineWidth', 2,'LineStyle', 'none');

[tempmean, se, ~, ~] = getMeanAndSE(yy2, 1);
bar(4, tempmean, 'FaceColor', 'g', 'FaceAlpha', 0.2);
errorbar(4, tempmean, se, 'g', 'LineWidth', 2,'LineStyle', 'none');

[~,p] = ttest(xx1, xx2); disp(p);
[~,p] = ttest(yy1, yy2); disp(p);

yl = get(gca, 'ylim'); xlim([0.5 4.5]); ylim(yl);
set(gca, 'xtick', 1:4, 'xticklabel', {'AVR - far', 'AVR - near', 'VVR - far', 'VVR - near'}, 'xticklabelrotation', 45);
set(gcf, 'Position', [111 111 222 222]);


% FOLD
xx3 = xx2 ./ xx1; xx3(isinf(xx3)) = [];
yy3 = yy2 ./ yy1; yy3(isinf(yy3)) = [];

figure; hold on
[tempmean, se, ~, ~] = getMeanAndSE(xx3, 1);
bar(1, tempmean, 'FaceColor', 'm', 'FaceAlpha', 0.2);
errorbar(1, tempmean, se, 'm', 'LineWidth', 2,'LineStyle', 'none');

[tempmean, se, ~, ~] = getMeanAndSE(yy3, 1);
bar(2, tempmean, 'FaceColor', 'g', 'FaceAlpha', 0.2);
errorbar(2, tempmean, se, 'g', 'LineWidth', 2,'LineStyle', 'none');

[~,p] = ttest2(xx3, yy3); disp(p);
[~,p] = ttest(xx3, 1); disp(p);
[~,p] = ttest(yy3, 1); disp(p);

yl = get(gca, 'ylim'); xlim([0.5 2.5]); ylim(yl);
set(gca, 'xtick', 1:5);
set(gcf, 'Position', [444 111 222 222]);







%% mean dfof
% 3A


figure; hold on
ao1 = nonAverageDfofValues_1{1,1};
bar(1, mean(ao1), 'm', 'FaceAlpha', 0.5);
errorbar(1, mean(ao1), std(ao1)/sqrt(length(ao1)), 'm', 'LineWidth', 1.5);

vs1 = nonAverageDfofValues_1{2,1};
bar(2, mean(vs1), 'g', 'FaceAlpha', 0.5);
errorbar(2, mean(vs1), std(vs1)/sqrt(length(vs1)), 'g', 'LineWidth', 1.5);

ao2 = nonAverageDfofValues_2{1,1};
bar(3, mean(ao2), 'm', 'FaceAlpha', 0.5);
errorbar(3, mean(ao2), std(ao2)/sqrt(length(ao2)), 'm', 'LineWidth', 1.5);

vs2 = nonAverageDfofValues_2{2,1};
bar(4, mean(vs2), 'g', 'FaceAlpha', 0.5);
errorbar(4, mean(vs2), std(vs2)/sqrt(length(vs2)), 'g', 'LineWidth', 1.5);

xlim([0 4]+0.5); ylabel('mean \DeltaF/F'); title('COMMON');
set(gca, 'xtick', 1:4, 'xticklabel', {'ao1', 'vs1', 'ao2', 'vs2'});
set(gcf, 'Position', [222 222 345 234]);

[~,p] = ttest2(ao1, vs1); disp(['ao1 - vs1: ' num2str(p)]);
[~,p] = ttest2(ao1, vs2); disp(['ao1 - ao2: ' num2str(p)]);
[~,p] = ttest2(ao1, vs2); disp(['ao1 - vs2: ' num2str(p)]);
[~,p] = ttest2(vs1, ao2); disp(['vs1 - ao2: ' num2str(p)]);
[~,p] = ttest2(vs1, vs2); disp(['vs1 - vs2: ' num2str(p)]);
[~,p] = ttest2(ao2, vs2); disp(['ao2 - vs2: ' num2str(p)]);




%% RBR consistency 
% 3E

%%%%%%%% COMMON %%%%%%%%
figure; hold on
useR = 3;
ao1 = rbrConsistency_1{1,1}; ao1 = ao1(:,useR);
bar(1, nanmean(ao1), 'm', 'FaceAlpha', 0.5);
errorbar(1, nanmean(ao1), nanstd(ao1)/sqrt(length(ao1)), 'm', 'LineWidth', 1.5);

vs1 = rbrConsistency_1{2,1}; vs1 = vs1(:,useR);
bar(2, nanmean(vs1), 'g', 'FaceAlpha', 0.5);
errorbar(2, nanmean(vs1), nanstd(vs1)/sqrt(length(vs1)), 'g', 'LineWidth', 1.5);

ao2 = rbrConsistency_2{1,1}; ao2 = ao2(:,useR);
bar(3, nanmean(ao2), 'm', 'FaceAlpha', 0.5);
errorbar(3, nanmean(ao2), nanstd(ao2)/sqrt(length(ao2)), 'm', 'LineWidth', 1.5);

vs2 = rbrConsistency_2{2,1}; vs2 = vs2(:,useR);
bar(4, nanmean(vs2), 'g', 'FaceAlpha', 0.5);
errorbar(4, nanmean(vs2), nanstd(vs2)/sqrt(length(vs2)), 'g', 'LineWidth', 1.5);

xlim([0 4]+0.5); ylabel('run consistency'); title('COMMON');
set(gca, 'xtick', 1:4, 'xticklabel', {'ao1', 'vs1', 'ao2', 'vs2'});
set(gcf, 'Position', [111 111 345 234]);

[~,p] = ttest2(ao1, vs1); disp(['ao1 - vs1: ' num2str(p)]);
[~,p] = ttest2(ao1, vs2); disp(['ao1 - ao2: ' num2str(p)]);
[~,p] = ttest2(ao1, vs2); disp(['ao1 - vs2: ' num2str(p)]);
[~,p] = ttest2(vs1, ao2); disp(['vs1 - ao2: ' num2str(p)]);
[~,p] = ttest2(vs1, vs2); disp(['vs1 - vs2: ' num2str(p)]);
[~,p] = ttest2(ao2, vs2); disp(['ao2 - vs2: ' num2str(p)]);




%% CUE CONSISTENCY avr1-vvr1
% 3H

figure; hold on
useR = 2; cueIdx = [6 20 34];

vs1 = cat(1, cueConsistency_1{2,1,useR}); vs1 = vs1(:, cueIdx);
bar(1:3, nanmean(vs1, 1), 'g', 'FaceAlpha', 0.5);
errorbar(1:3, nanmean(vs1,1), nanstd(vs1,1)/(sqrt(size(vs1,1))),...
    'g', 'LineStyle', 'none', 'LineWidth', 1.5);

ao1 = cat(1, cueConsistency_1{1,1,useR}); ao1 = ao1(:, cueIdx);
bar(1:3, nanmean(ao1, 1), 'm', 'FaceAlpha', 0.5);
errorbar(1:3, nanmean(ao1,1), nanstd(ao1,1)/(sqrt(size(ao1,1))),...
    'm', 'LineStyle', 'none', 'LineWidth', 1.5);

xlim([0 3]+0.5); ylabel('run consistency');
set(gca, 'xtick', 1:3);
set(gcf, 'Position', [111 111 345 234]);

[~,p] = ttest2(ao1(:,1), vs1(:,1)); disp(p);
[~,p] = ttest2(ao1(:,2), vs1(:,2)); disp(p);
[~,p] = ttest2(ao1(:,3), vs1(:,3)); disp(p);



%% CUE CONSISTENCY avr2-vvr2

figure; hold on
useR = 1; cueIdx = [10 26 44 58];

vs1 = cat(1, cueConsistency_2{2,1,useR}); vs1 = vs1(:, cueIdx);
bar(1:4, nanmean(vs1, 1), 'g', 'FaceAlpha', 0.5);
errorbar(1:4, nanmean(vs1,1), nanstd(vs1,1)/(sqrt(size(vs1,1))),...
    'g', 'LineStyle', 'none', 'LineWidth', 1.5);

ao1 = cat(1, cueConsistency_2{1,1,useR}); ao1 = ao1(:, cueIdx);
bar(1:4, nanmean(ao1, 1), 'm', 'FaceAlpha', 0.5);
errorbar(1:4, nanmean(ao1,1), nanstd(ao1,1)/(sqrt(size(ao1,1))),...
    'm', 'LineStyle', 'none', 'LineWidth', 1.5);

yl = get(gca, 'ylim'); plot([3.5 3.5], yl, 'y--', 'LineWidth', 2);
xlim([0 4]+0.5); ylabel('run consistency');
set(gca, 'xtick', 1:4);
set(gcf, 'Position', [111 111 345 234]);

[~,p] = ttest2(ao1(:,1), vs1(:,1)); disp(p);
[~,p] = ttest2(ao1(:,2), vs1(:,2)); disp(p);
[~,p] = ttest2(ao1(:,3), vs1(:,3)); disp(p);
[~,p] = ttest2(ao1(:,4), vs1(:,4)); disp(p);



%% CUE CONSISTENCY avr1-vvr1 - CUE-AROUND-REWARD / CUE-AWAY-FROM-REWARD diff
% 3I

zones = [6 20 34];
cuesCorr = cueConsistency_FOV_1;


%%%%%% BOTH CUES AROUND REWARD

figure; hold on
% audio 
xx = cat(1, cuesCorr{1,1,:});
xx1 = mean(xx(:, zones(1)), 2); xx2 = mean(xx(:, zones(2:3)), 2);
[tempmean, se, ~, ~] = getMeanAndSE(xx1, 1);
bar(1, tempmean, 'FaceColor', 'm', 'FaceAlpha', 0.2);
errorbar(1, tempmean, se, 'm', 'LineWidth', 2,'LineStyle', 'none');

[tempmean, se, ~, ~] = getMeanAndSE(xx2, 1);
bar(2, tempmean, 'FaceColor', 'm', 'FaceAlpha', 0.2);
errorbar(2, tempmean, se, 'm', 'LineWidth', 2,'LineStyle', 'none');

% visual
yy = cat(1, cuesCorr{2,1,:});
yy1 = mean(yy(:, zones(1)), 2); yy2 = mean(yy(:, zones(2:3)), 2);
[tempmean, se, ~, ~] = getMeanAndSE(yy1, 1);
bar(3, tempmean, 'FaceColor', 'g', 'FaceAlpha', 0.2);
errorbar(3, tempmean, se, 'g', 'LineWidth', 2,'LineStyle', 'none');

[tempmean, se, ~, ~] = getMeanAndSE(yy2, 1);
bar(4, tempmean, 'FaceColor', 'g', 'FaceAlpha', 0.2);
errorbar(4, tempmean, se, 'g', 'LineWidth', 2,'LineStyle', 'none');

[~,p] = ttest(xx1, xx2); disp(p);
[~,p] = ttest(yy1, yy2); disp(p);

yl = get(gca, 'ylim'); xlim([0.5 4.5]); ylim(yl);
set(gca, 'xtick', 1:4, 'xticklabel', {'AVR - far', 'AVR - near', 'VVR - far', 'VVR - near'}, 'xticklabelrotation', 45);
set(gcf, 'Position', [111 111 222 222]);


% DIFF
xx3 = xx2 - xx1;
yy3 = yy2 - yy1;

figure; hold on
[tempmean, se, ~, ~] = getMeanAndSE(xx3, 1);
bar(1, tempmean, 'FaceColor', 'm', 'FaceAlpha', 0.2);
errorbar(1, tempmean, se, 'm', 'LineWidth', 2,'LineStyle', 'none');

[tempmean, se, ~, ~] = getMeanAndSE(yy3, 1);
bar(2, tempmean, 'FaceColor', 'g', 'FaceAlpha', 0.2);
errorbar(2, tempmean, se, 'g', 'LineWidth', 2,'LineStyle', 'none');

[~,p] = ttest(xx3); disp(p);
[~,p] = ttest(yy3); disp(p);
[~,p] = ttest2(xx3, yy3); disp(p);

yl = get(gca, 'ylim'); xlim([0.5 2.5]); ylim(yl);
set(gca, 'xtick', 1:5);
set(gcf, 'Position', [444 111 222 222]);



%% CUE CONSISTENCY avr2-vvr2
% 3H

zones = [10 26 44 58]; xidx = 1:4;
cuesCorr = cueConsistency_FOV_2;

%%%%%%%% RAW
figure; hold on; 

% audio 
xx = cat(1, cuesCorr{1,1,:}); xx = xx(:, zones);
[tempmean, se, ~, ~] = getMeanAndSE(xx, 1);
bar(xidx, tempmean, 'FaceColor', 'm', 'FaceAlpha', 0.2);
errorbar(xidx, tempmean, se, 'm', 'LineWidth', 2,'LineStyle', 'none');

% visual
yy = cat(1, cuesCorr{2,1,:}); yy = yy(:, zones);
[tempmean, se, ~, ~] = getMeanAndSE(yy, 1);
bar(xidx + length(xidx) +1, tempmean, 'FaceColor', 'g', 'FaceAlpha', 0.2);
errorbar(xidx + length(xidx) +1, tempmean, se, 'g', 'LineWidth', 2,'LineStyle', 'none');

yl = get(gca, 'ylim'); xlim([0.5 9.5]); ylim(yl);
plot([3.5 3.5], yl, 'y-', 'LineWidth', 1.5);
plot([8.5 8.5], yl, 'y-', 'LineWidth', 1.5);
set(gca, 'xtick', 1:9);
set(gcf, 'Position', [111 111 345 222]);

[~, p] = ttest2( xx(:,1), yy(:,1) ); disp(p);
[~, p] = ttest2( xx(:,2), yy(:,2) ); disp(p);
[~, p] = ttest2( xx(:,3), yy(:,3) ); disp(p);
[~, p] = ttest2( xx(:,4), yy(:,4) ); disp(p);

[~, p] = ttest2( reshape(xx(:,1:2),[],1), reshape(xx(:,3:4),[],1) ); disp(p);
[~, p] = ttest2( reshape(yy(:,1:2),[],1), reshape(yy(:,3:4),[],1) ); disp(p);

[~, p] = ttest2( reshape(xx(:,1:2),[],1), reshape(xx(:,3),[],1) ); disp(p);
[~, p] = ttest2( reshape(yy(:,1:2),[],1), reshape(yy(:,3),[],1) ); disp(p);

[~, p] = ttest( reshape(xx(:,[1 1 2 2]),[],1), reshape(xx(:,[3 4 3 4]),[],1) ); disp(p);
[~, p] = ttest( reshape(yy(:,[1 1 2 2]),[],1), reshape(yy(:,[3 4 3 4]),[],1) ); disp(p);

[~, p] = ttest( reshape(xx(:,[1 2]),[],1), reshape(xx(:,[3 3]),[],1) ); disp(p);
[~, p] = ttest( reshape(yy(:,[1 2]),[],1), reshape(yy(:,[3 3]),[],1) ); disp(p);





%% CUE CONSISTENCY avr2-vvr2 - CUE-AROUND-REWARD / CUE-AWAY-FROM-REWARD diff
% 3I

zones = [10 26 44 58];
cuesCorr = cueConsistency_FOV_2;


%%%%%% BOTH CUES AROUND REWARD

figure; hold on
% audio 
xx = cat(1, cuesCorr{1,1,:});
xx1 = mean(xx(:, zones(1:2)), 2); xx2 = mean(xx(:, zones(3:4)), 2);
[tempmean, se, ~, ~] = getMeanAndSE(xx1, 1);
bar(1, tempmean, 'FaceColor', 'm', 'FaceAlpha', 0.2);
errorbar(1, tempmean, se, 'm', 'LineWidth', 2,'LineStyle', 'none');

[tempmean, se, ~, ~] = getMeanAndSE(xx2, 1);
bar(2, tempmean, 'FaceColor', 'm', 'FaceAlpha', 0.2);
errorbar(2, tempmean, se, 'm', 'LineWidth', 2,'LineStyle', 'none');

% visual
yy = cat(1, cuesCorr{2,1,:});
yy1 = mean(yy(:, zones(1:2)), 2); yy2 = mean(yy(:, zones(3:4)), 2);
[tempmean, se, ~, ~] = getMeanAndSE(yy1, 1);
bar(3, tempmean, 'FaceColor', 'g', 'FaceAlpha', 0.2);
errorbar(3, tempmean, se, 'g', 'LineWidth', 2,'LineStyle', 'none');

[tempmean, se, ~, ~] = getMeanAndSE(yy2, 1);
bar(4, tempmean, 'FaceColor', 'g', 'FaceAlpha', 0.2);
errorbar(4, tempmean, se, 'g', 'LineWidth', 2,'LineStyle', 'none');

[~,p] = ttest(xx1, xx2); disp(p);
[~,p] = ttest(yy1, yy2); disp(p);

yl = get(gca, 'ylim'); xlim([0.5 4.5]); ylim(yl);
set(gca, 'xtick', 1:4, 'xticklabel', {'AVR - far', 'AVR - near', 'VVR - far', 'VVR - near'}, 'xticklabelrotation', 45);
set(gcf, 'Position', [111 111 222 222]);


% DIFF
xx3 = xx2 - xx1;
yy3 = yy2 - yy1;

figure; hold on
[tempmean, se, ~, ~] = getMeanAndSE(xx3, 1);
bar(1, tempmean, 'FaceColor', 'm', 'FaceAlpha', 0.2);
errorbar(1, tempmean, se, 'm', 'LineWidth', 2,'LineStyle', 'none');

[tempmean, se, ~, ~] = getMeanAndSE(yy3, 1);
bar(2, tempmean, 'FaceColor', 'g', 'FaceAlpha', 0.2);
errorbar(2, tempmean, se, 'g', 'LineWidth', 2,'LineStyle', 'none');

[~,p] = ttest(xx3, 1); disp(p);
[~,p] = ttest(yy3, 1); disp(p);
[~,p] = ttest2(xx3, yy3); disp(p);

yl = get(gca, 'ylim'); xlim([0.5 2.5]); ylim(yl);
set(gca, 'xtick', 1:5);
set(gcf, 'Position', [444 111 222 222]);




%%
