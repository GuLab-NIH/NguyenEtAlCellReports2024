
% load data
% contact Lead Author for data


%% number of cells - AO vs VS
% panel 2F

figure; hold on

% AO1-VS1;
xx = cellRatioNcells{1}; k = size(xx, 1);
plot([ones(1,k); ones(1,k)*2], xx', 'Color', [1 1 1]/2);
yy = xx(:,1); clr = 'm';
bar(1, mean(yy), clr, 'FaceAlpha', 0.5);
errorbar(1, mean(yy), std(yy)/sqrt(length(yy)), '.', 'Color', clr, ...
    'MarkerFaceColor', clr, 'LineWidth', 1.5);
yy = xx(:,2); clr = 'g';
bar(2, mean(yy), clr, 'FaceAlpha', 0.5);
errorbar(2, mean(yy), std(yy)/sqrt(length(yy)), '.', 'Color', clr, ...
    'MarkerFaceColor', clr, 'LineWidth', 1.5);

% AO2-VS2;
xx = cellRatioNcells{2}; k = size(xx, 1);
xx = xx([1:4:k 4:4:k],:); k = size(xx, 1);
plot([ones(1,k)*3; ones(1,k)*4], xx', 'Color', [1 1 1]/2);
yy = xx(:,1); clr = 'm';
bar(3, mean(yy), clr, 'FaceAlpha', 0.5);
errorbar(3, mean(yy), std(yy)/sqrt(length(yy)), '.', 'Color', clr, ...
    'MarkerFaceColor', clr, 'LineWidth', 1.5);
yy = xx(:,2); clr = 'g';
bar(4, mean(yy), clr, 'FaceAlpha', 0.5);
errorbar(4, mean(yy), std(yy)/sqrt(length(yy)), '.', 'Color', clr, ...
    'MarkerFaceColor', clr, 'LineWidth', 1.5);
    
xlim([0.5 4.5]); title('Number of Active Cells'); ylabel('Number');
set(gca, 'xtick', 1:4, 'xticklabel', {'AO1', 'VS1', 'AO2', 'VS2'});
set(gcf, 'Position', [200 200 288 288]);





%% ANATOMY - CELL OVERLAP PERCENTAGE
% panel 2G

% prcOverlap = cell(4, 4);
    % AOn - AOn / AOo / VSn / ~~~
    % AOo - ~~~ / AOo / ~~~ / VSo
    % VSn - ~~~ / ~~~ / VSn / VSo
    % VSo - ~~~ / ~~~ / ~~~ / VSo
    
load('cellOverlap.mat');

indices = [6 5 1 14 9 16 15 11];
rbgValues = round(linspace(0, 255, 8))/255;
figure; hold on
for ii = 1:numel(indices)
    switch indices(ii)
        case {1, 5, 6}
            clr = 'm';
        case {11, 15, 16}
            clr = 'g';
        case {9, 14}
            clr = 'k';
    end
    
    yy = finalOverlap{indices(ii)}*100;
    bar(ii, mean(yy), clr, 'FaceAlpha', 0.2);
    errorbar(ii, mean(yy), std(yy)/sqrt(length(yy)), 'Color', clr, ...
        'MarkerSize', 20, 'MarkerFaceColor', clr, 'LineWidth', 2);
end
xlim([0.5 8.5]); ylim([0 50]); ylabel('% cells overlapped');
set(gca, 'xtick', []);
set(gcf, 'Position', [200 200 666 400]);


% test significant between each pair of columns
for ii = 1:(numel(indices)-1)
for jj = (ii+1):numel(indices)
    aa = indices(ii); bb = indices(jj);
    [~, p] = ttest2(finalOverlap{aa}, finalOverlap{bb}, 'vartype', 'unequal');
    disp([num2str(aa) '-' num2str(bb) '  --  ' num2str(p)]);
end
end

% xx = cat(2, finalOverlap{[6 5 1]});
xx = cat(2, finalOverlap{[16 15 11]});
[~, p] = ttest2(xx, finalOverlap{14}); disp(p);
[~, p] = ttest2(xx, finalOverlap{9}); disp(p);




%% ANATOMY - CELL OVERLAP PERCENTAGE + random
% panel S2F

% prcOverlap = cell(4, 4);
    % AOn - AOn / AOo / VSn / ~~~
    % AOo - ~~~ / AOo / ~~~ / VSo
    % VSn - ~~~ / ~~~ / VSn / VSo
    % VSo - ~~~ / ~~~ / ~~~ / VSo
    
% load('cellOverlap.mat');

indices = [6 5 1 14 9 16 15 11];
figure; hold on
for ii = 1:numel(indices)
    switch indices(ii)
        case {1, 5, 6}
            clr = 'm';
        case {11, 15, 16}
            clr = 'g';
        case {9, 14}
            clr = 'k';
    end
    
    yy = randomOverlap2{indices(ii)}*100;
%     bar(ii, mean(yy), clr, 'FaceAlpha', 0.2);
    errorbar(ii+.2, mean(yy), std(yy)/sqrt(length(yy)), 'Color', clr, ...
        'MarkerSize', 20, 'MarkerFaceColor', clr, 'LineWidth', 2);
    
    yy = finalOverlap{indices(ii)}*100;
%     bar(ii, mean(yy), clr, 'FaceAlpha', 0.2);
    errorbar(ii-.2, mean(yy), std(yy)/sqrt(length(yy)), 'Color', clr, ...
        'MarkerSize', 20, 'MarkerFaceColor', clr, 'LineWidth', 2);
end
xlim([0.5 8.5]); ylim([0 50]); ylabel('% cells overlapped');
set(gca, 'xtick', []);
set(gcf, 'Position', [200 200 666 400]);





%% RBR consistencies for all VRs

% load('allRbr.mat');     % e6.mat in local drive

for VR = 1:4
    switch VR
        case 1  % AVR1
            ttl  = 'avr1';
            clr = 'm';
        case 2  % AVR2
            ttl  = 'avr2';
            clr = 'm';
        case 3  % VVR1
            ttl  = 'vvr1';
            clr = 'g';
        case 4  % VVR2
            ttl  = 'vvr2';
            clr = 'g';
    end
    
    figure; hold on
    xx = allRbr{1,VR};
    yy = allRbr{2,VR};
    [~,p] = ttest2(xx, yy);
    
    xx = xx( xx < nanmean(xx)+ 3*nanstd(xx) & xx > nanmean(xx)- 3*nanstd(xx) );
    violin2(xx, 'x', 1, 'facecolor', clr, 'facealpha', 0.3);
    
    yy = yy( yy < nanmean(yy)+ 3*nanstd(yy) & yy > nanmean(yy)- 3*nanstd(yy) );
    violin2(yy, 'x', 2, 'facecolor', 'k', 'facealpha', 0.3);
    
    title([ttl ' - p: ' num2str(p)]);
    xlim([0 2]+0.5);
    set(gcf, 'position', [111*VR 111*VR 234 222]);
end





%% PLOTTING examples of Correct v.v Error runs

load('allFolders.mat');


% folders = folders_audio;
% templt = zeros(1,42); templt([1:10 15:24 29:38]) = 1;

folders = folders_audio_NEW_learned;
templt = zeros(1,66); templt([5:14 21:30 39:48 53:62]) = 1;

% folders = folders_visual;
% templt = zeros(1,42); templt([1:10 15:24 29:38]) = 1;

% folders = folders_visual_NEW_learned;
% templt = zeros(1,66); templt([5:14 21:30 39:48 53:62]) = 1;


for ff = 1:numel(folders)

    % correct runs index
    load([folders{ff} '\corrIncorr_20230701\isCorrectNoFirstLast.mat']);
    
    try
        % fields & dfof
        load([folders{ff} '\corrIncorr_20230701\PValueClassifier_KY2_6_sigC\allCells.mat']);
        fields_C = allCells.inFieldBins;
        dfof_C = allCells.dfofaveragesmooth; clear allCells;
        load([folders{ff} '\corrIncorr_20230701\PValueClassifier_KY2_6_sigIC\allCells.mat']);
        fields_E = allCells.inFieldBins;
        dfof_E = allCells.dfofaveragesmooth; clear allCells;

        % dfofM
        load([folders{ff} '\corrIncorr_20230701\dfofMInterpM.mat']);
        rbr_C = dfofMInterpM; clear dfofMInterpM;
        load([folders{ff} '\corrIncorr_20230701\dfofMInterpM_IC.mat']);
        rbr_E = dfofMInterpM_IC; clear dfofMInterpM_IC;
    catch er
        continue
    end
    
    
    for nn = 1:numel(rbr_C)
        if (size(rbr_C{nn},1) == 1) || (size(rbr_E{nn},1) == 1)
            continue
        end
        
        tempRBR = rbr_C{nn};
        corr_C = calculateRBR_singleCell(tempRBR); corr_C = corr_C.meantoOthers;
        if corr_C < 0.2
            continue
        end
        
        figure;
        subplot(2,2,1); hold on
        imagesc(tempRBR); set(gca, 'xtick', [1 length(templt)], 'xticklabel', {'start', 'end'});
        xlim([0 length(templt)]+0.5); ylim([0 size(tempRBR,1)]+0.5); title('Correct');
        caxis([0 max(tempRBR, [], 'all')]);

        subplot(2,2,2); hold on
        tempRBR = rbr_E{nn};
        corr_E = calculateRBR_singleCell(tempRBR); corr_E = corr_E.meantoOthers;
        imagesc(tempRBR); set(gca, 'xtick', [1 length(templt)], 'xticklabel', {'start', 'end'});
        xlim([0 length(templt)]+0.5); ylim([0 size(tempRBR,1)]+0.5); title('Error');
        caxis([0 max(tempRBR, [], 'all')]);

        subplot(2,2,3); hold on
        df = dfof_C(:, nn);
        plot(df, 'k', 'LineWidth', 2);
        plot(templt*1.1*max(df)); clear df;
        yl = get(gca, 'yLim');
        bins = fields_C{nn};
        if ~isempty(bins)
            bar(bins, ones(1, length(bins))*(0.05*yl(2)), 'FaceColor', 'r', 'EdgeColor', 'r', 'BarWidth', 1);
        end
        xlim([0 length(templt)] + 0.5); ylim([0 yl(2)]);
        title(corr_C);

        subplot(2,2,4); hold on
        df = dfof_E(:, nn);
        plot(df, 'k', 'LineWidth', 2);
        plot(templt*1.1*max(df)); clear df;
        yl = get(gca, 'yLim');
        bins = fields_E{nn};
        if ~isempty(bins)
            bar(bins, ones(1, length(bins))*(0.05*yl(2)), 'FaceColor', 'r', 'EdgeColor', 'r', 'BarWidth', 1);
        end
        xlim([0 length(templt)] + 0.5); ylim([0 yl(2)]);
        title(corr_E);

        set(gcf, 'Position', [200 200 555 333]);
        pause;
        close;

    end
end




%% Average speed for all environments
% S3J

load('allFolders.mat');

allSpeeds = cell(1,4);     % VR

for VR = 1:4
    clear folders;
    switch VR
        case 1  % AVR1
            folders = folders_audio;
        case 2  % AVR2
            folders = folders_audio_NEW_learned;
        case 3  % VVR1
            folders = folders_visual;
        case 4  % VVR2
            folders = folders_visual_NEW_learned;
    end
    
    for ff = 1:numel(folders)
        disp(['VR = ' num2str(VR) ' - progress: ' num2str(ff) '/' num2str(numel(folders))]);
        
        try
            load([folders{ff} '\corrIncorr_20230701\abfFakeC.mat']);
        catch er
        end
        
        ss = diff([0; abfFakeC.y]) ./ diff([0 abfFakeC.t])';
        ss(ss < 0.12) = [];
        
        allSpeeds{VR} = [allSpeeds{VR}; nanmean(ss)];
    end
end


figure; hold on
count = 1;
for VR = [1 3 2 4]
    switch VR
        case {1,2}  % AVRs
            clr = 'm';
        case {3,4}  % VVRs
            clr = 'g';
    end
    
    xx = allSpeeds{VR};
%     [mm, se] = getMeanAndSE(xx, 'all');
%     bar(count, mm, clr, 'FaceAlpha', 0.5);
%     errorbar(count, mm, se, clr);
    xx = xx(xx < nanmean(xx) + 3*nanstd(xx) & xx > nanmean(xx) - 3*nanstd(xx));
    violin2(xx, 'x', count, 'facecolor', clr, 'facealpha', 0.3);
    count = count + 1;
end
[~,p1] = ttest2( allSpeeds{1}, allSpeeds{3} );
[~,p2] = ttest2( allSpeeds{2}, allSpeeds{4} );
title(['p1 = ' num2str(p1) '; p2 = ' num2str(p2)]);
xlim([0 4]+0.5);
% set(gca, 'xtick', 1:4, 'xticklabel', {'avr1', 'avr2', 'vvr1', 'vvr2'});
set(gca, 'xtick', 1:4, 'xticklabel', {'avr1', 'vvr1', 'avr2', 'vvr2'});
set(gcf, 'position', [222 222 234 222]);



%%

