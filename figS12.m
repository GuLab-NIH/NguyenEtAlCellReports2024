
%% ACCELERATION VS DFOF


%% load FOVs

load('allFolders');
folders = folders_audio_NEW_learned;


%% spatial bins
clearvars -except folders;

distanceFrame = 20; % 20cm before and after reward
binWidthPreferred = 2;

meanFiringRate_spatial = cell( 1, 4 );    % dfof - speed - acceleration - reward

for ff = 1:numel(folders)
    load([folders{ff} '\abfFakeNEW.mat']);
    dd = dir([folders{ff} '\dfof_sig*']);
    load([folders{ff} '\' dd(1).name]); clear dd;
    
    reward = find(abfFakeNEW.rewardsIdx);
    yy = abfFakeNEW.y;
    rewardssss = abfFakeNEW.rewardsIdx;
    speeds = ([0; diff(yy)]);       % can't normalize due to track ends
    acclrs = ([0; diff(speeds)]);   % should not be normalized

    for rIdx = 1:length(reward)
        lowLim = floor(yy(reward(rIdx)) - distanceFrame);
        upperLim = ceil(yy(reward(rIdx)) + distanceFrame);
        CBR = find(yy(1:reward(rIdx)) < lowLim, 1, 'last');
        CAR = find(yy(reward(rIdx):end) > upperLim, 1);

        if isempty(CBR) || isempty(CAR)
            continue
        else
            CAR = CAR + reward(rIdx);   % get correct index for CAR
            temp_y = yy(CBR:CAR);
            [counts, yIdx] = histc(temp_y, lowLim : binWidthPreferred : upperLim);
            df = dfof_sig(CBR:CAR, :);
            speed = speeds(CBR:CAR);
            acclr = acclrs(CBR:CAR);
            rews = rewardssss(CBR:CAR);
            temp_dfof = nan(1, length(counts)-1); temp_speed = nan(1, length(counts)-1); temp_acclr = nan(1, length(counts)-1); temp_rew = nan(1, length(counts)-1);

            for ii = 1:(length(counts)-1)
                temp_dfof(ii) = nanmean(df(yIdx == ii, :), 'all');
                temp_speed(ii) = nanmean(speed(yIdx == ii));
                temp_acclr(ii) = nanmean(acclr(yIdx == ii));
% 
%                 temp_speed(ii) = nansum(speed(yIdx == ii));
%                 temp_acclr(ii) = nansum(acclr(yIdx == ii));
                temp_rew(ii) = sum(rews(yIdx == ii));
            end
        end

        meanFiringRate_spatial{1} = [meanFiringRate_spatial{1}; temp_dfof];     % dfof
        meanFiringRate_spatial{2} = [meanFiringRate_spatial{2}; temp_speed];    % speed
        meanFiringRate_spatial{3} = [meanFiringRate_spatial{3}; temp_acclr];    % acceleration
        meanFiringRate_spatial{4} = [meanFiringRate_spatial{4}; temp_rew];      % reward
    end
end

target = 2;
figure;
hold on
YY = meanFiringRate_spatial{target};
[tempmean1, ~, ~, ~] = getMeanAndSE(YY, 1);

ZZ = meanFiringRate_spatial{1};
[tempmean2, ~, ~, ~] = getMeanAndSE(ZZ, 1);

xidx = 1:length(tempmean2);

colororder({'k','r'})
yyaxis left
plot(xidx, tempmean1, 'k', 'LineWidth', 2);
switch target; case 2; ylabel('Speed'); case 3; ylabel('Acceleration'); otherwise; warning('ERROR'); end

yyaxis right
plot(xidx, tempmean2, 'r', 'LineWidth', 1.5);
ylabel('\DeltaF/F');

hold off
set(gca, 'xtick', 0:5:distanceFrame, 'xticklabel', -distanceFrame:10:distanceFrame);
set(gcf, 'Position', [222 222 345 222]);
xlabel('Distance to reward delivery (cm);');
switch target; case 2; title('Speed'); case 3; title('Acceleration'); otherwise; warning('ERROR'); end



%% ALL TRACKS
clearvars -except folders;

binWidthPreferred = 2;
r1 = 240/binWidthPreferred; r2 = 260/binWidthPreferred;

meanFiringRate_alongTrack = cell( 1, 3 );    % dfof - speed - acceleration - reward

nBin = 330 / binWidthPreferred;
for ff = 1:numel(folders)
    try
        load([folders{ff} '\abfFake.mat']);
        dd = dir([folders{ff} '\dfof_sig*']);
        load([folders{ff} '\' dd(1).name]); clear dd;
    catch
        continue
    end
    
    yy = abfFake.y;
    speeds = ([0; diff(yy)]);       % can't normalize due to track ends
    acclrs = ([0; diff(speeds)]);   % should not be normalized
    speeds(speeds < -50) = nan;
    acclrs(acclrs < -50) = nan;
    acclrs(acclrs > 50) = nan;
    
    [~, idx] = histc(yy, 0 : binWidthPreferred : 330); 
    
    temp_dfof = nan(1, nBin); temp_speed = nan(1, nBin); temp_acclr = nan(1, nBin); temp_rew = nan(1, nBin);
    for ii = 1:nBin
        indices = find(idx == ii);
        temp_dfof(ii) = nanmean(dfof_sig(indices, :), 'all');
        temp_speed(ii) = nanmean(speeds(indices));
        temp_acclr(ii) = nanmean(acclrs(indices));
    end


    meanFiringRate_alongTrack{1} = [meanFiringRate_alongTrack{1}; temp_dfof];     % dfof
    meanFiringRate_alongTrack{2} = [meanFiringRate_alongTrack{2}; temp_speed];    % speed
    meanFiringRate_alongTrack{3} = [meanFiringRate_alongTrack{3}; temp_acclr];    % acceleration

end

% FIGURE
target = 2;
figure;
hold on
YY = meanFiringRate_alongTrack{target};
[tempmean1, ~, ~, ~] = getMeanAndSE(YY, 1);
ZZ = meanFiringRate_alongTrack{1};
[tempmean2, ~, ~, ~] = getMeanAndSE(ZZ, 1);

xidx = 1:length(tempmean2);
colororder({'k','r'})
yyaxis left
plot(xidx, tempmean1, 'k', 'LineWidth', 2);
switch target; case 2; ylabel('Speed'); case 3; ylabel('Acceleration'); otherwise; warning('ERROR'); end

yyaxis right
plot(xidx, tempmean2, 'r', 'LineWidth', 1.5);
ylabel('\DeltaF/F');

yl = get(gca, 'ylim');
plot([r1 r1], yl, 'y--', 'LineWidth', 2);
plot([r2 r2], yl, 'y--', 'LineWidth', 2);
hold off
set(gca, 'xtick', [0 nBin], 'xticklabel', [0 210]);
set(gcf, 'Position', [200 200 500 300]);
xlabel('Distance (cm);'); xlim([0 nBin]+0.5)
switch target; case 2; title('Speed'); case 3; title('Acceleration'); otherwise; warning('ERROR'); end

% MORE FIGURE
xx = meanFiringRate_alongTrack{1};
[tempmean1, ~, ~, ~] = getMeanAndSE(xx, 1);
tempmean1 = normalize(tempmean1, 'range');

xx = meanFiringRate_alongTrack{2};
[tempmean2, ~, ~, ~] = getMeanAndSE(xx, 1);
tempmean2 = normalize(tempmean2, 'range');

xx = meanFiringRate_alongTrack{3};
[tempmean3, ~, ~, ~] = getMeanAndSE(xx, 1);
tempmean3 = normalize(tempmean3, 'range');

xidx = 1:length(tempmean2);
figure;
hold on
plot(xidx, tempmean1, 'r', 'LineWidth', 1.5);
plot(xidx, tempmean2, 'b', 'LineWidth', 1.5);
plot(xidx, tempmean3, 'g', 'LineWidth', 1.5);
plot([r1 r1], [0 1], 'y--', 'LineWidth', 2);
plot([r2 r2], [0 1], 'y--', 'LineWidth', 2);
hold off
set(gca, 'xtick', [0 nBin], 'xticklabel', [0 210]);
set(gcf, 'Position', [200 200 500 300]);
xlabel('Distance (cm);'); xlim([0 nBin]+0.5); title('Normalized everything');
legend('\DeltaF/F', 'Speed', 'Acceleration', 'orientation', 'horizontal', 'location', 'south');



%% correlation along track

stepWise = 2;
areasWith = 20;
binWidth = areasWith/stepWise;
Radius = round(binWidth/2);
binWidth = binWidth-1;
r1 = 240/stepWise; r2 = 260/stepWise;
nBin = 330 / stepWise;
nn = nBin - binWidth;

correlationAlongTrack = cell( 2, 2 );    % corr speed-dfof / corr acceleration-dfof  - pvalue speed-dfof / pvalue acceleration-dfof

for ff = 1:numel(folders)
    load([folders{ff} '\abfFakeNEW.mat']);
    dd = dir([folders{ff} '\dfof_sig*']);
    load([folders{ff} '\' dd(1).name]); clear dd;
    
    reward = find(abfFakeNEW.rewardsIdx);
    yy = abfFakeNEW.y;
    rewardssss = abfFakeNEW.rewardsIdx;
    speeds = ([0; diff(yy)]);       % can't normalize due to track ends
    acclrs = ([0; diff(speeds)]);   % should not be normalized
    speeds(speeds < -50) = nan;
    acclrs(acclrs < -50) = nan;
    acclrs(acclrs > 50) = nan;
    
    [~, idx] = histc(yy, 0 : stepWise : 330); 
    
    temp_dfof = nan(1, nBin); temp_speed = nan(1, nBin); temp_acclr = nan(1, nBin); temp_rew = nan(1, nBin);
    for ii = 1:nBin
        indices = find(idx == ii);
        temp_dfof(ii + Radius) = nanmean(dfof_sig(indices, :), 'all');
        temp_speed(ii + Radius) = nanmean(speeds(indices));
        temp_acclr(ii + Radius) = nanmean(acclrs(indices));
        temp_rew(ii + Radius) = sum(rewardssss(indices));
    end
    
    corsM = nan( nBin, 1 );
    ppsM = nan( nBin, 1 );
    kk = find(temp_acclr, 1);
    for ii = kk:nn
        [corsM(ii), ppsM(ii)] = corr( temp_acclr(ii:(ii+binWidth))', temp_dfof(ii:(ii+binWidth))' );
    end
    correlationAlongTrack{2, 1} = [correlationAlongTrack{2, 1} corsM];
    correlationAlongTrack{2, 2} = [correlationAlongTrack{2, 2} ppsM];
    
    corsM = nan( nBin, 1 );
    ppsM = nan( nBin, 1 );
    kk = find(temp_speed, 1);
    for ii = kk:nn
        [corsM(ii), ppsM(ii)] = corr( temp_speed(ii:(ii+binWidth))', temp_dfof(ii:(ii+binWidth))' );
    end
    correlationAlongTrack{1, 1} = [correlationAlongTrack{1, 1} corsM];
    correlationAlongTrack{1, 2} = [correlationAlongTrack{1, 2} ppsM];

end

figure;
hold on
for target = 1:2
    switch target; case 1; clr = 'c'; case 2; clr = 'k'; end;
    xx1 = nanmean(correlationAlongTrack{target,1}, 2);
    xx2 = nanmean(correlationAlongTrack{target,2}, 2);
    % errorbar(1:length(xx1), xx1, sd1, 'm', 'LineWidth', 2);
    ppp(target) = plot(xx1, clr, 'LineWidth', 2);
%     plot(xx2, 'y');
    yl = get(gca, 'ylim');
    plot([r1 r1], yl, 'r');
    plot([r2 r2], yl, 'r');
end
hold off
% legend('Correlation', 'P-value', 'orientation', 'horizontal', 'location', 'south');
legend(ppp, 'SPEED', 'ACCELERATION', 'orientation', 'horizontal', 'location', 'south');
set(gca, 'xtick', [0 nBin]+0.5, 'xticklabel', [0 210]); xlim([0 nBin]+0.5);
set(gcf, 'Position', [200 200 500 300]);
title(['Width = ' num2str(areasWith) 'cm']);
xlabel('Distance'); ylabel('Correlation with \DeltaF/F');
% switch target; case 1; sgtitle(['SPEED - width = ' num2str(areasWith)]); case 2; sgtitle(['ACCELERATION - width = ' num2str(areasWith)]); otherwise; warning('ERROR'); end


%% SIGNIFICANT TRANSIENT ON-OFF vs ACCELERATION
% 1 - where sig transient happends
% 2 - distribution of accelration/ speed during sig transient

allSpeeds = [];
allAcceleration = [];
allDfof = [];
correlationWithinSigTransient = cell(2,2);      % corr speed-dfof / corr acceleration-dfof  - pvalue speed-dfof / pvalue acceleration-dfof
for ff = 1:numel(folders)
    load([folders{ff} '\abfFake.mat']);
    dd = dir([folders{ff} '\dfof_sig*']);
    load([folders{ff} '\' dd(1).name]); clear dd;
    load([folders{ff} '\sigTransInfo.mat']);
    
    yy = abfFake.y;
    speeds = ([0; diff(yy)]);       % can't normalize due to track ends
    acclrs = ([0; diff(speeds)]);   % should not be normalized
    speeds(speeds < -50) = nan;
    acclrs(acclrs < -50) = nan;
    acclrs(acclrs > 50) = nan;
    
    allSpeeds = [allSpeeds; speeds];
    allAcceleration = [allAcceleration; acclrs];
    df = nanmean(dfof_sig,2);
    allDfof = [allDfof; df];
    
    for ii = 1:numel(t_on_off)
        sigTrans = t_on_off{ii};
        if ~isnan(sigTrans)
            for jj = 1:size(sigTrans,1)
                indices = sigTrans(jj,2):sigTrans(jj,3);
                
                [correlationWithinSigTransient{1,1}(end+1), correlationWithinSigTransient{1,2}(end+1)] =  ...
                    corr( speeds(indices), df(indices) );
                [correlationWithinSigTransient{2,1}(end+1), correlationWithinSigTransient{2,2}(end+1)] =  ...
                    corr( acclrs(indices), df(indices) );
            end
        end
    end
end


figure;
plot(allSpeeds, allDfof, 'r.');
xlabel('Speed'); ylabel('\DeltaF/F');
set(gcf, 'Position', [200 200 400 350]);

figure;
plot(allAcceleration, allDfof, 'r.');
xlabel('Acceleration'); ylabel('\DeltaF/F');
set(gcf, 'Position', [600 200 400 350]);

figure;
plot(allAcceleration, allSpeeds, 'r.');
xlabel('Acceleration'); ylabel('Speed');
set(gcf, 'Position', [600 200 400 350]);


for target = 1:2
    switch target; case 1; ttl = 'SPEED'; case 2; ttl = 'ACCELERATION'; end;
    xx = correlationWithinSigTransient{target,1};
    
    figure;
    hold on
    histogram(xx);
    yl = get(gca, 'ylim');
    plot([nanmean(xx) nanmean(xx)], yl, 'r', 'LineWidth', 1.5);
    hold off
    [~, p] = ttest(correlationWithinSigTransient{1,1});
    title(['Corr - ' ttl ' and Dfof - ttest-p = ' num2str(p)]); xlabel('Correlation'); ylabel('Number');
end



%% CORRELATION vs Mean value

% steps = [50 40 30 24 20 10];
steps = [30 24 20];

stepWise = 2;
trackLength = 330;
% r1 = 120/stepWise; r2 = 140/stepWise;
correlationSpeedVsDfof = cell( 4, numel(steps) );           % corr speed-dfof / average speed / max speed - steps
correlationAccelerationVsDfof = cell( 4, numel(steps) );    % corr acceleration-dfof / average acceleration / max acceleration - steps

for nStep = 1:numel(steps)
    
    binWidth = steps(nStep)/stepWise;
    radius = round(binWidth/2);
    binWidth = binWidth-1;
    nBin = trackLength / stepWise;
    nn = nBin - binWidth;

    for ff = 1:numel(folders)
        load([folders{ff} '\abfFake.mat']);
        dd = dir([folders{ff} '\dfof_sig*']);
        load([folders{ff} '\' dd(1).name]); clear dd;
        
        yy = abfFake.y;
        speeds = ([0; diff(yy)]);       % can't normalize due to track ends
        acclrs = ([0; diff(speeds)]);   % should not be normalized
        speeds(speeds < -50) = nan;
        acclrs(acclrs < -50) = nan;
        acclrs(acclrs > 50) = nan;
%         acclrs(acclrs > prctile(acclrs, 95)) = nan;

        [~, idx] = histc(yy, 0 : stepWise : trackLength); 

        temp_dfof = nan(1, nBin); temp_speed = nan(1, nBin); temp_acclr = nan(1, nBin);
        for ii = 1:nBin
            indices = find(idx == ii);
            temp_dfof(ii) = nanmean(dfof_sig(indices, :), 'all');
            temp_speed(ii) = nanmean(speeds(indices));
            temp_acclr(ii) = nanmean(acclrs(indices));
        end

        corsM = nan( nBin, 1 );
        aveAcce = nan(nBin, 1);
        maxAcce = nan(nBin, 1);
        meanDF = nan(nBin, 1);
        kk = find(temp_acclr, 1);
        for ii = kk:nn
            [corsM(ii), ~] = corr( temp_acclr(ii:(ii+binWidth))', temp_dfof(ii:(ii+binWidth))' );
            aveAcce(ii ) = nanmean( temp_acclr(ii:(ii+binWidth)) );
            maxAcce(ii ) = max( temp_acclr(ii:(ii+binWidth)) );
            meanDF(ii ) = nanmean( temp_dfof(ii:(ii+binWidth)) );
        end
        correlationAccelerationVsDfof{1, nStep} = [correlationAccelerationVsDfof{1, nStep} corsM];
        correlationAccelerationVsDfof{2, nStep} = [correlationAccelerationVsDfof{2, nStep} aveAcce];
        correlationAccelerationVsDfof{3, nStep} = [correlationAccelerationVsDfof{3, nStep} maxAcce];
        correlationAccelerationVsDfof{4, nStep} = [correlationAccelerationVsDfof{4, nStep} meanDF];

        corsM = nan( nBin, 1 );
        aveSpeed = nan( nBin, 1 );
        maxSpeed = nan( nBin, 1 );
        meanDF = nan(nBin, 1);
        kk = find(temp_speed, 1);
        for ii = kk:nn
            [corsM(ii), ~] = corr( temp_speed(ii:(ii+binWidth))', temp_dfof(ii:(ii+binWidth))' );
            aveSpeed(ii) = nanmean( temp_speed(ii:(ii+binWidth)) );
            maxSpeed(ii) = max( temp_speed(ii:(ii+binWidth)) );
            meanDF(ii) = nanmean( temp_dfof(ii:(ii+binWidth)) );
        end
        correlationSpeedVsDfof{1, nStep} = [correlationSpeedVsDfof{1, nStep} corsM];
        correlationSpeedVsDfof{2, nStep} = [correlationSpeedVsDfof{2, nStep} aveSpeed];
        correlationSpeedVsDfof{3, nStep} = [correlationSpeedVsDfof{3, nStep} maxSpeed];
        correlationSpeedVsDfof{4, nStep} = [correlationSpeedVsDfof{4, nStep} meanDF];

    end
end


%%
n = numel(steps);
for types = 1%:2
    switch types
        case 1
            useData = correlationAccelerationVsDfof;
            ttl1 = 'ACCELERATION';
        case 2
            useData = correlationSpeedVsDfof;
            ttl1 = 'SPEED';
    end
    fig = figure;
    for target = 2:3
        switch target; case 2; ttl2 = 'mean'; case 3; ttl2 = 'max'; end;
        for nStep = 1:numel(steps);
            subplot(2,n, (target-2)*n + nStep);
            hold on
            xx = useData{1, nStep};
            xx = reshape(xx, 1, []);
            yy = useData{target, nStep};
            yy = reshape(yy, 1, []);
            plot(yy, xx, 'g.');
            yl = get(gca, 'ylim');
            plot([prctile(yy, 95) prctile(yy, 95)], yl, 'r');
            hold off
            indices = intersect( find(~isnan(xx)), find(~isnan(yy)) );
            [c, p] = corr( xx(indices)', yy(indices)' );
            c = c > 0;
            title([ttl2 ' - ' num2str(steps(nStep)) ' - ' num2str(round(c),2) ' - ' num2str(round(p),2)]);
        end
    end
    sgtitle('Different section length');
    han=axes(fig,'visible','off');
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han, 'Correlation', 'FontSize', 15);
    xlabel(han, ttl1, 'FontSize', 15);
end


useData = correlationAccelerationVsDfof;
n = 3;
xx = useData{1, n};
xx = reshape(xx, 1, []);
yy = useData{2, n};
yy = reshape(yy, 1, []);
indices = intersect( find(~isnan(xx)), find(~isnan(yy)) );
xx = xx(indices);
yy = yy(indices);

figure;
hold on
plot(yy, xx, 'r.');
[b,dev,stats] = glmfit(xx, yy, 'normal');
xl = get(gca, 'xlim');
mm = xl(1):range(xl)/10:xl(2);
nn = mm*b(2) + b(1);
plot(nn, mm, 'y', 'LineWidth', 3);
hold off


%%
nStep = 2;
xx = correlationAccelerationVsDfof{1, nStep};
yy = correlationAccelerationVsDfof{2, nStep};

figure;
plot(yy, xx, 'k.'); xlim([-0.2 0.6]);

edges = [-0.2 -0.1:0.02:0.2 0.4 0.6];
xidx = (edges(2:end) + edges(1:end-1))/2;
nAccBin = length(edges)-1;

nFov = size(xx,2);
meanDfof = nan(8, nFov);
for fov = 1:nFov
    xx1 = xx(:,fov);
    yy1 = yy(:,fov);
    
    [~, idx] = histc(yy1, edges);
    for ii = 1:nAccBin
        indices = find(idx == ii);
        meanDfof(ii, fov) = nanmean(xx1(indices));
    end
end
figure;
plot(xidx, meanDfof, '.')
% set(gca, 'xtick', [0:nAccBin]+0.5, 'xticklabel', edges);
xlabel('Acceleration'); ylabel('Correlation');
% xlim([0 nAccBin+1]);

figure;
plot(xidx, meanDfof)
% set(gca, 'xtick', [0:nAccBin]+0.5, 'xticklabel', edges);
xlabel('Acceleration'); ylabel('Correlation');
% xlim([0 nAccBin+1]);

figure;
hold on
plot(xidx, nanmean(meanDfof, 2), 'LineWidth', 3);
plot(xidx, nanmean(meanDfof(:, randperm(30, 20)), 2));
plot(xidx, nanmean(meanDfof(:, randperm(30, 20)), 2));
plot(xidx, nanmean(meanDfof(:, randperm(30, 20)), 2));
plot(xidx, nanmean(meanDfof(:, randperm(30, 20)), 2));
plot(xidx, nanmean(meanDfof(:, randperm(30, 20)), 2));
plot(xidx, nanmean(meanDfof(:, randperm(30, 20)), 2));
plot(xidx, nanmean(meanDfof(:, randperm(30, 20)), 2));
plot(xidx, nanmean(meanDfof(:, randperm(30, 20)), 2));
hold off
xlabel('Average acceleration per spatial bin'); ylabel('Correlation accleration-dfof');
set(gcf, 'Position', [200 200 700 350]);



%% PLOT ACCELERATION VS CORRELATION-ACCEL-DFOF

nStep = 2;
xx = correlationAccelerationVsDfof{1, nStep};
yy = correlationAccelerationVsDfof{2, nStep};
% edges = [-0.2 -0.1:0.02:0.2 0.3 0.4 0.5];
% edges = [-0.4:0.04:0.84];
edges = linspace(-0.4, 0.84, 40);
xidx = (edges(2:end) + edges(1:end-1))/2;
nAccBin = length(edges)-1;

nFov = size(xx,2);
meanDfof = nan(8, nFov);
meanYY = nan(8, nFov);
for fov = 1:nFov
    xx1 = xx(:,fov);
    yy1 = yy(:,fov);
    
    [~, idx] = histc(yy1, edges);
%     [~, ~, idx] = histcounts(yy1, prctile(yy1,0:5:100));
    for ii = 1:nAccBin
        indices = find(idx == ii);
        meanDfof(ii, fov) = nanmean(xx1(indices));
        meanYY(ii, fov) = nanmean(yy1(indices));
    end
end

figure; hold on
errorbar(xidx*10, nanmean(meanDfof, 2), nanstd(meanDfof, 0, 2)/sqrt(30), 'c', 'LineWidth', 1.5);
% errorbar(nanmean(meanYY, 2)*10, nanmean(meanDfof, 2), nanstd(meanDfof, 0, 2)/sqrt(30), 'c', 'LineWidth', 1.5);
yl = get(gca, 'ylim');
plot([1.1 1.1], yl, 'k--');
xlabel('Acceleration (cm/s^{2})'); ylabel('Corr. accleration-\DeltaF/F');
set(gcf, 'Position', [123 123 345 222]);
ylim(yl);


xx = reshape(xx, 1, []);
yy = reshape(yy, 1, []);

n = 20;
edges = prctile(yy, 0:(100/n):100);
xidx = (edges(2:end) + edges(1:end-1))/2;
[~,~,binIdx] = histcounts(yy, edges);
mmXX = nan(1,n); seXX = nan(1,n);
mmYY = nan(1,n); seYY = nan(1,n);
for ii = 1:n
    [mmXX(ii), seXX(ii), ~, ~] = getMeanAndSE( xx(binIdx == ii), 'all');
    [mmYY(ii), seYY(ii), ~, ~] = getMeanAndSE( yy(binIdx == ii)*10, 'all');
end

figure; hold on
errorbar(mmYY, mmXX, seXX, 'c', 'LineWidth', 1.5);
% errorbar(xidx*10, mmXX, seXX, 'c', 'LineWidth', 1.5);
% errorbar(mmXX*10, mmYY, seYY, 'c', 'LineWidth', 1.5);
yl = get(gca, 'ylim');
plot([1.1 1.1], yl, 'k--');
xlabel('Acceleration (cm/s^{2})'); ylabel('Corr. accleration-\DeltaF/F');
set(gcf, 'Position', [123 123 345 222]);


%%
xx = correlationSpeedVsDfof{1, nStep};
yy = correlationSpeedVsDfof{2, nStep};
edges = [-0.2 0:0.5:8];
xidx = (edges(2:end) + edges(1:end-1))/2;
nAccBin = length(edges)-1;

nFov = size(xx,2);
meanDfof = nan(8, nFov);
for fov = 1:nFov
    xx1 = xx(:,fov);
    yy1 = yy(:,fov);
    
    [~, idx] = histc(yy1, edges);
    for ii = 1:nAccBin
        indices = find(idx == ii);
        meanDfof(ii, fov) = nanmean(xx1(indices));
    end
end

figure;
plot(xidx, nanmean(meanDfof, 2));
xlabel('Speed'); ylabel('Correlation');

%% Histogram of accelerations and speed

stepWise = 2;
nBin = 210/stepWise;
r1 = 120/stepWise; r2 = 140/stepWise;

allAccelerations = [];
allAccelerationsSmooth = [];
allSpeeds = [];
for ff = 1:numel(folders)
    load([folders{ff} '\abfFake.mat']);

    yy = abfFake.y;
    speeds = ([0; diff(yy)]);       % can't normalize due to track ends
    acclrs = ([0; diff(speeds)]);   % should not be normalized
    speeds(speeds < -50) = nan;
    acclrs(acclrs < -50) = nan;
    acclrs(acclrs > 50) = nan;

    allAccelerations = [allAccelerations; acclrs];
    allSpeeds = [allSpeeds; speeds];
    
    [~, idx] = histc(yy, 0 : stepWise : 210); 

    temp_acclr = nan(1, nBin);
    for ii = 1:nBin
        indices = find(idx == ii);
        temp_acclr(ii) = nanmean(acclrs(indices));
    end
    allAccelerationsSmooth = [allAccelerationsSmooth; temp_acclr];
end

figure;
nStd = 1;
hold on
xx = allAccelerations;
% xx = allAccelerationsSmooth; xx = reshape(xx, 1, []);
histogram(xx);
yl = get(gca, 'ylim');
p(1) = plot([nanmean(xx) nanmean(xx)], yl, 'r', 'LineWidth', 1.5);
p(2) = plot([nanmean(xx)-nStd*nanstd(xx) nanmean(xx)-nStd*nanstd(xx)], yl, 'y', 'LineWidth', 1.5);
plot([nanmean(xx)+nStd*nanstd(xx) nanmean(xx)+nStd*nanstd(xx)], yl, 'y', 'LineWidth', 1.5);
hold off
% xlim([-1 1]);
set(gcf, 'Position', [200 200 600 380]);
legend(p, 'mean', '1-std', 'orientation', 'horizontal');
xlabel('ACCELERATIONS');

figure;
nStd = 1;
hold on
xx = allSpeeds;
histogram(xx);
% boxplot(xx, 'orientation', 'horizontal');
yl = get(gca, 'ylim');
p(1) = plot([nanmean(xx) nanmean(xx)], yl, 'r', 'LineWidth', 1.5);
p(2) = plot([nanmean(xx)-nStd*nanstd(xx) nanmean(xx)-nStd*nanstd(xx)], yl, 'y', 'LineWidth', 1.5);
plot([nanmean(xx)+nStd*nanstd(xx) nanmean(xx)+nStd*nanstd(xx)], yl, 'y', 'LineWidth', 1.5);
hold off
% xlim([-1 1]);
set(gcf, 'Position', [200 200 600 380]);
legend(p, 'mean', '1-std', 'orientation', 'horizontal');
xlabel('SPEEDS');

% figure;
% xx = allAccelerationsSmooth; xx = reshape(xx, 1, []); thresh = nanmean(xx)+nStd*nanstd(xx);
% xx = allAccelerationsSmooth;
% xx(xx > thresh) = nan;
% hold on
% errorbar(1:nBin, nanmean(xx, 1), nanstd(xx, 1), 'b');
% plot(1:nBin, nanmean(xx, 1), 'LineWidth', 2);
% yl = get(gca, 'ylim');
% plot([r1 r1], yl, 'r');
% plot([r2 r2], yl, 'r');
% hold off
% set(gca, 'xtick', [0 nBin]+0.5, 'xticklabel', [0 210]); xlim([0 nBin]+0.5);
% set(gcf, 'Position', [200 200 600 380]);
% xlim([0 nBin] + 0.5); ylim([-0.2 0.3]);
% xlabel('Distance'); ylabel('Acceleration');



%% CORRELATION vs Mean value
% ID05

stepWise = 2;
trackLength = 330;
areaWidth = 24;

ID05_folders{1} = 'H:\3 - audioVR_NEW_learning\ID20211105\20220201\loc1\audioVR_NEW\pcaica';
ID05_folders{2} = 'H:\3 - audioVR_NEW_learning\ID20211105\20220202\loc3\audioVR_NEW\pcaica';
ID05_folders{3} = 'H:\3 - audioVR_NEW_learning\ID20211105\20220203\loc2\audioVR_NEW\pcaica';
ID05_folders{4} = 'H:\3 - audioVR_NEW_learning\ID20211105\20220204\loc3\audioVR_NEW\pcaica';
ID05_folders{5} = 'H:\3 - audioVR_NEW_learning\ID20211105\20220205\loc1\audioVR_NEW\pcaica';
ID05_folders{6} = 'H:\3 - audioVR_NEW_learning\ID20211105\20220205\loc2\audioVR_NEW\pcaica';

correlationAccelerationVsDfof_ID05 = cell(4,1);
binWidth = areaWidth/stepWise;
radius = round(binWidth/2);
binWidth = binWidth-1;
nBin = trackLength / stepWise;
nn = nBin - binWidth;

for ff = 1:numel(ID05_folders)
    load([ID05_folders{ff} '\abfFake.mat']);
    dd = dir([ID05_folders{ff} '\dfof_sig*']);
    load([ID05_folders{ff} '\' dd(1).name]); clear dd;

    yy = abfFake.y;
    speeds = ([0; diff(yy)]);       % can't normalize due to track ends
    acclrs = ([0; diff(speeds)]);   % should not be normalized
    speeds(speeds < -50) = nan;
    acclrs(acclrs < -50) = nan;
    acclrs(acclrs > 50) = nan;
%         acclrs(acclrs > prctile(acclrs, 95)) = nan;

    [~, idx] = histc(yy, 0 : stepWise : trackLength); 

    temp_dfof = nan(1, nBin); temp_speed = nan(1, nBin); temp_acclr = nan(1, nBin);
    for ii = 1:nBin
        indices = find(idx == ii);
        temp_dfof(ii) = nanmean(dfof_sig(indices, :), 'all');
        temp_speed(ii) = nanmean(speeds(indices));
        temp_acclr(ii) = nanmean(acclrs(indices));
    end

    corsM = nan( nBin, 1 );
    aveAcce = nan(nBin, 1);
    maxAcce = nan(nBin, 1);
    meanDF = nan(nBin, 1);
    kk = find(temp_acclr, 1);
    for ii = kk:nn
        [corsM(ii), ~] = corr( temp_acclr(ii:(ii+binWidth))', temp_dfof(ii:(ii+binWidth))' );
        aveAcce(ii ) = nanmean( temp_acclr(ii:(ii+binWidth)) );
        maxAcce(ii ) = max( temp_acclr(ii:(ii+binWidth)) );
        meanDF(ii ) = nanmean( temp_dfof(ii:(ii+binWidth)) );
    end
    correlationAccelerationVsDfof_ID05{1} = [correlationAccelerationVsDfof_ID05{1} corsM];
    correlationAccelerationVsDfof_ID05{2} = [correlationAccelerationVsDfof_ID05{2} aveAcce];
    correlationAccelerationVsDfof_ID05{3} = [correlationAccelerationVsDfof_ID05{3} maxAcce];
    correlationAccelerationVsDfof_ID05{4} = [correlationAccelerationVsDfof_ID05{4} meanDF];
end




%%
nStep = 2;
xx = correlationAccelerationVsDfof_ID05{1};
yy = correlationAccelerationVsDfof_ID05{2};

figure;
plot(yy, xx, 'k.'); xlim([-0.2 0.6]);

edges = [-0.2 -0.1:0.02:0.2 0.4 0.6]/10;
xidx = (edges(2:end) + edges(1:end-1))/2;
nAccBin = length(edges)-1;

nFov = size(xx,2);
meanDfof = nan(8, nFov);
for fov = 1:nFov
    xx1 = xx(:,fov);
    yy1 = yy(:,fov);
    
    [~, idx] = histc(yy1, edges);
    for ii = 1:nAccBin
        indices = find(idx == ii);
        meanDfof(ii, fov) = nanmean(xx1(indices));
    end
end
figure;
plot(xidx, meanDfof, '.')
% set(gca, 'xtick', [0:nAccBin]+0.5, 'xticklabel', edges);
xlabel('Acceleration'); ylabel('Correlation');
% xlim([0 nAccBin+1]);

figure;
plot(xidx, meanDfof)
% set(gca, 'xtick', [0:nAccBin]+0.5, 'xticklabel', edges);
xlabel('Acceleration'); ylabel('Correlation');
% xlim([0 nAccBin+1]);

figure;
plot(xidx, nanmean(meanDfof, 2));
xlabel('Acceleration'); ylabel('Correlation');



%% ACCELERATION REMOVED
% Dfof after
% Correlation along track after

AccThreshs = [inf 0.8 1 1.1 1.2 1.5];
% AccThreshs = 1:10;

stepWise = 5;
trackLength = 330;
nBin = trackLength / stepWise;
areasWith = 25;
binWidth = areasWith/stepWise;
Radius = round(binWidth/2);
binWidth = binWidth-1;
nn = nBin - binWidth;

newDfof = cell( 1, numel(AccThreshs) );
newAcc = cell( 1, numel(AccThreshs) );
newCors = cell( 1, numel(AccThreshs) );
newPrctiles = cell( 1, numel(AccThreshs) );

for nThresh = 1:numel(AccThreshs)
    thresh = AccThreshs(nThresh);
%     thresh = 0.8;
    
    for ff = 1:numel(folders)
%     for ff = randperm(88, 15)
        load([folders{ff} '\abfFake.mat']);
        dd = dir([folders{ff} '\dfof_sig*']);
        load([folders{ff} '\' dd(1).name]); clear dd;

        yy = abfFake.y;
        speeds = ([0; diff(yy)]);       % can't normalize due to track ends
        acclrs = ([0; diff(speeds)]);   % should not be normalized
        speeds(speeds < -50) = nan;
        acclrs(acclrs < -50) = nan;
        acclrs(acclrs > 50) = nan;
        
        indexes = find(acclrs > thresh);
        pct = length(indexes) / length(acclrs(~isnan(acclrs)));
        newPrctiles{nThresh} = [newPrctiles{nThresh} pct];
        acclrs(indexes) = nan;
        dfof_sig(indexes, :) = nan;
        
        [~, idx] = histc(yy, 0 : stepWise : trackLength); 

        temp_dfof = nan(nBin, size(dfof_sig, 2) ); temp_acclr = nan(1, nBin);
        for ii = 1:nBin
            indices = find(idx == ii);
            temp_acclr(ii) = nanmean(acclrs(indices));
            for nn = 1:size(dfof_sig, 2)
                temp_dfof(ii, nn) = nanmean(dfof_sig(indices, nn), 'all');
            end
        end
        
        newDfof{nThresh} = [newDfof{nThresh} temp_dfof];
        newAcc{nThresh} = [newAcc{nThresh}; temp_acclr];
        
        temp_dfof = nanmean(temp_dfof, 2)';
        corsM = nan( 1, nBin );
        kk = find(temp_acclr, 1);
        for ii = kk:nn
            [corsM(ii+Radius), ~] = corr( temp_acclr(ii:(ii+binWidth))', temp_dfof(ii:(ii+binWidth))' );
        end
        newCors{nThresh} = [newCors{nThresh}; corsM];
    end
end


for data = [3 2 1]
    clear pp;
    switch data
        case 1
            useData = newAcc;
            ylbl = 'ACCELERATION';
        case 2
            useData = newDfof;
            ylbl = '\DeltaF/F';
        case 3
            useData = newCors;
            ylbl = 'CORRELATION';
    end
    figure;
    hold on
    for nThresh = 1:numel(AccThreshs)
        xx = useData{nThresh};
        pp(nThresh) = plot(nanmean(xx, 1), 'LineWidth', 1.5);
    end
    if data == 1
        xx = cat(2, newPrctiles{:})*100;
        ttl = [num2str(mean(xx)) ' +- ' num2str(std(xx)) '% removed'];
    end
    hold off
%     legend(pp, 'control', '0.5', '0.8', '1', '1.12', '1.5', 'orientation', 'horizontal', 'location', 'south');
    set(gca, 'xtick', [0 nBin]+0.5, 'xticklabel', [0 210]); xlim([0 nBin]+0.5);
    set(gcf, 'Position', [200 200 600 350]);
    xlim([0 nBin] + 0.5);
    xlabel('Distance'); ylabel(ylbl); if data == 1; title(ttl); end
end
% saveas(figure(1), 'fig1.jpg');
% saveas(figure(2), 'fig2.jpg');
% saveas(figure(3), 'fig3.jpg');



%% SPEED REMOVED
% Dfof after
% Correlation along track after

% SpdThreshs = [-inf 0.05 0.2 0.5 1 1.2];
SpdThreshs = 1:1;

stepWise = 2;
trackLength = 210;
nBin = trackLength / stepWise;
areasWith = 24;
binWidth = areasWith/stepWise;
Radius = round(binWidth/2);
binWidth = binWidth-1;
nn = nBin - binWidth;

newDfof = cell( 1, numel(SpdThreshs) );
newSpeeds = cell( 1, numel(SpdThreshs) );
newCors = cell( 1, numel(SpdThreshs) );
newPrctiles = cell( 1, numel(SpdThreshs) );

for nThresh = 1:numel(SpdThreshs)
%     thresh = SpdThreshs(nThresh);
    thresh = -inf;
    
    for ff = 1:numel(folders)
%     for ff = randperm(88, 15)
        load([folders{ff} '\abfFake.mat']);
        dd = dir([folders{ff} '\dfof_sig*']);
        load([folders{ff} '\' dd(1).name]); clear dd;

        yy = abfFake.y;
        speeds = ([0; diff(yy)]);       % can't normalize due to track ends
        speeds(speeds < -50) = nan;
        
        indexes = find(speeds < thresh);
        pct = length(indexes) / length(speeds(~isnan(speeds)));
        newPrctiles{nThresh} = [newPrctiles{nThresh} pct];
        speeds(indexes) = nan;
        dfof_sig(indexes, :) = nan;
        
        [~, idx] = histc(yy, 0 : stepWise : trackLength); 

        temp_dfof = nan(1, nBin); temp_speed = nan(1, nBin);
        for ii = 1:nBin
            indices = find(idx == ii);
            temp_dfof(ii) = nanmean(dfof_sig(indices, :), 'all');
            temp_speed(ii) = nanmean(speeds(indices));
        end
        
        newDfof{nThresh} = [newDfof{nThresh}; temp_dfof];
        newSpeeds{nThresh} = [newSpeeds{nThresh}; temp_speed];
        
        corsM = nan( 1, nBin );
        kk = find(temp_speed, 1);
        for ii = kk:nn
            [corsM(ii+Radius), ~] = corr( temp_speed(ii:(ii+binWidth))', temp_dfof(ii:(ii+binWidth))' );
        end
        newCors{nThresh} = [newCors{nThresh}; corsM];
    end
end

for data = [3 2 1]
    clear pp;
    switch data
        case 1
            useData = newSpeeds;
            ylbl = 'SPEED';
        case 2
            useData = newDfof;
            ylbl = '\DeltaF/F';
        case 3
            useData = newCors;
            ylbl = 'CORRELATION';
    end
    figure;
    hold on
    for nThresh = 1:numel(SpdThreshs)
        xx = useData{nThresh};
        pp(nThresh) = plot(nanmean(xx, 1), 'LineWidth', 1.5);
    end
    if data == 1
        xx = cat(2, newPrctiles{:})*100;
        ttl = [num2str(mean(xx)) ' +- ' num2str(std(xx)) '% removed'];
    end
    hold off
%     legend(pp, 'control', '0.05', '0.2', '0.5', '1', '1.2', 'orientation', 'horizontal', 'location', 'south');
    set(gca, 'xtick', [0 nBin]+0.5, 'xticklabel', [0 210]); xlim([0 nBin]+0.5);
    set(gcf, 'Position', [200 200 600 350]);
    xlim([0 nBin] + 0.5);
    xlabel('Distance'); ylabel(ylbl); if data == 1; title(ttl); end
end
% saveas(figure(1), 'fig1.jpg');
% saveas(figure(2), 'fig2.jpg');
% saveas(figure(3), 'fig3.jpg');


%% BOTH ACCELERATION AND SPEED REMOVED
% Dfof after
% Correlation along track after

% left = accleration-threshold
% right = speed-threshold
% shuffles = [inf -inf;
%     1.12 -inf;
%     inf 0.04;
%     1.12 0.04;
%     0.8 1;
%     1.12 0.1];
shuffles = [inf -inf;
    1 0.12;
    1.1 0.12];
nShuffle = size(shuffles, 1);

stepWise = 5;
trackLength = 330;
nBin = trackLength / stepWise;
areasWith = 25;
binWidth = areasWith/stepWise;
Radius = round(binWidth/2);
binWidth = binWidth-1;
nn = nBin - binWidth;

newDfof = cell( 1, nShuffle );
newSpeeds = cell( 1, nShuffle );
newAccs = cell( 1, nShuffle );
newCorsAcc = cell( 1, nShuffle );
newCorsSpeed = cell( 1, nShuffle );
newPrctiles = cell( 1, nShuffle );

for nThresh = 1:nShuffle
    threshA = shuffles(nThresh, 1);
    threshS = shuffles(nThresh, 2);
    
%     if nThresh == 1
        fovs = 1:numel(folders);
%     else
%         fovs = randperm(88, 15);
%     end
    for ff = fovs
        load([folders{ff} '\abfFake.mat']);
        dd = dir([folders{ff} '\dfof_sig*']);
        load([folders{ff} '\' dd(1).name]); clear dd;

        yy = abfFake.y;
        speeds = ([0; diff(yy)]);       % can't normalize due to track ends
        acclrs = ([0; diff(speeds)]);   % should not be normalized
        speeds(speeds < -50) = nan;
        acclrs(acclrs < -50) = nan;
        acclrs(acclrs > 50) = nan;
        
        indexesA = find(acclrs > threshA);
        indexesS = find(speeds < threshS);
        indexes = unique([indexesA; indexesS]);
        pct = length(indexes) / length(speeds(~isnan(speeds)));
        newPrctiles{nThresh} = [newPrctiles{nThresh} pct];
        speeds(indexes) = nan;
        acclrs(indexes) = nan;
        dfof_sig(indexes, :) = nan;
        
        [~, idx] = histc(yy, 0 : stepWise : trackLength); 

        temp_dfof = nan( 1, nBin); temp_speed = nan(1, nBin); temp_acclr = nan(1, nBin);
        for ii = 1:nBin
            indices = find(idx == ii);
            temp_speed(ii) = nanmean(speeds(indices));
            temp_acclr(ii) = nanmean(acclrs(indices));
            temp_dfof(ii) = nanmean(dfof_sig(indices, :), 'all');

        end
        
        newDfof{nThresh} = [newDfof{nThresh}; temp_dfof];
        newSpeeds{nThresh} = [newSpeeds{nThresh}; temp_speed];
        newAccs{nThresh} = [newAccs{nThresh}; temp_acclr];
        
        temp_dfof = nanmean(temp_dfof, 1);
        
%         corsM = nan( 1, nBin );
%         kk = find(temp_speed, 1);
%         for ii = kk:nn
%             [corsM(ii+Radius), ~] = corr( temp_speed(ii:(ii+binWidth))', temp_dfof(ii:(ii+binWidth))' );
%         end
%         newCorsSpeed{nThresh} = [newCorsSpeed{nThresh}; corsM];
%         
%         corsM = nan( 1, nBin );
%         kk = find(temp_acclr, 1);
%         for ii = kk:nn
%             [corsM(ii+Radius), ~] = corr( temp_acclr(ii:(ii+binWidth))', temp_dfof(ii:(ii+binWidth))' );
%         end
%         newCorsAcc{nThresh} = [newCorsAcc{nThresh}; corsM];
    end
end


for data = 3%1:3
    switch data
        case 1
            useData = newCorsAcc;
            ttl = 'ACCELERATION';
            ylbl = 'Correlation';
        case 2
            useData = newCorsSpeed;
            ttl = 'SPEED';
            ylbl = 'Correlation';
        case 3
            useData = newDfof;
            ttl = '\DeltaF/F';
            ylbl = '\DeltaF/F';
    end
    clear pp;
    lgnds = [];
    figure;
    hold on
    for nThresh = 1:nShuffle
        xx = useData{nThresh};
        pp(nThresh) = plot(nanmean(xx, 1), 'LineWidth', 1.2);
        if nThresh == 1
            lgnds{nThresh} = 'control';
        else
            lgnds{nThresh} = [num2str(shuffles(nThresh,1)) ' - ' num2str(shuffles(nThresh,2)) ' - ' num2str(round(mean(newPrctiles{nThresh}*100),2))];
        end
    end
    plot([0 nBin] + 0.5, [0 0], 'k--');
    hold off
    legend(pp, lgnds, 'orientation', 'vertical', 'location', 'eastoutside');
    set(gca, 'xtick', [0 nBin]+0.5, 'xticklabel', [0 trackLength]); xlim([0 nBin]+0.5);
    set(gcf, 'Position', [200 200 800 350]);
    xlim([0 nBin] + 0.5);
    xlabel('Distance'); ylabel(ylbl, 'FontSize', 15); title(ttl);
end




%% BOTH ACCELERATION AND SPEED REMOVED
% Dfof after
% Correlation along track after
% REALISTIC VERSION

% left = accleration-threshold
% right = speed-threshold
% shuffles = [inf -inf;
%     1.12 -inf;
%     inf 0.04;
%     1.12 0.04;
%     0.8 1;
%     1.12 0.1];
shuffles = [inf -inf;
    1 0.12;
    1.1 0.12];
nShuffle = size(shuffles, 1);

stepWise = 5;
trackLength = 330;
nBin = trackLength / stepWise;
areasWith = 25;
binWidth = areasWith/stepWise;
Radius = round(binWidth/2);
binWidth = binWidth-1;
nn = nBin - binWidth;

newDfof = cell( 1, nShuffle );

for nThresh = 1:nShuffle
    threshA = shuffles(nThresh, 1);
    threshS = shuffles(nThresh, 2);
    
    fovs = 1:numel(folders);
    for ff = fovs
        load([folders{ff} '\abfFake.mat']);
        dd = dir([folders{ff} '\dfof_sig*']);
        load([folders{ff} '\' dd(1).name]); clear dd;

        % GET ACCELERATION
        yy = abfFake.y;
        speeds = ([0; diff(yy)]);       % can't normalize due to track ends
        acclrs = ([0; diff(speeds)]);   % should not be normalized
        acclrs(acclrs < -50) = nan;
        acclrs(acclrs > 50) = nan;
        
        % REMOVE WHERE ACCELERATION ABOVE THRESHOLD
        indexes = find(acclrs > threshA);
        abf = abfFake;
        abf.t(indexes) = [];
        abf.y(indexes) = [];
        abf.imageIndex = 1:length(abf.y);
        dfof_sig(indexes, :) = [];
        
        %%% GET DFOFM
        sz1 = size(dfof_sig, 1);    % n data points
        sz2 = size(dfof_sig, 2);    % n cells
        dfofM = getRunByRunActivity(1, sz2, 1, sz1, 1, sz1, ...
                    5, 0, trackLength, threshS, dfof_sig, abf);
                
        dfofaveragesmooth = nan( trackLength/5, sz2 );
        
        for nCell = 1:size(dfof_sig, 2)
            if ~isempty(dfofM{nCell})
                % get dfofaveragesmooth
                dfofaveragesmooth(:, nCell) = nanmean( dfofM{nCell}, 1 )';
            end
        end
        
        newDfof{nThresh} = [newDfof{nThresh} dfofaveragesmooth];        
    end
end


useData = newDfof;
ttl = '\DeltaF/F';
ylbl = '\DeltaF/F';

clear pp;
lgnds = [];
figure;
hold on
for nThresh = 1:nShuffle
    xx = useData{nThresh};
    pp(nThresh) = plot(nanmean(xx, 2), 'LineWidth', 1.2);
    if nThresh == 1
        lgnds{nThresh} = 'control';
    else
        lgnds{nThresh} = [num2str(shuffles(nThresh,1)) ' - ' num2str(shuffles(nThresh,2))];
    end
end
plot([0 nBin] + 0.5, [0 0], 'k--');
hold off
legend(pp, lgnds, 'orientation', 'vertical', 'location', 'eastoutside');
set(gca, 'xtick', [0 nBin]+0.5, 'xticklabel', [0 trackLength]); xlim([0 nBin]+0.5);
set(gcf, 'Position', [200 200 800 350]);
xlim([0 nBin] + 0.5);
xlabel('Distance'); ylabel(ylbl, 'FontSize', 15); title(ttl);



%%
