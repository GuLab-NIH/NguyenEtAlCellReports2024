
% load data
% contact Lead Author for data


%% plot distribution of correlation scores

figure
for n=1:length(c)
    subplot(3,3,n);
    [x,p]= ksdensity(c{n},'width',0.05);
hold on
plot(p,x,'k')
[x,p]= ksdensity(cs{n},'width',0.05);
hold on
plot(p,x,'k--')
[r1,p1]=kstest2(c{n},cs{n},'tail','larger');
[r2,p2]=kstest2(c{n},cs{n});
% [r,p]=kstest2(c{n},cs{n});
title(['p1=',num2str(p1),' p2=',num2str(p2)]);
end

%from top left to top right: avr2 twice, avr2 and vvr2, vvr2 twice
%from bottom left to top right: avr1 twice, avr1 and vvr1, vvr1 twice


%% plotting visual twice and visual audio

figure
%visual new twice
[x,p]= ksdensity(c{3},'width',0.05);
hold on
plot(p,x,'g')
%visual and overlap
[x,p]= ksdensity(c{9},'width',0.05);
hold on
plot(p,x,'m')
[r,p]=kstest2(c{3},c{9},'tail','smaller');
title(['visual twice and visual overlap',num2str(p)]);

% saveas(gcf,'visualTwiceVisualOV.fig');

%% what info do cell encode if they are very similar


allCells=[];
allFieldsInfo=[];

clear allDfof
load('AnewTwice.mat');
allCells{1}=allDfof;
allFieldsInfo{1}=fieldsInfo;

clear allDfof
load('AnewVnew.mat');
allCells{2}=allDfof;
allFieldsInfo{2}=fieldsInfo;

clear allDfof
load('VnewTwice.mat');
allCells{3}=allDfof;
allFieldsInfo{3}=fieldsInfo;

clear allDfof
load('AoldTwice.mat');
allCells{4}=allDfof;
allFieldsInfo{4}=fieldsInfo;

clear allDfof
load('AoldVold.mat');
allCells{5}=allDfof;
allFieldsInfo{5}=fieldsInfo;

clear allDfof
load('VoldTwice.mat');
allCells{6}=allDfof;
allFieldsInfo{6}=fieldsInfo;

clear allDfof
load('OVTwice.mat');
allCells{7}=allDfof;
allFieldsInfo{7}=fieldsInfo;

clear allDfof
load('OVAnew.mat');
allCells{8}=allDfof;
allFieldsInfo{8}=fieldsInfo;

clear allDfof
load('OVVnew.mat');
allCells{9}=allDfof;
allFieldsInfo{9}=fieldsInfo;

save('allCells.mat','allCells');
save('allFieldsInfo.mat','allFieldsInfo');


%%
load('c.mat');
load('cs.mat');
load('allCells.mat');
load('allFieldsInfo.mat');

%for all FOVs, extrac the cell pairs > 95 percentile of shuffle
similarCells={};
similarCells_fields={};
idxSimilar=[];
for n=1:length(c);
    thresh=prctile(cs{n},95);
    i=find(c{n}>thresh);
    idxSimilar{n}=i;
    similarCells{n}=allCells{n}(i);
    similarCells_fields{n}=allFieldsInfo{n}(i);
end

similarCellsFOV={};
similarCellsFOV_fields={};
for n=1:length(similarCells);
    similarCellsFOV{n}=[];
    similarCellsFOV_fields{n}=[];
    for m=1:length(similarCells{n});
        similarCellsFOV{n}(:,m,1)=similarCells{n}{m}(:,1);
        similarCellsFOV{n}(:,m,2)=similarCells{n}{m}(:,2);
        
        similarCellsFOV_fields{n}(:,m,1)=similarCells_fields{n}{m}(:,1);
        similarCellsFOV_fields{n}(:,m,2)=similarCells_fields{n}{m}(:,2);
    end
end




%% plot

load('idxSimilar.mat');
load('allFieldsInfo.mat');

figure;
for n=1:length(idxSimilar);
    subplot(3,3,n);
    
    simCells = cat(2, allFieldsInfo{n}{idxSimilar{n}});
    diffCells = cat(2, allFieldsInfo{n}{setdiff(1:length(allFieldsInfo{n}), idxSimilar{n})});
    
    simCells = simCells(:, 1:2:end);
    diffCells = diffCells(:, 1:2:end);
    
    hold on
    plot(mean(simCells, 2), 'm');
    plot(mean(diffCells, 2), 'k');
    hold on
end
saveas(gcf,'allFields.fig');


%% MORE PLOTS

% load('AoldTwice.mat'); iidx = [2 4]; xlb = 'AVR1'; ylb = 'AVR1';
% load('AnewTwice.mat'); iidx = [2 4]; xlb = 'AVR2'; ylb = 'AVR2';
% load('VnewTwice.mat'); iidx = [1 2]; xlb = 'VVR2'; ylb = 'VVR2';
% load('VoldTwice.mat'); iidx = [1 2]; xlb = 'VVR1'; ylb = 'VVR1';

figure; hold on
xx = cueScores(:, iidx(1));
yy = cueScores(:, iidx(2));
isnanIndices = unique([ find(isnan(xx)); find(isnan(yy)) ]);
xx(isnanIndices) = []; yy(isnanIndices) = [];

[bins, ~, ind] = histcounts(xx, prctile(xx, 0:10:100));
aves = nan(1, length(bins)); mm = nan(1, length(bins)); se = nan(1, length(bins));
for ii = 1:length(bins)
    aves(ii) = mean(xx(ind==ii));
    mm(ii) = mean(yy(ind==ii));
    se(ii) = std(yy(ind==ii))/(numel(yy(ind==ii)));
end
errorbar(aves, mm, se, 'b', 'linewidth', 2);

[rho, p] = corr( xx, yy ); disp(p);
title([num2str(round(rho,4)) '-' num2str(round(p,4))]);
xlabel(xlb); ylabel(ylb);
set(gcf, 'position', [333 222 288 234]);



%% CELL SEQUENCES

for VR = 1:4
    switch VR
        case 1
            load('AoldTwice.mat'); ttl = 'AO1';
        case 2
            load('AnewTwice.mat'); ttl = 'AO2';
        case 3
            load('VnewTwice.mat'); ttl = 'VS2';
        case 4
            load('VoldTwice.mat'); ttl = 'VS1';
    end
    
    df = cat(2, allDfof{:});
    fields = cat(2, fieldsInfo{:});
    
    seq = df.*fields;
    seq1 = seq(:,1:2:size(seq,2));
    seq2 = seq(:,2:2:size(seq,2));

    iidx = find(sum(seq1,1)>0 & sum(seq2,1)>0);
    plotSequenceByMax2( seq1(:,iidx), seq2(:,iidx) ); colormap((parula)); sgtitle(ttl);
%     figure; plotSequenceByMax( seq1(:,iidx), seq2(:,iidx) ); colormap((parula)); title(ttl);
    set(gcf, 'Position', [111*VR 111*VR 456 333]);
end



%% CELL SEQUENCES - for cue cells

for VR = 1:4
    switch VR
        case 1
            load('AoldTwice.mat'); ttl = 'AO1';
        case 2
            load('AnewTwice.mat'); ttl = 'AO2';
        case 3
            load('VnewTwice.mat'); ttl = 'VS2';
        case 4
            load('VoldTwice.mat'); ttl = 'VS1';
    end
    
    df = cat(2, allDfof{:});
    df1 = df(:,1:2:size(df,2));
    df2 = df(:,2:2:size(df,2));
    fields = cat(2, fieldsInfo{:});
    ff1 = fields(:,1:2:size(fields,2));
    ff2 = fields(:,2:2:size(fields,2));
    
    iidx = find(cueScores(:,1) > 0.05);
    
    seq1 = df1(:, iidx) .* ff1(:, iidx);
    seq2 = df2(:, iidx) .* ff2(:, iidx);

    iidx = find(sum(seq1,1)>0 & sum(seq2,1)>0);
    plotSequenceByMax2( seq1(:,iidx), seq2(:,iidx) ); colormap((parula)); sgtitle(ttl);
    set(gcf, 'Position', [111*VR 111*VR 456 333]);
end

