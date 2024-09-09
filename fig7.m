
%load data: contact lead author
%%
allIdx=unique([aoCellIdx;vsCellIdx;aoRewardIdx;vsRewardIdx]);
otherIdx=setdiff([1:1:size(dfofUse,1)],allIdx);
otherIdx=[otherIdx intersect([aoCellIdx;aoRewardIdx],[vsCellIdx;vsRewardIdx])'];
save('otherIdx.mat','otherIdx')
dfofOther=dfofUse(otherIdx,:);
dfofFOther=fieldDistriUse(otherIdx,:);

figure; plotSequenceByMax( dfofOther.*dfofFOther );
% figure; plotSequenceByMax( dfofR);
title('allOthercells');
saveas(gcf,'fieldMatrix_allOtherCells.fig')

%% pie

i1=length(intersect(aoCellIdx,aoRewardIdx));
i2=length(setdiff(aoCellIdx,aoRewardIdx));
i3=length(setdiff(aoRewardIdx,aoCellIdx));

i4=length(intersect(vsCellIdx,vsRewardIdx));
i5=length(setdiff(vsCellIdx,vsRewardIdx));
i6=length(setdiff(vsRewardIdx,vsCellIdx));

i7=length(intersect([aoCellIdx;aoRewardIdx],[vsCellIdx;vsRewardIdx]));
i8=length(intersect(aoCellIdx,vsCellIdx));
i9=length(intersect(aoRewardIdx,vsRewardIdx));
i10=length(intersect(aoRewardIdx,vsCellIdx));
i11=length(intersect(vsRewardIdx,aoCellIdx));

%i10 and i11 are not zero. We include them into the other cell chart

nTotal=3333;%i10+i11=i7
nAO=length(unique([aoCellIdx;aoRewardIdx]))-i7;
nVS=length(unique([vsCellIdx;vsRewardIdx]))-i7;
figure;
ax = gca();
% pie(ax, [(nTotal-nAO-nVS) nAO nVS], [1 1 1], {'other', 'audio', 'visual'}');
pie(ax, [(nTotal-nAO-nVS) nAO nVS], {'other', 'audio', 'visual'}');

ax.Colormap = [...
    0.8 0.8 0.8;
    1 0 1;
    0 1 0];
set(gcf, 'position', [222 222 333 333]);

saveas(gcf,'pie.fig');

nAO/nTotal
nVS/nTotal
1-nAO/nTotal-nVS/nTotal

%% field distribution

i=find(nRuns(:,1)>1);
% allDfofsUse=allDfofs(i);
% %similar to the previousone, need to interplate allDfof
% for n=1:length(allDfofsUse);
%     for m=1:size(allDfofsUse{n},1);
%         allDfofsUse{n}(m,:)=naninterp(allDfofsUse{n}(m,:));
%     end
% end
allFieldsUse=allFields(i);

idxOtherFOV={};
for n=1:length(allFieldsUse);
    AVcellIdx=unique([idxAOFOV{n};idxVSFOV{n};idxARFOV{n};idxVRFOV{n}]);
    idxOtherFOV{n}=setdiff([1:1:size(allFieldsUse{n},1)],AVcellIdx);
end

fieldsOther = cellfun(@(x,y) x(y,:), allFieldsUse, idxOtherFOV, 'uni', 0);
fields = cellfun(@(x) nanmean(x,1), fieldsOther, 'uni', 0);
fields = cat(1, fields{:});

%shuffle fields
nShuffle=1000;
shuffleFields={};
for n=1:length(fieldsOther);
     shuffleFields{n}=[];
    A=fieldsOther{n};
    for i=1:nShuffle;
            S=[];
        for m=1:size(A,1);
            aa=A(m,:);
            S(m,:)=aa(randperm(length(aa)));
        end
     shuffleFields{n}=[shuffleFields{n};S];
    end
  
end

shuffles = cellfun(@(x) nanmean(x,1), shuffleFields, 'uni', 0);
shuffles = cat(1, shuffles{:});


xidx = 1:140;
templateAO = [29:38 43:52 63:72 111:120];
templateVS = [11:16 83:88 97:102 129:134];
rew = [39:42 91:94];

figure; hold on
[tempmean, ~, seUpper, seLower] = getMeanAndSE(fields, 1);
plot(xidx, tempmean(xidx), 'k', 'LineWidth', 2);
patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'k', 'FaceAlpha', 0.2, 'LineStyle', 'none');

[tempmean, ~, seUpper, seLower] = getMeanAndSE(shuffles, 1);
plot(xidx, tempmean(xidx), 'r', 'LineWidth', 2);
patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'r', 'FaceAlpha', 0.2, 'LineStyle', 'none');


slideWindow = 2;
for ii = 1:(140-slideWindow)
    aa = fields(:,ii:(ii+slideWindow)); aa = reshape(aa, 1, []);
    bb = shuffles(:,ii:(ii+slideWindow)); bb = reshape(bb, 1, []);
%     [hh,~] = ttest2( aa, bb, 'Tail', 'right' );
%     if hh
%         text(ii+slideWindow/2, 0.5, '*', 'Color', 'k');
%     end
%     [hh,~] = ttest2( aa, bb, 'Tail', 'left' );
%     if hh
%         text(ii+slideWindow/2, 0.5, '*', 'Color', 'r');
%     end
    [hh,~] = ttest2( aa, bb );
    if hh
        text(ii+slideWindow/2, 0.5, '*', 'Color', 'k');
    end
end

yl = get(gca, 'ylim');
bar(templateAO-0.5, ones(1, length(templateAO))*yl(2), 'BaseValue', yl(1), 'EdgeColor', 'none', 'BarWidth', 1, 'FaceColor', 'm', 'FaceAlpha', 0.1);
bar(templateVS-0.5, ones(1, length(templateVS))*yl(2), 'BaseValue', yl(1), 'EdgeColor', 'none', 'BarWidth', 1, 'FaceColor', 'g', 'FaceAlpha', 0.1);
bar(rew-0.5, ones(1, length(rew))*yl(2), 'BaseValue', yl(1), 'EdgeColor', 'none', 'BarWidth', 1, 'FaceColor', 'y', 'FaceAlpha', 0.3);
xlim([0 140]+0.5);
set(gca, 'xtick', ([0 200 450 700]+2.5)/5, 'xticklabel', {'0', 'AO-r', 'VS-r', '700'});
set(gcf, 'Position', [222 222 555 288]);
title('FIELDS'); xlabel('Distance along track (cm)'); ylabel('% of cells with fields');

saveas(gcf,'fieldDistri.fig')

%% did not use this


fovIndices = find(nRuns(:,1)>1);

fields = cellfun(@(x,y) mean(x(y,:),1), allFieldsUse,idxOtherFOV, 'uni', 0);
% fields = cat(1, fields{:});
% fields=fields(fovIndices);

figure; hold on

% ao
% iidx = [31:36 45:50 65:70 113:118]; xidx = 1;
iidx = [29:38 39:42 43:52 63:72 111:120]; xidx = 1; %add reward
xx = cellfun(@(x) nanmean(x(iidx),'all'), fields, 'uni', 0);
xx = cell2mat(xx);
% xx = fields(:, iidx);
bar(xidx, nanmean(xx), 'm', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(xx, 'all'), nanstd(xx, [], 'all')/sqrt(size(xx,1)), 'Color', 'm', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'm', 'LineWidth', 2, 'LineStyle', 'none');

% vs
% iidx = [11:16 83:88 97:102 129:134]; xidx = 2;
iidx = [9:18 81:90 91:94 95:104 127:136]; xidx = 2;
yy = cellfun(@(x) nanmean(x(iidx),'all'), fields, 'uni', 0);
yy = cell2mat(yy);
% yy = fields(:, iidx);
bar(xidx, nanmean(yy), 'g', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(yy, 'all'), nanstd(yy, [], 'all')/sqrt(size(yy,1)), 'Color', 'g', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'g', 'LineWidth', 2, 'LineStyle', 'none');
    
xlim([0.5 2.5]); 
set(gca, 'xtick', 1:2, 'xticklabel', {'AO', 'VS'});
set(gcf, 'Position', [123 123 234 222]);

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);
title(['Mean % cells',num2str(p)]); ylabel('% cells');
saveas(gcf,'percentageOfCellsWithFields.fig')

%% percent of bins iwth field: in VS and AO zones
aodx = [29:38 39:42 43:52 63:72 111:120];
vsdx = [9:18 81:90 91:94 95:104 127:136];
F=allFieldsUse;
F=cell2mat(F');
F=F(otherIdx,:);

AOVSF=[];

for n=1:size(F,1);
    a=F(n,:);
    a1=a(aodx);
    a2=a(vsdx);
    AOVSF(n,1)=nansum(a1)/length(a1);
    AOVSF(n,2)=nansum(a2)/length(a2);
end

xx=AOVSF(:,1); xidx=1;

figure
hold on
bar(xidx, nanmean(xx), 'm', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(xx, 'all'), nanstd(xx, [], 'all')/sqrt(size(xx,1)), 'Color', 'm', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'm', 'LineWidth', 2, 'LineStyle', 'none');

yy=AOVSF(:,2); xidx=2;

bar(xidx, nanmean(yy), 'g', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(yy, 'all'), nanstd(yy, [], 'all')/sqrt(size(yy,1)), 'Color', 'g', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'g', 'LineWidth', 2, 'LineStyle', 'none');

xlim([0.5 2.5]); 
set(gca, 'xtick', 1:2, 'xticklabel', {'AO', 'VS'});
set(gcf, 'Position', [123 123 234 222]);

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);
title(['field coverage',num2str(p)]); ylabel('% cells');
saveas(gcf,'fieldCoverage.fig')

% figure
% for n=1:size(AOVSF,1);
%     hold on
%     plot(AOVSF(n,:),'color',[0.8 0.8 0.8]);
% end

%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate reward consistency on my own

%do not re run
p=pwd;
load('folders_combo.mat')
rewardConsistency={};
ra=[39:42];
rv=[91:94];
for n=1:length(folders_combo)
    disp(n)
    cd(folders_combo{n})
    load('roiIdxUse.mat')
    load('RunByRun_sig\dfofMInterpM_sig.mat');
    load('successfulRuns_multipleRewards.mat');
    %runs that have two rewards
    irun=successfulRuns_multipleRewards{1}.*successfulRuns_multipleRewards{2};
    irun=find(irun); %only use runs with two rewards
    dfofMInterpM_sigUse=dfofMInterpM_sig(find(roiIdxUse));
    c=[];
    if irun<=1;
        c=nan(length(dfofMInterpM_sigUse),2);
    else
    for m=1:length(dfofMInterpM_sigUse);
        ram=dfofMInterpM_sigUse{m}(irun,ra);
        rvm=dfofMInterpM_sigUse{m}(irun,rv);

        cram=[];
        crvm=[];
        for i=1:size(ram,1)-1;
            for ii=i+1:size(ram,1);
                cram(end+1)=corr(ram(i,:)',ram(ii,:)');
                crvm(end+1)=corr(rvm(i,:)',rvm(ii,:)');
            end
        end
        c(m,1)=nanmean(cram);
        c(m,2)=nanmean(crvm);
    end
    rewardConsistency{n}=c;
    end
end
cd(p);
save('rewardConsistency.mat','rewardConsistency')

%% recalculate cue consistency on my own because somehow Duc's data showed different numbers (3332 cells instead of 3333 cells).

%do not re run
pp=pwd;
load('folders_combo.mat')
cueConsistenciesYG={};
aoIdx = [29:38;43:52;63:72;111:120];
vsIdx = [9:18;81:90;95:104;127:136];
aovsIdx=[aoIdx;vsIdx];
for n=1:length(folders_combo)

    disp(n)
    cd(folders_combo{n})
    load('roiIdxUse.mat')
    load('RunByRun_sig\dfofMInterpM_sig.mat');
  load('successfulRuns_multipleRewards.mat');
    %runs that have two rewards
    irun=successfulRuns_multipleRewards{1}.*successfulRuns_multipleRewards{2};
    irun=find(irun);%only use runs with two rewards
    dfofMInterpM_sigUse=dfofMInterpM_sig(find(roiIdxUse));
    c=[];
    if irun<=1;
        c=nan(length(dfofMInterpM_sigUse),140);
    else
  
    for m=1:length(dfofMInterpM_sigUse);
cThisCell=nan(1,140);

        for p=1:size(aovsIdx,1);
            this=dfofMInterpM_sigUse{m}(irun,aovsIdx(p,:));
            cThisCue=[];
        for i=1:size(this,1)-1;
            for ii=i+1:size(this,1);
                cThisCue(end+1)=corr(this(i,:)',this(ii,:)');
            end
        end

        
       cThisCue=nanmean(cThisCue);
       cThisCell(1,aovsIdx(p,:))=cThisCue;
        end
         c=[c;cThisCell];
    end
cueConsistenciesYG{n}=c;
    end
end
cd(pp);
save('cueConsistenciesYG.mat','cueConsistenciesYG')
   

%% dfof along the track
dfofOther=dfofUse(otherIdx,:);

shuffles=[];
for m=1:size(dfofOther,1);
     this=dfofOther(m,:);
    s=[];
for n=1:1000;     
        s(n,:)=this(randperm(length(this)));
    end
    shuffles(m,:)=nanmean(s,1);
    % shuffles=[shuffles;s];
end


xidx = 1:140;
templateAO = [29:38 43:52 63:72 111:120];
templateVS = [11:16 83:88 97:102 129:134];
rew = [39:42 91:94];



figure; hold on
[tempmean, ~, seUpper, seLower] = getMeanAndSE(dfofOther, 1);
plot(xidx, tempmean(xidx), 'k', 'LineWidth', 2);
patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'k', 'FaceAlpha', 0.2, 'LineStyle', 'none');

[tempmean, ~, seUpper, seLower] = getMeanAndSE(shuffles, 1);
plot(xidx, tempmean(xidx), 'r', 'LineWidth', 2);
patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'r', 'FaceAlpha', 0.2, 'LineStyle', 'none');


slideWindow = 2;
for ii = 1:(140-slideWindow)
    aa = dfofOther(:,ii:(ii+slideWindow)); aa = reshape(aa, 1, []);
    bb = shuffles(:,ii:(ii+slideWindow)); bb = reshape(bb, 1, []);
%     [hh,~] = ttest2( aa, bb, 'Tail', 'right' );
%     if hh
%         text(ii+slideWindow/2, 0.5, '*', 'Color', 'k');
%     end
%     [hh,~] = ttest2( aa, bb, 'Tail', 'left' );
%     if hh
%         text(ii+slideWindow/2, 0.5, '*', 'Color', 'r');
%     end
    [hh,~] = ttest2( aa, bb );
    if hh
        text(ii+slideWindow/2, 0.1, '*', 'Color', 'k');
    end
end

yl = get(gca, 'ylim');
bar(templateAO-0.5, ones(1, length(templateAO))*yl(2), 'BaseValue', yl(1), 'EdgeColor', 'none', 'BarWidth', 1, 'FaceColor', 'm', 'FaceAlpha', 0.1);
bar(templateVS-0.5, ones(1, length(templateVS))*yl(2), 'BaseValue', yl(1), 'EdgeColor', 'none', 'BarWidth', 1, 'FaceColor', 'g', 'FaceAlpha', 0.1);
bar(rew-0.5, ones(1, length(rew))*yl(2), 'BaseValue', yl(1), 'EdgeColor', 'none', 'BarWidth', 1, 'FaceColor', 'y', 'FaceAlpha', 0.3);
xlim([0 140]+0.5);
set(gca, 'xtick', ([0 200 450 700]+2.5)/5, 'xticklabel', {'0', 'AO-r', 'VS-r', '700'});
set(gcf, 'Position', [222 222 555 288]);
title('DFOF'); xlabel('Distance along track (cm)'); ylabel('% of cells with fields');


saveas(gcf,'dfofdistribution.fig');

%% recalculate rbr consistency on my own because somehow Duc's data showed different numbers (3332 cells instead of 3333 cells).

% %do not re run
% pp=pwd;
% load('folders_combo.mat')
% rbrConsistenciesYG={};
% % aoIdx = [29:38;43:52;63:72;111:120];
% % vsIdx = [9:18;81:90;95:104;127:136];
% % aovsIdx=[aoIdx;vsIdx];
% for n=1:length(folders_combo)
% 
%     disp(n)
%     cd(folders_combo{n})
%     load('roiIdxUse.mat')
%     load('RunByRun_sig\corrInfo.mat');
%   rbrConsistenciesYG{n}=corrInfo.meantoOthers(find(roiIdxUse));
% 
% end
% cd(pp);
% save('rbrConsistenciesYG.mat','rbrConsistenciesYG')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT - cue consistency and rewrd consistency
% panel 7F
load('rewardConsistency.mat')
fovIndices = find(nRuns(:,1)>1);

rewardConsistencyUse=rewardConsistency(fovIndices);
rewardConsistencyUse=cell2mat(rewardConsistencyUse');

%remove the 4th fov after this because the number of cells do not match
figure; hold on
% cuesCorr = cat(1, cueConsistencies{1:2,fovIndices});
% cuesCorr = cellfun(@(x) x, cueConsistencies(useR,:), 'uni', 0);

% ao
iidx = [34 48 68 116];
% xidx = [2 3 4 7];
xidx = [2 3 4 5 9];
% xx = cat(1, cueConsistencies{2,fovIndices});
xx=cell2mat(cueConsistenciesYG(fovIndices)');
xx = xx(:, iidx);
%include reward
xxr=[];
xxr(:,1)=xx(:,1);
xxr(:,2)=rewardConsistencyUse(:,1);
xxr(:,3:5)=xx(:,2:4);
xxr=xxr(otherIdx,:);

bar(xidx, nanmean(xxr, 1), 'm', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(xxr, 1), nanstd(xxr, 1)/sqrt(size(xxr,1)), 'Color', 'm', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'm', 'LineWidth', 2, 'LineStyle', 'none');
    
% vs
iidx = [14 86 100 132];
xidx = [1 6 7 8 10];
yy=cell2mat(cueConsistenciesYG(fovIndices)');
yy = yy(:, iidx);

yyr=[];
yyr(:,1)=yy(:,1);
yyr(:,2)=rewardConsistencyUse(:,2);
yyr(:,3:5)=yy(:,2:4);
yyr=yyr(otherIdx,:);
bar(xidx, nanmean(yyr, 1), 'g', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(yyr, 1), nanstd(yyr, 1)/sqrt(size(yyr,1)), 'Color', 'g', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'g', 'LineWidth', 2, 'LineStyle', 'none');
    
% zz = cuesCorr(:, [34 48 68 116 14 86 100 132]);
zz=[xxr yyr];
zz = nanmean(zz, 2);
mm = nanmean(zz);
se = nanstd(zz) / sqrt(numel(zz));
plot([0.5 10.5], [mm mm], 'k');
plot([0.5 10.5], [mm+se mm+se], 'k--');
plot([0.5 10.5], [mm-se mm-se], 'k--');

[~,p] = ttest(xxr(:,1), mm); disp(p);
[~,p] = ttest(xxr(:,2), mm); disp(p);
[~,p] = ttest(xxr(:,3), mm); disp(p);
[~,p] = ttest(xxr(:,4), mm); disp(p);
[~,p] = ttest(xxr(:,5), mm); disp(p);

[~,p] = ttest(yyr(:,1), mm); disp(p);
[~,p] = ttest(yyr(:,2), mm); disp(p);
[~,p] = ttest(yyr(:,3), mm); disp(p);
[~,p] = ttest(yyr(:,4), mm); disp(p);
[~,p] = ttest(yyr(:,5), mm); disp(p);

% yl = get(gca, 'ylim');
% bar([2 3 4 7], ones(1, 4)*yl(2), 'BaseValue', yl(1), 'EdgeColor', 'none', 'BarWidth', 0.8, 'FaceColor', 'm', 'FaceAlpha', 0.1);
% bar([1 5 6 8], ones(1, 4)*yl(2), 'BaseValue', yl(1), 'EdgeColor', 'none', 'BarWidth', 0.8, 'FaceColor', 'g', 'FaceAlpha', 0.1);

xlim([0.5 10.5]); xlabel('Cue Order along track');
set(gca, 'xtick', [2.5 5.5], 'xticklabel', {'AO-r', 'VS-r'});
set(gcf, 'Position', [123 123 333 222]);

saveas(gcf,'cueRewardConsistency.fig')

%% 
%% run-consistency in AUDIO vs VISUAL
% panel 7G

fovIndices = find(nRuns(:,1)>1);

figure; hold on

% ao
iidx = [34 48 68 116]; xidx = 1;
xx=cell2mat(cueConsistenciesYG(fovIndices)');
xx = xx(:, iidx);
%include reward
xxr=[];
xxr(:,1)=xx(:,1);
xxr(:,2)=rewardConsistencyUse(:,1);
xxr(:,3:5)=xx(:,2:4);
xxr=xxr(otherIdx,:);
bar(xidx, nanmean(xxr, 'all'), 'm', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(xxr, 'all'), nanstd(xxr, [], 'all')/sqrt(size(xxr,1)), 'Color', 'm', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'm', 'LineWidth', 2, 'LineStyle', 'none');

% vs
iidx = [14 86 100 132]; xidx = 2;
yy=cell2mat(cueConsistenciesYG(fovIndices)');
yy = yy(:, iidx);

yyr=[];
yyr(:,1)=yy(:,1);
yyr(:,2)=rewardConsistencyUse(:,2);
yyr(:,3:5)=yy(:,2:4);
yyr=yyr(otherIdx,:);
bar(xidx, nanmean(yyr, 'all'), 'g', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(yyr, 'all'), nanstd(yyr, [], 'all')/sqrt(size(yyr,1)), 'Color', 'g', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'g', 'LineWidth', 2, 'LineStyle', 'none');
    
xlim([0.5 2.5]); 
set(gca, 'xtick', 1:2, 'xticklabel', {'AO', 'VS'});
set(gcf, 'Position', [123 123 222 222]);

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);
title(['Runs Consistency',num2str(p)]);
saveas(gcf,'RBRConsistencyAtAOVS.fig')

%% Mean-DFOF in AUDIO vs VISUAL

figure; hold on

% ao
iidx = [29:38 39:42 43:52 63:72 111:120]; xidx = 1;
xx = dfofUse; xx = xx(otherIdx, iidx);
bar(xidx, nanmean(xx, 'all'), 'm', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(xx, 'all'), nanstd(xx, [], 'all')/sqrt(size(xx,1)), 'Color', 'm', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'm', 'LineWidth', 2, 'LineStyle', 'none');

% vs
iidx = [9:18 81:90 91:94 95:104 127:136]; xidx = 2;
yy = dfofUse; yy = yy(otherIdx, iidx);
bar(xidx, nanmean(yy, 'all'), 'g', 'FaceAlpha', 0.5);
errorbar(xidx, nanmean(yy, 'all'), nanstd(yy, [], 'all')/sqrt(size(yy,1)), 'Color', 'g', ...
        'MarkerSize', 20, 'MarkerFaceColor', 'g', 'LineWidth', 2, 'LineStyle', 'none');
    
xlim([0.5 2.5]); title('Mean \DeltaF/F'); ylabel('\DeltaF/F');
set(gca, 'xtick', 1:2, 'xticklabel', {'AO', 'VS'});
set(gcf, 'Position', [123 123 234 222]);

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);
title(['dfof',num2str(p)])

saveas(gcf,'dfofAOVS.fig')
%% calculate the different responses to auditory and visual cues on cell basis
iidx = [29:38 39:42 43:52 63:72 111:120]; xidx = 1;%ao response
xx = dfofUse; xx = xx(otherIdx, iidx);
x=nanmean(xx,2);
iidx = [9:18 81:90 91:94 95:104 127:136]; xidx = 2;%vs response
yy = dfofUse; yy = yy(otherIdx, iidx);
y=nanmean(yy,2);

figure,
subplot(222)
[x1,p1]=ksdensity(a);
plot(x1,p1)
hold on
line([0 max(x1)],[0 0])
ylim([min(p1) max(p1)])

subplot(221)
a=x-y;
a=a(randperm(length(a)));
plot(a,'.','Color',[0.8 0.8 0.8]);
hold on
line([0 length(a)],[0 0])
ylim([min(p1) max(p1)])

subplot(223)
plot(x,y,'k.')
line([0 0.5],[0 0.5])

saveas(gcf,'sensorydFOFDIFFRENCE.fig')

