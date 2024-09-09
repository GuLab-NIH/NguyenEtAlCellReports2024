%load data: contact lead author
%%
fields=dfofUse;
% load('fieldDistriUse.mat');
% dfofRF=fieldDistriUse;
% fields=dfofUse.*dfofRF;
d=[];%dfof within cues before the reward compared to other cues
r=[39:42 91:94];
nr=setdiff([1:1:140],r);
rdfof = fields(:, r); rdfof = nanmean(rdfof, 2);
nrdfof = fields(:, nr); nrdfof = nanmean(nrdfof, 2);
d = (rdfof - nrdfof) ./ (rdfof + nrdfof);

%%
load('shF.mat');
%%
allF=cell2mat(shF');
dShuffle=[];
rdfof = allF(:, r); rdfof = nanmean(rdfof, 2);
nrdfof = allF(:, nr); nrdfof = nanmean(nrdfof, 2);
dShuffle = (rdfof - nrdfof) ./ (rdfof + nrdfof);
dShuffle=dShuffle(abs(dShuffle)<=1);
threshReward=prctile(dShuffle,99);
thresh1=prctile(dShuffle,1);

idx=find(d>threshReward);

dfofR=fields(idx,:);
[p1,x1]=ksdensity(dShuffle,'width',0.05);
[p2,x2]=ksdensity(d,'width',0.05);
figure,
plot(x1,p1/max(p1))
hold on
plot(x2,p2/max(p2));
hold on
line([threshReward threshReward],[0 1])
x3=x2(x2>=threshReward);
p3=p2(x2>=threshReward)/max(p2);
hold on
plot([threshReward threshReward x3],[0 p3(1)+0.0029524 p3])
saveas(gcf,'rewardScoreDistribution.fig')
save('threshReward.mat','threshReward')

%% plot cells based on fields
load('fieldDistriUse.mat');
dfofRF=fieldDistriUse(idx,:);

figure; plotSequenceByMax( dfofR.*dfofRF );
% figure; plotSequenceByMax( dfofR);
title('all reward cells');

saveas(gcf,'allRewardFields.fig')

%% see whether dfofR has sensory bias
ra=[39:42];
rv=[91:94];

%calculate their activity preferneces

radfof = dfofR(:, ra); radfof = nanmean(radfof, 2);
rvdfof = dfofR(:, rv); rvdfof = nanmean(rvdfof, 2);
dr = (radfof - rvdfof) ./ (radfof + rvdfof);

figure,

subplot(242)
[p,x]=ksdensity(dr,'width',0.1);%show bimodal distribution
plot(x,p/max(p)*95)
hold on
histogram(dr,10)
% hold on
% histogram(dr)
xlim([min(x) max(x)])


subplot(241)
plot(dr,'.')
ylim([min(x) max(x)])

ia=find(dr>0);
iv=find(dr<0);

aoRewardIdx=ia;
vsRewardIdx=iv;
save('aoRewardIdx.mat','aoRewardIdx');
save('vsRewardIdx.mat','vsRewardIdx');
dfofra=dfofR(ia,:);
dfofrv=dfofR(iv,:);
dfofraF=dfofRF(ia,:);
dfofrvF=dfofRF(iv,:);

subplot(243)
imagesc(dfofra);
title('auditory r');
subplot(244)
imagesc(dfofrv);
title('visual r');
subplot(245)
 plotSequenceByMax( dfofra.*dfofraF );
title('auditory r');
subplot(246)
 plotSequenceByMax( dfofrv.*dfofrvF );
title('visual r');

subplot(247)
plot(mean(dfofra,1))
hold on
plot(mean(dfofrv,1))

saveas(gcf,'rewardCells2_inReward.fig')
%%
% %separate reward cells in two cues
% figure,
% subplot(242)
% [p,x]=ksdensity(dr,'width',0.1);%show bimodal distribution
% plot(x,p/max(p)*95)
% hold on
% histogram(dr,10)
% % hold on
% % histogram(dr)
% xlim([min(x) max(x)])
% 
% [a,b]=findpeaks(-p);
% thresh1=x(b(1));
% thresh2=x(b(3));
% 
% subplot(243)
% ia=find(dr>thresh2);
% nAO=length(ia);
% iv=find(dr<thresh1);
% nVS=length(iv);
% iav=find(dr>thresh1&dr<thresh2);
% nTotal=length(dr);
% ax = gca();
% % pie(ax, [(nTotal-nAO-nVS) nAO nVS], [1 1 1], {'other', 'audio', 'visual'}');
% pie(ax, [(nTotal-nAO-nVS) nAO nVS], {'both', 'audio', 'visual'}');
% 
% ax.Colormap = [...
%     0.8 0.8 0.8;
%     1 0 1;
%     0 1 0];
% set(gcf, 'position', [222 222 333 333]);
% 
% 
% dfofra=dfofR(ia,:);
% dfofrv=dfofR(iv,:);
% dfofrav=dfofR(iav,:);
% 
% dfofraF=dfofRF(ia,:);
% dfofrvF=dfofRF(iv,:);
% dfofravF=dfofRF(iav,:);
% 
% subplot(241)
% plot(dr,'.')
% ylim([min(x) max(x)])
% 
% subplot(244)
%  plotSequenceByMax( dfofra.*dfofraF );
% title('auditory r');
% subplot(245)
%  plotSequenceByMax( dfofrv.*dfofrvF );
% title('visual r');
% 
% subplot(246)
%  plotSequenceByMax( dfofrav.*dfofravF );
% title('audiovisual r');



%% reward cell indices per fov
cd ..\
load('fig7.mat')

cd("fig_YG_DFOF_completeRandom_99\")
%% 

load('threshReward.mat')

i=find(nRuns(:,1)>1);
allDfofsUse=allDfofs(i);
%similar to the previousone, need to interplate allDfof
for n=1:length(allDfofsUse);
    for m=1:size(allDfofsUse{n},1);
        allDfofsUse{n}(m,:)=naninterp(allDfofsUse{n}(m,:));
    end
end
allFieldsUse=allFields(i);

r=[39:42 91:94];
nr=setdiff([1:1:140],r);

dd = cellfun(@(x) (nanmean(x(:,r),2) - nanmean(x(:,nr),2)) ...
    ./ (nanmean(x(:,r),2) + nanmean(x(:,nr),2)), allDfofsUse, 'uni', 0);

ra=[39:42];
rv=[91:94];

ddSensory = cellfun(@(x) (nanmean(x(:,ra),2) - nanmean(x(:,rv),2)) ...
    ./ (nanmean(x(:,ra),2) + nanmean(x(:,rv),2)), allDfofsUse, 'uni', 0);


idxARFOV={};
idxVRFOV={};
for n=1:length(dd);
  idxARFOV{n}=find(dd{n}>threshReward&ddSensory{n}>0);
idxVRFOV{n}=find(dd{n}>threshReward&ddSensory{n}<0);

end

save('idxARFOV.mat','idxARFOV');
save('idxVRFOV.mat','idxVRFOV')

ddAO = cellfun(@(x, y,z) x(y > threshReward&z>0,:), allFieldsUse, dd, ddSensory,'uni', 0);
ddAO = cellfun(@(x) nanmean(x,1), ddAO, 'uni', 0);
ddAO = cat(1, ddAO{:});

ddVS = cellfun(@(x, y,z) x(y > threshReward&z<0,:), allFieldsUse, dd, ddSensory,'uni', 0);
ddVS = cellfun(@(x) nanmean(x,1), ddVS, 'uni', 0);
ddVS = cat(1, ddVS{:});
% 

xidx = 1:140;
templateAO = [29:38 43:52 63:72 111:120];
templateVS = [11:16 83:88 97:102 129:134];
rew = [39:42 91:94];


% COMPARE BOTH
figure; hold on
[tempmean, ~, seUpper, seLower] = getMeanAndSE(ddAO, 1);
plot(xidx, tempmean(xidx), 'm', 'LineWidth', 2);
patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'm', 'FaceAlpha', 0.2, 'LineStyle', 'none');

[tempmean, ~, seUpper, seLower] = getMeanAndSE(ddVS, 1);
plot(xidx, tempmean(xidx), 'g', 'LineWidth', 2);
patch([xidx flip(xidx)], [seLower(xidx) flip(seUpper(xidx))], 'g', 'FaceAlpha', 0.2, 'LineStyle', 'none');

slideWindow = 0;
for ii = 1:(140-slideWindow)
    aa = ddAO(:,ii:(ii+slideWindow)); aa = reshape(aa, 1, []);
    bb = ddVS(:,ii:(ii+slideWindow)); bb = reshape(bb, 1, []);
%     [hh,~] = ttest2( aa, bb, 'Tail', 'right' );
%     if hh
%         text(ii+slideWindow/2, 0.5, '*', 'Color', 'm');
%     end
%     [hh,~] = ttest2( aa, bb, 'Tail', 'left' );
%     if hh
%         text(ii+slideWindow/2, 0.5, '*', 'Color', 'g');
%     end
    [hh,~] = ttest2( aa, bb );
    if hh
        text(ii+slideWindow/2, 0.6, '*', 'Color', 'k');
    end
end

yl = get(gca, 'ylim');
bar(templateAO-0.5, ones(1, length(templateAO))*yl(2), 'BaseValue', yl(1), 'EdgeColor', 'none', 'BarWidth', 1, 'FaceColor', 'm', 'FaceAlpha', 0.1);
bar(templateVS-0.5, ones(1, length(templateVS))*yl(2), 'BaseValue', yl(1), 'EdgeColor', 'none', 'BarWidth', 1, 'FaceColor', 'g', 'FaceAlpha', 0.1);
bar(rew-0.5, ones(1, length(rew))*yl(2), 'BaseValue', yl(1), 'EdgeColor', 'none', 'BarWidth', 1, 'FaceColor', 'y', 'FaceAlpha', 0.3);
xlim([0 140]+0.5); ylim(yl);
set(gca, 'xtick', ([0 200 450 700]+2.5)/5, 'xticklabel', {'0', 'AO-r', 'VS-r', '700'});
set(gcf, 'Position', [222 222 555 288]);
title('FIELDS'); xlabel('Distance along track (cm)'); ylabel('% of cells with fields');

saveas(gcf,'fieldCoverage_reward.fig');

%%
%% rbr consistencies

aoCellIdx=aoRewardIdx;
vsCellIdx=vsRewardIdx;
i=find(nRuns(:,1)>1);

rbrConsistenciesUse=rbrConsistencies(i);

ddAO = cat(1, rbrConsistenciesUse{:});
% ddAO=ddAO(iUse);
ddAO = ddAO(aoCellIdx);

ddVS = cat(1, rbrConsistenciesUse{:});
% ddVS=ddVS(iUse);
ddVS = ddVS(vsCellIdx);

figure; hold on

% ao
xidx = 1; clr = 'm';
[tempmean, se, ~,~] = getMeanAndSE(ddAO, 'all');
bar( xidx, tempmean, 'FaceColor', clr, 'FaceAlpha', 0.5);
errorbar( xidx, tempmean, se, clr);

% ao
xidx = 2; clr = 'g';
[tempmean, se, ~,~] = getMeanAndSE(ddVS, 'all');
bar( xidx, tempmean, 'FaceColor', clr, 'FaceAlpha', 0.5);
errorbar( xidx, tempmean, se, clr);

[~,p] = ttest2(ddAO, ddVS); disp(p);

set(gca, 'xtick', 1:2, 'xticklabel', {'Audio indices', 'Visual indices'}, 'xticklabelrotation', 30);
set(gcf, 'Position', [222 222 234 333]);

title(['rbr consistency, p=',num2str(p)])
xlim([0.5 2.5])

saveas(gcf,'rbrconsistencyReward.fig');


%% MEAN DFOF

allAverageDfofsUse=allAverageDfofs(i);


ddAO = cat(2, allAverageDfofsUse{:});
% ddAO=ddAO(iUse);
ddAO = ddAO(aoCellIdx);

ddVS = cat(2, allAverageDfofsUse{:});
% ddVS=ddVS(iUse);
ddVS = ddVS(vsCellIdx);


figure; hold on

% ao
xidx = 1; clr = 'm';
[tempmean, se, ~,~] = getMeanAndSE(ddAO, 'all');
bar( xidx, tempmean, 'FaceColor', clr, 'FaceAlpha', 0.5);
errorbar( xidx, tempmean, se, clr);

% ao
xidx = 2; clr = 'g';
[tempmean, se, ~,~] = getMeanAndSE(ddVS, 'all');
bar( xidx, tempmean, 'FaceColor', clr, 'FaceAlpha', 0.5);
errorbar( xidx, tempmean, se, clr);

[~,p] = ttest2(ddAO, ddVS); disp(p);

set(gca, 'xtick', 1:2, 'xticklabel', {'Audio indices', 'Visual indices'}, 'xticklabelrotation', 30);
set(gcf, 'Position', [222 222 234 333]);

title(['mean dfof, p=',num2str(p)])
xlim([0.5 2.5])

saveas(gcf,'meandfofReward.fig');


% %% PLOT - cue consistencies direct comparison
% cueConsistenciesUse=cueConsistencies(:,i);
% % 10-bin
% 
% i=find(nRuns(:,1)>1);
% allDfofsUse=allDfofs(i);
% 
% aoArea = [29:38 43:52 63:72 111:120];
% vsArea = [9:18 81:90 95:104 127:136];
% 
% dd = cellfun(@(x) (mean(x(:,aoArea),2) - mean(x(:,vsArea),2)) ...
%     ./ (mean(x(:,aoArea),2) + mean(x(:,vsArea),2)), allDfofsUse, 'uni', 0);
% 
% ddAO = cellfun(@(x, y) x(y > thresh95,:), cueConsistenciesUse([1 3],:), [dd; dd], 'uni', 0);
% ddAO = cat(1, ddAO{:});
% ddAO = ddAO(:, [34 48 68 116 14 86 100 132]);
% 
% ddVS = cellfun(@(x, y) x(y < thresh5,:), cueConsistenciesUse([1 3],:), [dd; dd], 'uni', 0);
% ddVS = cat(1, ddVS{:});
% ddVS = ddVS(:, [34 48 68 116 14 86 100 132]);
% 
% 
% figure; hold on
% 
% % ao
% xidx = [2 3 4 7]-.2; xx = ddAO(:, 1:4);
% bar(xidx, nanmean(xx, 1), 'm', 'FaceAlpha', 0.7, 'Barwidth', 0.4);
% errorbar(xidx, nanmean(xx, 1), nanstd(xx, 1)/sqrt(size(xx,1)), 'Color', 'm', ...
%         'MarkerSize', 20, 'MarkerFaceColor', 'm', 'LineWidth', 2, 'LineStyle', 'none');
% 
% % ao
% xidx = [2 3 4 7]+.2; xx = ddVS(:, 1:4);
% bar(xidx, nanmean(xx, 1), 'g', 'FaceAlpha', 0.7, 'Barwidth', 0.4);
% errorbar(xidx, nanmean(xx, 1), nanstd(xx, 1)/sqrt(size(xx,1)), 'Color', 'g', ...
%         'MarkerSize', 20, 'MarkerFaceColor', 'g', 'LineWidth', 2, 'LineStyle', 'none');
% 
% % vs
% xidx = [1 5 6 8]-.2; yy = ddAO(:, 5:8);
% bar(xidx, nanmean(yy, 1), 'm', 'FaceAlpha', 0.7, 'Barwidth', 0.4);
% errorbar(xidx, nanmean(yy, 1), nanstd(yy, 1)/sqrt(size(yy,1)), 'Color', 'm', ...
%         'MarkerSize', 20, 'MarkerFaceColor', 'm', 'LineWidth', 2, 'LineStyle', 'none');
% 
% % vs
% xidx = [1 5 6 8]+.2; yy = ddVS(:, 5:8);
% bar(xidx, nanmean(yy, 1), 'g', 'FaceAlpha', 0.7, 'Barwidth', 0.4);
% errorbar(xidx, nanmean(yy, 1), nanstd(yy, 1)/sqrt(size(yy,1)), 'Color', 'g', ...
%         'MarkerSize', 20, 'MarkerFaceColor', 'g', 'LineWidth', 2, 'LineStyle', 'none');
% 
% yl = get(gca, 'ylim');
% xidx = [2 3 4 7 1 5 6 8];
% phh=[];
% for ii = 1:8
%     if ii < 4.5
%         [hh,phh(ii)] = ttest2( ddAO(:, ii), ddVS(:, ii), 'tail', 'right' );
%     else
%         [hh,phh(ii)] = ttest2( ddAO(:, ii), ddVS(:, ii), 'tail', 'left' );
%     end
%     if hh
%         text(xidx(ii), yl(2), '*', 'horizontalAlignment', 'center');
%     else
%         text(xidx(ii), yl(2), 'n.s.', 'horizontalAlignment', 'center');
%     end
% end
% xlabel('Cue Order along track');
% set(gca, 'xtick', [2.5 5.5], 'xticklabel', {'AO-r', 'VS-r'});
% set(gcf, 'Position', [123 123 333 222]);
% xlim([0.5 8.5])
% 
% %
% saveas(gcf,'cueConsistency.fig');

