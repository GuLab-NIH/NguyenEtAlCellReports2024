
%load data: contact lead author
%%
load('otherIdx.mat')
i=find(nRuns(:,1)>1);
allFieldsUse=allFields(i);

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
yy=AOVSF(:,2); xidx=2;

a=xx;
b=yy;
c={};
c{1}='AO';%common
c{2}='VS';%unique
color=[1 0 1;0 1 0];
% color=repmat(color,2,1);%all use the same color
%data
All=[];
for n=1:size(a,2);
    All=[All;a(:,n)];
end
for n=1:size(b,2);
    All=[All;b(:,n)];
end

%category
AllC={};
for n=1:size(a,2),
    AllC=[AllC;repmat(c(n),size(a,1),1)];
end

for n=1:size(b,2),
    AllC=[AllC;repmat(c(n+1),size(b,1),1)];
end
violinplot(All, AllC,'ShowData',false,'ViolinColor',color,'ShowMean',true,'ShowWhiskers',false,'ShowBox',false);
xlim tight;

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);
title(['field coverage',num2str(p)]); ylabel('% cells');
saveas(gcf,'fieldCoverage_violin.fig')

%% rbr consistency
load('cueConsistenciesYG.mat');
fovIndices = find(nRuns(:,1)>1);
load('rewardConsistency.mat')
rewardConsistencyUse=rewardConsistency(fovIndices);
rewardConsistencyUse=cell2mat(rewardConsistencyUse');

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
xx=nanmean(xxr,2);

% vs
iidx = [14 86 100 132]; xidx = 2;
yy=cell2mat(cueConsistenciesYG(fovIndices)');
yy = yy(:, iidx);

yyr=[];
yyr(:,1)=yy(:,1);
yyr(:,2)=rewardConsistencyUse(:,2);
yyr(:,3:5)=yy(:,2:4);
yyr=yyr(otherIdx,:);
yy=nanmean(yyr,2);

a=xx;
b=yy;
c={};
c{1}='AO';%common
c{2}='VS';%unique
color=[1 0 1;0 1 0];
% color=repmat(color,2,1);%all use the same color
%data
All=[];
for n=1:size(a,2);
    All=[All;a(:,n)];
end
for n=1:size(b,2);
    All=[All;b(:,n)];
end

%category
AllC={};
for n=1:size(a,2),
    AllC=[AllC;repmat(c(n),size(a,1),1)];
end

for n=1:size(b,2),
    AllC=[AllC;repmat(c(n+1),size(b,1),1)];
end
violinplot(All, AllC,'ShowData',false,'ViolinColor',color,'ShowMean',true,'ShowWhiskers',false,'ShowBox',false);
xlim tight;

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);
title(['Runs Consistency',num2str(p)]);
saveas(gcf,'RBRConsistencyAtAOVS_violin.fig')

%% Mean-DFOF in AUDIO vs VISUAL

figure; hold on

% ao
iidx = [29:38 39:42 43:52 63:72 111:120]; xidx = 1;
xx = dfofUse; xx = xx(otherIdx, iidx);
xx=nanmean(xx,2);
% vs
iidx = [9:18 81:90 91:94 95:104 127:136]; xidx = 2;
yy = dfofUse; yy = yy(otherIdx, iidx);
yy=nanmean(yy,2);

a=xx;
b=yy;
c={};
c{1}='AO';%common
c{2}='VS';%unique
color=[1 0 1;0 1 0];
% color=repmat(color,2,1);%all use the same color
%data
All=[];
for n=1:size(a,2);
    All=[All;a(:,n)];
end
for n=1:size(b,2);
    All=[All;b(:,n)];
end

%category
AllC={};
for n=1:size(a,2),
    AllC=[AllC;repmat(c(n),size(a,1),1)];
end

for n=1:size(b,2),
    AllC=[AllC;repmat(c(n+1),size(b,1),1)];
end
violinplot(All, AllC,'ShowData',false,'ViolinColor',color,'ShowMean',true,'ShowWhiskers',false,'ShowBox',false);
xlim tight;

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);
title(['dfof',num2str(p)])

saveas(gcf,'dfofAOVS_violin.fig')

%% cue reward consistency

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
xx=cell2mat(cueConsistenciesYG(fovIndices)');
xx = xx(:, iidx);
%include reward
xxr=[];
xxr(:,1)=xx(:,1);
xxr(:,2)=rewardConsistencyUse(:,1);
xxr(:,3:5)=xx(:,2:4);
xxr=xxr(otherIdx,:);

xidx = [2 3 4 5 9];

% vs
iidx = [14 86 100 132];
yy=cell2mat(cueConsistenciesYG(fovIndices)');
yy = yy(:, iidx);

yyr=[];
yyr(:,1)=yy(:,1);
yyr(:,2)=rewardConsistencyUse(:,2);
yyr(:,3:5)=yy(:,2:4);
yyr=yyr(otherIdx,:);

yidx = [1 6 7 8 10];

AllColumn=[];

for n=1:length(xidx);
    % a=xxr(:,n);
    % a=a(~isnan(a))
    AllColumn(:,xidx(n))=xxr(:,n);
    % a=yyr(:,n);
    % a=a(~isnan(a))
    AllColumn(:,yidx(n))=(yyr(:,n));
end

All=[];
for n=1:size(AllColumn,2);
    a=AllColumn(:,n);
    a=a(~isnan(a));
    All=[All;a];
end

c={};
c{1}='v1';
c{2}='a1';
c{3}='r1';
c{4}='a2';
c{5}='a3';
c{6}='v2';
c{7}='r2';
c{8}='v3';
c{9}='a4';
c{10}='v4';

color=[0 1 0;1 0 1;0 0 0;1 0 1;1 0 1;0 1 0;0 0 0;0 1 0;1 0 1;0 1 0];
%category
AllC={};
for n=1:size(AllColumn,2);
    a=AllColumn(:,n);
    sa=length(find(~isnan(a)));
    AllC=[AllC;repmat(c(n),sa,1)];

end

grouporder={'v1','a1','r1','a2','a3','v2','r2','v3','a4','v4'};
figure
violinplot(All, AllC,'ShowData',false,'ViolinColor',color,'ShowMean',true,'ShowWhiskers',false,'ShowBox',false,'GroupOrder',grouporder);

zz=[xxr yyr];
zz = nanmean(zz, 2);
mm = nanmean(zz);
se = nanstd(zz) / sqrt(numel(zz));
plot([0.5 10.5], [mm mm], 'k--');
% plot([0.5 10.5], [mm+se mm+se], 'k--');
% plot([0.5 10.5], [mm-se mm-se], 'k--');

xlim tight;

saveas(gcf,'cueRewardConsistency_violin.fig')




