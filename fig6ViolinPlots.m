%load data: contact lead author

%%
aoCellIdx=aoRewardIdx;
vsCellIdx=vsRewardIdx;
i=find(nRuns(:,1)>1);
%% rbr consistencies

rbrConsistenciesUse=rbrConsistencies(i);

ddAO = cat(1, rbrConsistenciesUse{:});
% ddAO=ddAO(iUse);
ddAO = ddAO(aoCellIdx);

ddVS = cat(1, rbrConsistenciesUse{:});
% ddVS=ddVS(iUse);
ddVS = ddVS(vsCellIdx);

figure; hold on
a=ddAO;
b=ddVS;
c={};
c{1}='ao';%common
c{2}='vs';%unique
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

[~,p] = ttest2(ddAO, ddVS); disp(p);

title(['rbr consistency, p=',num2str(p)])

saveas(gcf,'rbrconsistencyReward_violin.fig');

%% MEAN DFOF
allAverageDfofsUse=allAverageDfofs(i);


ddAO = cat(2, allAverageDfofsUse{:});
% ddAO=ddAO(iUse);
ddAO = ddAO(aoCellIdx);

ddVS = cat(2, allAverageDfofsUse{:});
% ddVS=ddVS(iUse);
ddVS = ddVS(vsCellIdx);

figure; hold on
a=ddAO';
b=ddVS';
c={};
c{1}='ao';%common
c{2}='vs';%unique
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

[~,p] = ttest2(ddAO, ddVS); disp(p);
title(['mean dfof, p=',num2str(p)])

saveas(gcf,'meandfofReward_violin.fig');

%% average field coverage (did not plot previously)
allFieldsUse=allFields(i);
allFieldsUse=allFields(i);
allFieldsPer={};
for n=1:length(allFieldsUse);
    allFieldsPer{n}=[];
    for m=1:size(allFieldsUse{n},1);
        allFieldsPer{n}(m)=sum(allFieldsUse{n}(m,:))/size(allFieldsUse{n},2);
    end
end

ddAO = cat(2, allFieldsPer{:});
ddAO = ddAO(aoCellIdx)*100;

ddVS = cat(2, allFieldsPer{:});
ddVS = ddVS(vsCellIdx)*100;


figure; hold on

% % UNIQUE cells
% xx = cat(1, ddAO); xx = xx *100; % *100 for percentage
% yy = cat(1, ddVS); yy = yy *100; % *100 for percentage

a=ddAO';
b=ddVS';
c={};
c{1}='ao';%common
c{2}='vs';%unique
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

[~,p] = ttest2(ddAO, ddVS); disp(p);
title(['field coverage, p=',num2str(p)])

saveas(gcf,'fieldCoverage_violin.fig');
