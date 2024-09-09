%load data: contact lead author 
%% MEAN DFOF
i=find(nRuns(:,1)>1);

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

% saveas(gcf,'meandfof_violin.fig');

% %% average field coverage
% i=find(nRuns(:,1)>1);
% allDfofsUse=allDfofs(i);
% %similar to the previousone, need to interplate allDfof
% for n=1:length(allDfofsUse);
%     for m=1:size(allDfofsUse{n},1);
%         allDfofsUse{n}(m,:)=naninterp(allDfofsUse{n}(m,:));
%     end
% end
% allFieldsUse=allFields(i);
% 
% aoArea = [29:38 43:52 63:72 111:120];
% vsArea = [9:18 81:90 95:104 127:136];
% 
% dd = cellfun(@(x) (mean(x(:,aoArea),2) - mean(x(:,vsArea),2)) ...
%     ./ (mean(x(:,aoArea),2) + mean(x(:,vsArea),2)), allDfofsUse, 'uni', 0);
% 
% ddAO = cellfun(@(x, y) x(y > thresh95,:), allFieldsUse, dd, 'uni', 0);
% ddAO = cat(1, ddAO{:});
% ddAO = nanmean(ddAO, 2);
% 
% ddVS = cellfun(@(x, y) x(y < thresh5,:), allFieldsUse, dd, 'uni', 0);
% ddVS = cat(1, ddVS{:});
% ddVS = nanmean(ddVS, 2);
% 
% figure; hold on
% 
% % UNIQUE cells
% xx = cat(1, ddAO); xx = xx *100; % *100 for percentage
% yy = cat(1, ddVS); yy = yy *100; % *100 for percentage
% 
% a=xx;
% b=yy;
% c={};
% c{1}='ao';%common
% c{2}='vs';%unique
% color=[1 0 1;0 1 0];
% % color=repmat(color,2,1);%all use the same color
% %data
% All=[];
% for n=1:size(a,2);
%     All=[All;a(:,n)];
% end
% for n=1:size(b,2);
%     All=[All;b(:,n)];
% end
% 
% %category
% AllC={};
% for n=1:size(a,2),
%     AllC=[AllC;repmat(c(n),size(a,1),1)];
% end
% 
% for n=1:size(b,2),
%     AllC=[AllC;repmat(c(n+1),size(b,1),1)];
% end
% violinplot(All, AllC,'ShowData',false,'ViolinColor',color,'ShowMean',true,'ShowWhiskers',false,'ShowBox',false);
% xlim tight;
% 
% [~,p] = ttest2(xx, yy);
% title(['field coverage, p=',num2str(p)])
% 
% saveas(gcf,'fieldCoverage_violin.fig');
% 

%% average field coverage per cell: the above calcualtion is to complicated, we just want to calculate this on a per cell basis 
i=find(nRuns(:,1)>1);
% allDfofsUse=allDfofs(i);
% %similar to the previousone, need to interplate allDfof
% for n=1:length(allDfofsUse);
%     for m=1:size(allDfofsUse{n},1);
%         allDfofsUse{n}(m,:)=naninterp(allDfofsUse{n}(m,:));
%     end
% end
allFieldsUse=allFields(i);
allFieldsPer={};
for n=1:length(allFieldsUse);
    allFieldsPer{n}=[];
    for m=1:size(allFieldsUse{n},1);
        allFieldsPer{n}(m)=sum(allFieldsUse{n}(m,:))/size(allFieldsUse{n},2);
    end
end



% aoArea = [29:38 43:52 63:72 111:120];
% vsArea = [9:18 81:90 95:104 127:136];

ddAO = cat(2, allFieldsPer{:});
ddAO = ddAO(aoCellIdx)*100;

ddVS = cat(2, allFieldsPer{:});
ddVS = ddVS(vsCellIdx)*100;
% 
% dd = cellfun(@(x) (mean(x(:,aoArea),2) - mean(x(:,vsArea),2)) ...
%     ./ (mean(x(:,aoArea),2) + mean(x(:,vsArea),2)), allDfofsUse, 'uni', 0);
% 
% ddAO = cellfun(@(x, y) x(y > thresh95,:), allFieldsUse, dd, 'uni', 0);
% ddAO = cat(1, ddAO{:});
% ddAO = nanmean(ddAO, 2);
% 
% ddVS = cellfun(@(x, y) x(y < thresh5,:), allFieldsUse, dd, 'uni', 0);
% ddVS = cat(1, ddVS{:});
% ddVS = nanmean(ddVS, 2);

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

saveas(gcf,'fieldCoverage_violin_useThis.fig');
%this final plot is actually equivalent to "fieldCoverage_violin.fig' above

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

saveas(gcf,'rbrconsistency_violin.fig');

%% PLOT - cue consistencies direct comparison
cueConsistenciesUse=cueConsistencies(:,i);
% 10-bin
aoArea = [29:38 43:52 63:72 111:120];
vsArea = [9:18 81:90 95:104 127:136];

dd = cellfun(@(x) (mean(x(:,aoArea),2) - mean(x(:,vsArea),2)) ...
    ./ (mean(x(:,aoArea),2) + mean(x(:,vsArea),2)), allDfofsUse, 'uni', 0);

ddAO = cellfun(@(x, y) x(y > thresh95,:), cueConsistenciesUse([1 3],:), [dd; dd], 'uni', 0);
ddAO = cat(1, ddAO{:});
ddAO = ddAO(:, [34 48 68 116 14 86 100 132]);

ddVS = cellfun(@(x, y) x(y < thresh5,:), cueConsistenciesUse([1 3],:), [dd; dd], 'uni', 0);
ddVS = cat(1, ddVS{:});
ddVS = ddVS(:, [34 48 68 116 14 86 100 132]);

xidx = [2 3 4 7 1 5 6 8];
xx=[];%AO
yy=[];%VS
for n=1:length(xidx);
    xx(:,xidx(n))=ddAO(:,n);
    yy(:,xidx(n))=ddVS(:,n);
end

All=[];
for n=1:size(xx,2);
    a=xx(:,n);
    a=a(~isnan(a));
    All=[All;a];
    v=yy(:,n);
    v=v(~isnan(v));
    All=[All;v];
end

c={};
c{1}='a1';%ao cue1
c{2}='v1';%vs cue1
c{3}='a2';%ao cue2
c{4}='v2';%vs cue2
c{5}='a3';%ao cue3
c{6}='v3';%vs cue3
c{7}='a4';%ao cue4
c{8}='v4';%vs cue4
c{9}='a5';%ao cue5
c{10}='v5';%vs cue5
c{11}='a6';%ao cue6
c{12}='v6';%vs cue6
c{13}='a7';%ao cue7
c{14}='v7';%vs cue7
c{15}='a8';%ao cue8
c{16}='v8';%vs cue8

color=[1 0 1;0 1 0];
color=repmat(color,8,1);%all use the same color

%category
AllC={};
for n=1:size(xx,2),
    a=xx(:,n);
    sa=length(find(~isnan(a)));
    AllC=[AllC;repmat(c(2*(n-1)+1),sa,1)];
    v=yy(:,n);
    sv=length(find(~isnan(v)));
    AllC=[AllC;repmat(c(2*n),sv,1)];
end
grouporder={'a1','v1','a2','v2','a3','v3','a4','v4','a5','v5','a6','v6','a7','v7','a8','v8'};

figure
violinplot(All, AllC,'ShowData',false,'ViolinColor',color,'ShowMean',true,'ShowWhiskers',false,'ShowBox',false,'GroupOrder',grouporder);
xlim tight;
saveas(gcf,'cueConsistency_violin.fig');