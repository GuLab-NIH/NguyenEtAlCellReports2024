load('fig5b.mat');

% data structure:
% 
% Y:\labMembers\DN\_PROJECT_\finalData\figures\20230705_v10\fig5b.mat
% 
% 4x2xN
% 2 is common-unique cells: 1 is common, 2 is unique
% 
% 
% N is number of FOVs
% 
% cueConsistency_FOV_1 is avr1-vvr1
% cueConsistency_FOV_2 is avr2-vvr2
% 
% 1 is the cells in AVR, cue consistnecy calculated using 10-bin area
% 
% 2 is the cells in VVR, cue consistnecy calculated using 6-bin area
% 
% 3 is the cells in AVR, cue consistnecy calculated using 6-bin area
% 
% 4 is the cells in VVR, cue consistnecy calculated using 10-bin area
% 
% in figures, used 1 and 4.
% 
% zones = [6 20 34]; xidx = 1:3;
% cuesCorr = cueConsistency_FOV_1;
% 
% xx = cat(1, cuesCorr{1,1,:}); xx = xx(:, zones);
% zones = [6 20 34];
% zones = [10 26 44 58];

%% AVR1
A=cueConsistency_FOV_1;
zones = [6 20 34];
AC = cat(1, A{1,1,:});
% i=setdiff([1:1:31],[2 3 31]);%at least 40 cells
% AC = cat(1, A{1,1,i});

AC = AC(:, zones);
B=mean(AC(:,[end-1 end]),2);
ACM=[mean(AC(:,1),2) B];
% ACMR=ACM(:,2)-ACM(:,1);


AU = cat(1, A{1,2,:});
% AU = cat(1, A{1,2,i});
AU = AU(:, zones);
B=mean(AU(:,[end-1 end]),2);
AUM=[mean(AU(:,1),2) B];
% AUMR=AUM(:,2)-AUM(:,1);

figure

a=ACM;
b=AUM;
c={};
c{1}='ca';%commoncell away
c{2}='cn';%commoncell near
c{3}='ua';%unique away
c{4}='un';%unique near
color=[1 0 1]; 
color=repmat(color,4,1);%all use the same color

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
    AllC=[AllC;repmat(c(n+2),size(b,1),1)];
end
   
violinplot(All, AllC,'ShowData',false,'ViolinColor',color,'ShowMean',true,'ShowWhiskers',false,'ShowBox',false);
p=[];
for n=1:size(ACM,2);
[~,p(n)]=ttest2(ACM(:,n),AUM(:,n));
end
title([num2str(p)])
xlim tight
saveas(gcf,'avr1CommonUnique_violin.fig')

%% AVR2
A=cueConsistency_FOV_2;
zones = [10 26 44 58];
% i=setdiff([1:1:31],[1 2 3 4]);%at least 40 cells
AC = cat(1, A{1,1,:});
% AC = cat(1, A{1,1,i});

AC = AC(:, zones);
B=mean(AC(:,[end-1 end]),2);
ACM=[mean(AC(:,[1 2]),2) B];

AU = cat(1, A{1,2,:});
% AU = cat(1, A{1,2,i});
AU = AU(:, zones);
B=mean(AU(:,[end-1 end]),2);
AUM=[mean(AU(:,[1 2]),2) B];

figure

a=ACM;
b=AUM;
c={};
c{1}='ca';%commoncell away
c{2}='cn';%commoncell near
c{3}='ua';%unique away
c{4}='un';%unique near
color=[1 0 1];
color=repmat(color,4,1);%all use the same color

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
    AllC=[AllC;repmat(c(n+2),size(b,1),1)];
end
   
violinplot(All, AllC,'ShowData',false,'ViolinColor',color,'ShowMean',true,'ShowWhiskers',false,'ShowBox',false);
p=[];
for n=1:size(ACM,2);
[~,p(n)]=ttest2(ACM(:,n),AUM(:,n));
end
title([num2str(p)])
xlim tight
saveas(gcf,'avr2CommonUnique_violin.fig')

%% VVR1
A=cueConsistency_FOV_1;
zones = [6 20 34];
AC = cat(1, A{4,1,:});
AC = AC(:, zones);
B=mean(AC(:,[end-1 end]),2);
ACM=[mean(AC(:,1),2) B];

AU = cat(1, A{4,2,:});
AU = AU(:, zones);
B=mean(AU(:,[end-1 end]),2);
AUM=[mean(AU(:,1),2) B];

figure
a=ACM;
b=AUM;
c={};
c{1}='ca';%commoncell away
c{2}='cn';%commoncell near
c{3}='ua';%unique away
c{4}='un';%unique near
color=[0 1 0];
color=repmat(color,4,1);%all use the same color

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
    AllC=[AllC;repmat(c(n+2),size(b,1),1)];
end
   
violinplot(All, AllC,'ShowData',false,'ViolinColor',color,'ShowMean',true,'ShowWhiskers',false,'ShowBox',false);
p=[];
for n=1:size(ACM,2);
[~,p(n)]=ttest2(ACM(:,n),AUM(:,n));
end
title([num2str(p)])
xlim tight
saveas(gcf,'vvr1CommonUnique_violin.fig')

%% VVR2
A=cueConsistency_FOV_2;
zones = [10 26 44 58];
AC = cat(1, A{4,1,:});
AC = AC(:, zones);
B=mean(AC(:,[end-1 end]),2);
ACM=[mean(AC(:,[1 2]),2) B];

AU = cat(1, A{4,2,:});
AU = AU(:, zones);
B=mean(AU(:,[end-1 end]),2);
AUM=[mean(AU(:,[1 2]),2) B];

figure
a=ACM;
b=AUM;
c={};
c{1}='ca';%commoncell away
c{2}='cn';%commoncell near
c{3}='ua';%unique away
c{4}='un';%unique near
color=[0 1 0];
color=repmat(color,4,1);%all use the same color

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
    AllC=[AllC;repmat(c(n+2),size(b,1),1)];
end
   
violinplot(All, AllC,'ShowData',false,'ViolinColor',color,'ShowMean',true,'ShowWhiskers',false,'ShowBox',false);
p=[];
for n=1:size(ACM,2);
[~,p(n)]=ttest2(ACM(:,n),AUM(:,n));
end
title([num2str(p)])
xlim tight
saveas(gcf,'vvr2CommonUnique_violin.fig')


%% all consistencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
figure

subplot(221)
A=rbrConsistency_1;
a=A{1,1}(:,3);
b=A{1,2}(:,3);
c={};
c{1}='c';%common
c{2}='u';%unique
color=[1 0 1];
color=repmat(color,2,1);%all use the same color

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

[~,p]=ttest2(A{1,1}(:,3),A{1,2}(:,3));
title(['avr1',num2str(p)])


subplot(222)
A=rbrConsistency_2;

a=A{1,1}(:,3);
b=A{1,2}(:,3);
c={};
c{1}='c';%common
c{2}='u';%unique
color=[1 0 1];
color=repmat(color,2,1);%all use the same color

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

[~,p]=ttest2(A{1,1}(:,3),A{1,2}(:,3));
title(['avr2',num2str(p)])


subplot(223)
A=rbrConsistency_1;
a=A{2,1}(:,3);
b=A{2,2}(:,3);
c={};
c{1}='c';%common
c{2}='u';%unique
color=[0 1 0];
color=repmat(color,2,1);%all use the same color

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

[~,p]=ttest2(A{2,1}(:,3),A{2,2}(:,3));
title(['vvr1',num2str(p)])


subplot(224)
A=rbrConsistency_2;
a=A{2,1}(:,3);
b=A{2,2}(:,3);
c={};
c{1}='c';%common
c{2}='u';%unique
color=[0 1 0];
color=repmat(color,2,1);%all use the same color

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

[~,p]=ttest2(A{2,1}(:,3),A{2,2}(:,3));
title(['vvr2',num2str(p)])

saveas(gcf,'commonUniqueConsistencyAll_violin.fig')

%% plot fields in cue/out cue: similar to above but move bars: avr together, vvr together
% AVR1-VVR1: common cells
% % inArea = [3:8 17:22 31:36]; %6 bins
inArea=[1:10 15:24 29:38];
rew = 25:28;
outArea = setdiff(1:42, [inArea rew]);
allFields = cellfun(@(x) nanmean(x,2), allFields_1, 'uni', 0);
allFields = cellfun(@(x) nanmean(x(inArea)) / nanmean(x(outArea)), allFields, 'uni', 0);
env = 1;

figure; 
subplot(121)
%avr1 common unique
hold on
xx = cat(2, allFields{1,1,:});
yy = cat(2, allFields{1,2,:});
xx(isinf(xx))=nan;
a=xx';
b=yy';
c={};
c{1}='avr1c';%common
c{2}='avr1u';%unique
color=[1 0 1];
color=repmat(color,2,1);%all use the same color
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

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);
title(['avr1 common unique',num2str(p)])

% vvr1 common unique

subplot(122)
hold on
xx = cat(2, allFields{2,1,:});
yy = cat(2, allFields{2,2,:});
xx(isinf(xx))=nan;

a=xx';
b=yy';
c={};
c{1}='vvr1c';%common
c{2}='vvr1u';%unique
color=[0 1 0];
color=repmat(color,2,1);%all use the same color
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

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);
title(['vvr1 common unique',num2str(p)])

saveas(gcf,'avr1vvr1_commonUniqueCompare_rearrangeYG_violin.fig')


%% plot fields in cue/out cue
% AVR2-VVR2: common cells
% % inArea = [3:8 17:22 31:36]; %6 bins
inArea=[5:14 21:30 39:48 53:62];%10bins

rew = 49:52;
outArea = setdiff(1:66, [inArea rew]);
allFields = cellfun(@(x) nanmean(x,2), allFields_2, 'uni', 0);
allFields = cellfun(@(x) nanmean(x(inArea)) / nanmean(x(outArea)), allFields, 'uni', 0);
env = 2;

figure; 
subplot(121)
%avr1 common unique
hold on
xx = cat(2, allFields{1,1,:});
yy = cat(2, allFields{1,2,:});
xx(isinf(xx))=nan;

a=xx';
b=yy';
c={};
c{1}='avr2c';%common
c{2}='avr2u';%unique
color=[1 0 1];
color=repmat(color,2,1);%all use the same color
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

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);
title(['avr2 common unique',num2str(p)])

% vvr1 common unique

subplot(122)
hold on
xx = cat(2, allFields{2,1,:});
yy = cat(2, allFields{2,2,:});
xx(isinf(xx))=nan;

a=xx';
b=yy';
c={};
c{1}='vvr2c';%common
c{2}='vvr2u';%unique
color=[0 1 0];
color=repmat(color,2,1);%all use the same color
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

xx = reshape(xx, 1, []); yy = reshape(yy, 1, []);
[~,p] = ttest(xx, yy); disp(p);
title(['vvr2 common unique',num2str(p)])

saveas(gcf,'avr2vvr2_commonUniqueCompare_rearrangeYG_violin.fig')


