%load data: contact lead author
%%
figure; 
subplot(121)
hold on
xx1 = overlapPerc{1,1} * 100; 
xx2 = overlapPerc{1,2} * 100; 

a=xx1';
b=xx2';
c={};
c{1}='AVR';%common
c{2}='VVR';%unique
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

[~,pa] = ttest( xx1, xx2 ); disp(pa);
title(['audio p=',num2str(pa)])

%%
subplot(122)
hold on
yy1 = overlapPerc{2,1} * 100; 
yy2 = overlapPerc{2,2} * 100;  

a=yy1';
b=yy2';
c={};
c{1}='AVR';%common
c{2}='VVR';%unique
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

[~,pv] = ttest( yy1, yy2 ); disp(pv);
title(['visual p=',num2str(pv)])

saveas(gcf,'overlapPerc_cueReward_sameNormAsFig3_violin.fig')