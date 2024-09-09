%load data: contact lead author
%%
r1=[25:28];
nr1=setdiff([1:1:42],r1);
r2=[49:52];
nr2=setdiff([1:1:66],r2);

percentMouse=[];
percentMouseS=[];
for mouse=1:6;
df2 = cueDfofs{mouse, 1}';
df1 = cueDfofs{mouse, 2}';

rdf1 = df1(:, r1); rdf1 = nanmean(rdf1, 2);
nrdf1 = df1(:, nr1); nrdf1 = nanmean(nrdf1, 2);
d1 = (rdf1 - nrdf1) ./ (rdf1 + nrdf1);

rdf2 = df2(:, r2); rdf2 = nanmean(rdf2, 2);
nrdf2 = df2(:, nr2); nrdf2 = nanmean(nrdf2, 2);
d2 = (rdf2 - nrdf2) ./ (rdf2 + nrdf2);

i1=find(d1>thresh99d1);
i2=find(d2>thresh99d2);
ii=intersect(i1,i2);
percent=length(ii)/length(d1);

percentMouse(mouse)=percent;

% random
percentS=[];
nShuffle=1000;
for n=1:nShuffle;
    d2Fake=d2(randperm(length(d2)));
    i2Fake=find(d2Fake>thresh99d2);
    iiFake=intersect(i1,i2Fake);
    percentS(n)=length(iiFake)/length(d1);
end
percentMouseS(:,mouse)=percentS;
end

figure

subplot(261)
pr=mean(percentMouse);
s=mean(percentMouseS,2);
[p,x]=ksdensity(s);
plot(x,p);
hold on
line([pr pr],[0 max(p)])
pp=length(find(s>pr))/nShuffle;
title(['p=',num2str(pp)])

subplot(262)

mouseRealS=[];
mouseRealS(:,1)=percentMouse;
mouseRealS(:,2)=mean(percentMouseS,1);

for n=1:6;
    hold on
    plot(mouseRealS(n,:));
end
xlim([0.5 2.5])

[~,p]=ttest(mouseRealS(:,1),mouseRealS(:,2));
title(['p=',num2str(p)])

for n=1:6;
    subplot(2,6, 6+n)
    s=percentMouseS(:,n);
    [p,x]=ksdensity(s);
    plot(x,p);
    hold on
    line([percentMouse(n) percentMouse(n)],[0 max(p)]);
    pp=length(find(s>percentMouse(n)))/nShuffle;
title(['p=',num2str(pp)])
end
 

saveas(gcf,'avr1avr2Reward_mouse.fig')

%%
load('fig6_avr1vvr1_10bin.mat');
r1=[25:28];
nr1=setdiff([1:1:42],r1);
r2=[25:28];
nr2=setdiff([1:1:42],r1);

percentMouse=[];
percentMouseS=[];
load('thresha1v1_rewardZone.mat')
for mouse=1:6;
df2 = cueDfofs{mouse, 1}';
df1 = cueDfofs{mouse, 2}';

rdf1 = df1(:, r1); rdf1 = nanmean(rdf1, 2);
nrdf1 = df1(:, nr1); nrdf1 = nanmean(nrdf1, 2);
d1 = (rdf1 - nrdf1) ./ (rdf1 + nrdf1);

rdf2 = df2(:, r2); rdf2 = nanmean(rdf2, 2);
nrdf2 = df2(:, nr2); nrdf2 = nanmean(nrdf2, 2);
d2 = (rdf2 - nrdf2) ./ (rdf2 + nrdf2);

i1=find(d1>thresh99d1);
i2=find(d2>thresh99d2);
ii=intersect(i1,i2);
percent=length(ii)/length(d1);

percentMouse(mouse)=percent;

% random
percentS=[];
nShuffle=1000;
for n=1:nShuffle;
    d2Fake=d2(randperm(length(d2)));
    i2Fake=find(d2Fake>thresh99d2);
    iiFake=intersect(i1,i2Fake);
    percentS(n)=length(iiFake)/length(d1);
end
percentMouseS(:,mouse)=percentS;
end

figure

subplot(261)
pr=mean(percentMouse);
s=mean(percentMouseS,2);
[p,x]=ksdensity(s);
plot(x,p);
hold on
line([pr pr],[0 max(p)])
pp=length(find(s>pr))/nShuffle;
title(['p=',num2str(pp)])

subplot(262)

mouseRealS=[];
mouseRealS(:,1)=percentMouse;
mouseRealS(:,2)=mean(percentMouseS,1);

for n=1:6;
    hold on
    plot(mouseRealS(n,:));
end
xlim([0.5 2.5])
[~,p]=ttest(mouseRealS(:,1),mouseRealS(:,2));
title(['p=',num2str(p)])

for n=1:6;
    subplot(2,6, 6+n)
    s=percentMouseS(:,n);
    [p,x]=ksdensity(s);
    plot(x,p);
    hold on
    line([percentMouse(n) percentMouse(n)],[0 max(p)]);
    pp=length(find(s>percentMouse(n)))/nShuffle;
title(['p=',num2str(pp)])
end

saveas(gcf,'avr1vvr1Reward_mouse.fig')

%%
load('fig6_vvr1vvr2_10bin.mat');
r1=[25:28];
nr1=setdiff([1:1:42],r1);
r2=[49:52];
nr2=setdiff([1:1:66],r2);


percentMouse=[];
percentMouseS=[];
load('threshv1v2_rewardZone.mat')
for mouse=1:6;
df2 = cueDfofs{mouse, 1}';
df1 = cueDfofs{mouse, 2}';

rdf1 = df1(:, r1); rdf1 = nanmean(rdf1, 2);
nrdf1 = df1(:, nr1); nrdf1 = nanmean(nrdf1, 2);
d1 = (rdf1 - nrdf1) ./ (rdf1 + nrdf1);

rdf2 = df2(:, r2); rdf2 = nanmean(rdf2, 2);
nrdf2 = df2(:, nr2); nrdf2 = nanmean(nrdf2, 2);
d2 = (rdf2 - nrdf2) ./ (rdf2 + nrdf2);

i1=find(d1>thresh99d1);
i2=find(d2>thresh99d2);
ii=intersect(i1,i2);
percent=length(ii)/length(d1);

percentMouse(mouse)=percent;

% random
percentS=[];
nShuffle=1000;
for n=1:nShuffle;
    d2Fake=d2(randperm(length(d2)));
    i2Fake=find(d2Fake>thresh99d2);
    iiFake=intersect(i1,i2Fake);
    percentS(n)=length(iiFake)/length(d1);
end
percentMouseS(:,mouse)=percentS;
end

figure

subplot(261)
pr=mean(percentMouse);
s=mean(percentMouseS,2);
[p,x]=ksdensity(s);
plot(x,p);
hold on
line([pr pr],[0 max(p)])
pp=length(find(s>pr))/nShuffle;
title(['p=',num2str(pp)])

subplot(262)

mouseRealS=[];
mouseRealS(:,1)=percentMouse;
mouseRealS(:,2)=mean(percentMouseS,1);

for n=1:6;
    hold on
    plot(mouseRealS(n,:));
end
xlim([0.5 2.5])
[~,p]=ttest(mouseRealS(:,1),mouseRealS(:,2));
title(['p=',num2str(p)])

for n=1:6;
    subplot(2,6, 6+n)
    s=percentMouseS(:,n);
    [p,x]=ksdensity(s);
    plot(x,p);
    hold on
    line([percentMouse(n) percentMouse(n)],[0 max(p)]);
    pp=length(find(s>percentMouse(n)))/nShuffle;
title(['p=',num2str(pp)])
end

saveas(gcf,'vvr1vvr2Reward_mouse.fig')
%%
load('fig6_avr2vvr2_10bin.mat');
r1=[49:52];
nr1=setdiff([1:1:66],r2);
r2=[49:52];
nr2=setdiff([1:1:66],r2);

percentMouse=[];
percentMouseS=[];
load('thresha2v2_rewardZone.mat')
for mouse=1:6;
df2 = cueDfofs{mouse, 1}';
df1 = cueDfofs{mouse, 2}';

rdf1 = df1(:, r1); rdf1 = nanmean(rdf1, 2);
nrdf1 = df1(:, nr1); nrdf1 = nanmean(nrdf1, 2);
d1 = (rdf1 - nrdf1) ./ (rdf1 + nrdf1);

rdf2 = df2(:, r2); rdf2 = nanmean(rdf2, 2);
nrdf2 = df2(:, nr2); nrdf2 = nanmean(nrdf2, 2);
d2 = (rdf2 - nrdf2) ./ (rdf2 + nrdf2);

i1=find(d1>thresh99d1);
i2=find(d2>thresh99d2);
ii=intersect(i1,i2);
percent=length(ii)/length(d1);

percentMouse(mouse)=percent;

% random
percentS=[];
nShuffle=1000;
for n=1:nShuffle;
    d2Fake=d2(randperm(length(d2)));
    i2Fake=find(d2Fake>thresh99d2);
    iiFake=intersect(i1,i2Fake);
    percentS(n)=length(iiFake)/length(d1);
end
percentMouseS(:,mouse)=percentS;
end

figure

subplot(261)
pr=mean(percentMouse);
s=mean(percentMouseS,2);
[p,x]=ksdensity(s);
plot(x,p);
hold on
line([pr pr],[0 max(p)])
pp=length(find(s>pr))/nShuffle;
title(['p=',num2str(pp)])

subplot(262)

mouseRealS=[];
mouseRealS(:,1)=percentMouse;
mouseRealS(:,2)=mean(percentMouseS,1);

for n=1:6;
    hold on
    plot(mouseRealS(n,:));
end
xlim([0.5 2.5])
[~,p]=ttest(mouseRealS(:,1),mouseRealS(:,2));
title(['p=',num2str(p)])

for n=1:6;
    subplot(2,6, 6+n)
    s=percentMouseS(:,n);
    [p,x]=ksdensity(s);
    plot(x,p);
    hold on
    line([percentMouse(n) percentMouse(n)],[0 max(p)]);
    pp=length(find(s>percentMouse(n)))/nShuffle;
    % [r,pp]=ttest(s,percentMouse(n));
title(['p=',num2str(pp)])
end
saveas(gcf,'avr2vvr2Reward_mouse.fig')

%%   for paper figures

r1=[25:28];
nr1=setdiff([1:1:42],r1);
r2=[49:52];
nr2=setdiff([1:1:66],r2);

temp1=zeros(1,42);
temp1(r1)=1;
temp2=zeros(1,66);
temp2(r2)=1;

figure,

subplot(211)
r=r1;
line([1 r(1)],[0 0])
hold on
line([r(1) r(end)],[1 1])
hold on
line([r(end) 42],[0 0])
hold on
line([r(1) r(1)],[0 1])
hold on
line([r(end) r(end)],[0 1])

xlim([1 66])
axis off

subplot(212)
r=r2;
line([1 r(1)],[0 0])
hold on
line([r(1) r(end)],[1 1])
hold on
line([r(end) 66],[0 0])
hold on
line([r(1) r(1)],[0 1])
hold on
line([r(end) r(end)],[0 1])
xlim([1 66])
axis off

saveas(gcf,'rewardTemp.fig')

%example of cells with thresholds
load('fig6_avr1avr2_10bin.mat');
r1=[25:28];
nr1=setdiff([1:1:42],r1);
r2=[49:52];
nr2=setdiff([1:1:66],r2);

percentMouse=[];
percentMouseS=[];
load('thresha1a2_rewardZone.mat')
for mouse=6;
df2 = cueDfofs{mouse, 1}';
df1 = cueDfofs{mouse, 2}';

rdf1 = df1(:, r1); rdf1 = nanmean(rdf1, 2);
nrdf1 = df1(:, nr1); nrdf1 = nanmean(nrdf1, 2);
d1 = (rdf1 - nrdf1) ./ (rdf1 + nrdf1);

rdf2 = df2(:, r2); rdf2 = nanmean(rdf2, 2);
nrdf2 = df2(:, nr2); nrdf2 = nanmean(nrdf2, 2);
d2 = (rdf2 - nrdf2) ./ (rdf2 + nrdf2);

i1=find(d1>thresh99d1);
i2=find(d2>thresh99d2);
ii=intersect(i1,i2);
percent=length(ii)/length(d1);

percentMouse(mouse)=percent;

% random
percentS=[];
nShuffle=1000;
for n=1:nShuffle;
    d2Fake=d2(randperm(length(d2)));
    i2Fake=find(d2Fake>thresh99d2);
    iiFake=intersect(i1,i2Fake);
    percentS(n)=length(iiFake)/length(d1);
end
percentMouseS(:,mouse)=percentS;

figure
subplot(121)
plot(d1,d2,'k.')
hold on
plot(d1(i1),d2(i1),'m.')
hold on
plot(d1(i2),d2(i2),'m.')
hold on
plot(d1(ii),d2(ii),'c.')
hold on
line([thresh99d1 thresh99d1],[-1 1])
hold on
line([-1 1],[thresh99d2 thresh99d2])

subplot(122)
  s=percentMouseS(:,mouse);
    [p,x]=ksdensity(s);
    plot(x,p);
    hold on
    line([percentMouse(mouse) percentMouse(mouse)],[0 max(p)]);
    pp=length(find(s>percentMouse(mouse)))/nShuffle;
title(['p=',num2str(pp)])

saveas(gcf,'exampleOfMoreRewardCells.fig')
end

%%
%cell example
figure,
subplot(211)
a=df1(ii(18),:);
plot(a)
m=max(a);
hold on
r=r1;
line([1 r(1)],[0 0])
hold on
line([r(1) r(end)],[m m])
hold on
line([r(end) 42],[0 0])
hold on
line([r(1) r(1)],[0 m])
hold on
line([r(end) r(end)],[0 m])

xlim([1 42])


subplot(212)
a=df2(ii(9),:);
plot(a)
m=max(a);
hold on
r=r2;
line([1 r(1)],[0 0])
hold on
line([r(1) r(end)],[m m])
hold on
line([r(end) 66],[0 0])
hold on
line([r(1) r(1)],[0 m])
hold on
line([r(end) r(end)],[0 m])

xlim([1 66])
saveas(gcf,'exampleCellsMoreRewardCells.fig')

for mouse=2;
df2 = cueDfofs{mouse, 1}';
df1 = cueDfofs{mouse, 2}';

rdf1 = df1(:, r1); rdf1 = nanmean(rdf1, 2);
nrdf1 = df1(:, nr1); nrdf1 = nanmean(nrdf1, 2);
d1 = (rdf1 - nrdf1) ./ (rdf1 + nrdf1);

rdf2 = df2(:, r2); rdf2 = nanmean(rdf2, 2);
nrdf2 = df2(:, nr2); nrdf2 = nanmean(nrdf2, 2);
d2 = (rdf2 - nrdf2) ./ (rdf2 + nrdf2);

i1=find(d1>thresh99d1);
i2=find(d2>thresh99d2);
ii=intersect(i1,i2);
percent=length(ii)/length(d1);

percentMouse(mouse)=percent;

% random
percentS=[];
nShuffle=1000;
for n=1:nShuffle;
    d2Fake=d2(randperm(length(d2)));
    i2Fake=find(d2Fake>thresh99d2);
    iiFake=intersect(i1,i2Fake);
    percentS(n)=length(iiFake)/length(d1);
end
percentMouseS(:,mouse)=percentS;

figure
subplot(121)
plot(d1,d2,'k.')
hold on
plot(d1(i1),d2(i1),'m.')
hold on
plot(d1(i2),d2(i2),'m.')
hold on
plot(d1(ii),d2(ii),'c.')
hold on
line([thresh99d1 thresh99d1],[-1 1])
hold on
line([-1 1],[thresh99d2 thresh99d2])

subplot(122)
  s=percentMouseS(:,mouse);
    [p,x]=ksdensity(s);
    plot(x,p);
    hold on
    line([percentMouse(mouse) percentMouse(mouse)],[0 max(p)]);
    pp=length(find(s>percentMouse(mouse)))/nShuffle;
title(['p=',num2str(pp)])

saveas(gcf,'exampleOfLessRewardCells.fig')
end

%%

%example of cells with thresholds _ different sensory
load('fig6_avr1vvr1_10bin.mat');
r1=[25:28];
nr1=setdiff([1:1:42],r1);
r2=[25:28];
nr2=setdiff([1:1:42],r1);

percentMouse=[];
percentMouseS=[];
load('thresha1v1_rewardZone.mat')

%% plot
for mouse=5;
df2 = cueDfofs{mouse, 1}';
df1 = cueDfofs{mouse, 2}';

rdf1 = df1(:, r1); rdf1 = nanmean(rdf1, 2);
nrdf1 = df1(:, nr1); nrdf1 = nanmean(nrdf1, 2);
d1 = (rdf1 - nrdf1) ./ (rdf1 + nrdf1);

rdf2 = df2(:, r2); rdf2 = nanmean(rdf2, 2);
nrdf2 = df2(:, nr2); nrdf2 = nanmean(nrdf2, 2);
d2 = (rdf2 - nrdf2) ./ (rdf2 + nrdf2);

i1=find(d1>thresh99d1);
i2=find(d2>thresh99d2);
ii=intersect(i1,i2);
percent=length(ii)/length(d1);

percentMouse(mouse)=percent;

% random
percentS=[];
nShuffle=1000;
for n=1:nShuffle;
    d2Fake=d2(randperm(length(d2)));
    i2Fake=find(d2Fake>thresh99d2);
    iiFake=intersect(i1,i2Fake);
    percentS(n)=length(iiFake)/length(d1);
end
percentMouseS(:,mouse)=percentS;

figure
subplot(121)
plot(d1,d2,'k.')
hold on
plot(d1(i1),d2(i1),'m.')
hold on
plot(d1(i2),d2(i2),'m.')
hold on
plot(d1(ii),d2(ii),'c.')
hold on
line([thresh99d1 thresh99d1],[-1 1])
hold on
line([-1 1],[thresh99d2 thresh99d2])

subplot(122)
  s=percentMouseS(:,mouse);
    [p,x]=ksdensity(s);
    plot(x,p);
    hold on
    line([percentMouse(mouse) percentMouse(mouse)],[0 max(p)]);
    pp=length(find(s>percentMouse(mouse)))/nShuffle;
title(['p=',num2str(pp)])

end

saveas(gcf,'exampleOfMoreRewardCells_diffSensory.fig')
