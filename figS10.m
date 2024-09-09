%load data: contact lead author
%%
fields=dfofUse;
%% %all pairs of fake reward zones: 4 bins per zone, any distance throught any location of the track
%get the pairs
%create locations of every 4 bins
totalBins=140;
% all4Bins=[0 0];
% for n=1:140/4;
%     all4Bins(n+1,:)=[all4Bins(n,2)+1 n*4];
% end
% all4Bins=all4Bins(2:end,:);

all4Bins=[];
for n=1:140-3;
    all4Bins(n,:)=[n n+3];
end

all4BinPairs=[];%the first 2 columns are the start and end of the first reward zone, the second columns are the start and end of the second reward zone.

for n=1:size(all4Bins,1)-4;
    z1=all4Bins(n,:);
    for m=n+4:size(all4Bins,1);%to ensure that the bins are not overlapping
        z2=all4Bins(m,:);
        all4BinPairs(end+1,:)=[z1 z2];
    end
end

%remove the real zone
remove1=[39 42 91 94];
[~,i1,~]=intersect(all4BinPairs,remove1,'rows');
remove2=[91 94 39 42];
[~,i2,~]=intersect(all4BinPairs,remove2,'rows');

i=setdiff([1:1:size(all4BinPairs,1)],[i1 i2]);
all4BinPairsNR=all4BinPairs(i,:);%remove the current zones

%% %all pairs of fake reward zones: 4 bins per zone, any distance throught any location of the track
%identify cells at all distance combinations
%%%%%%%%%%%%%%%%%%do not re run! takes a long time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nTotalAll=[];
nAOAll=[];
nVSAll=[];
for shufflePairs=1:size(all4BinPairsNR,1)
    clear idx
    clear aoRewardIdx
    clear vsRewardIdx

    disp(shufflePairs)
ra=[all4BinPairsNR(shufflePairs,1):all4BinPairsNR(shufflePairs,2)];
rv=[all4BinPairsNR(shufflePairs,3):all4BinPairsNR(shufflePairs,4)];
r=[ra rv];

%classify reward cells and auditory vs visual reward cells: below is same as this code: Z:\labMembers\DN\_PROJECT_\finalData\figures\20230705_v10\fig_YG_DFOF_completeRandom_99\reward_useShuffleThresh2\rewardCell2_inReward.m
%reward scores
d=[];%dfof within cues before the reward compared to other cues
nr=setdiff([1:1:140],r);
rdfof = fields(:, r); rdfof = nanmean(rdfof, 2);
nrdfof = fields(:, nr); nrdfof = nanmean(nrdfof, 2);
d = (rdfof - nrdfof) ./ (rdfof + nrdfof);

allF=cell2mat(shF');
dShuffle=[];
rdfof = allF(:, r); rdfof = nanmean(rdfof, 2);
nrdfof = allF(:, nr); nrdfof = nanmean(nrdfof, 2);
dShuffle = (rdfof - nrdfof) ./ (rdfof + nrdfof);
% 
allFP=allF(abs(dShuffle)<=1,:);
%
dShuffle=dShuffle(abs(dShuffle)<=1);
threshReward=prctile(dShuffle,99);
thresh1=prctile(dShuffle,1);

allFP=allFP(dShuffle>threshReward,:);

idx=find(d>threshReward);
rewardIdx=idx;

dfofR=fields(idx,:);
%calculate their activity preferneces

radfof = dfofR(:, ra); radfof = nanmean(radfof, 2);
rvdfof = dfofR(:, rv); rvdfof = nanmean(rvdfof, 2);
dr = (radfof - rvdfof) ./ (radfof + rvdfof);

nshuffle=100;
allShuffles={};
gaussianWindow=gausswin(7,1); %smooth the data within reward zone
for n=1:nshuffle;
    a=[];
    for m=1:size(dfofR,1);
        thisdfof=dfofR(m,:);
        thisdfofR=thisdfof(r);
        thisdfofRR=thisdfofR(randperm(length(thisdfofR)));
 thisdfofRR=gaussianSmoothWithNan( thisdfofRR,gaussianWindow);
        thisdfof(ra)=thisdfofRR(1:length(ra));
        thisdfof(rv)=thisdfofRR(length(ra)+1:length(r));  
        

        a(m,:)=thisdfof;
    end
    allShuffles{n,1}=a;
end
all=cell2mat(allShuffles);

radfof = all(:, ra); radfof = nanmean(radfof, 2);
rvdfof = all(:, rv); rvdfof = nanmean(rvdfof, 2);
ds = (radfof - rvdfof) ./ (radfof + rvdfof);

avrewardThresh1=prctile(ds,1);
avrewardThresh99=prctile(ds,99);

ia=find(dr>avrewardThresh99);
iv=find(dr<avrewardThresh1);

%these are idx in reward cells. we need to find idx in all cells
aoRewardIdx=rewardIdx(ia);
vsRewardIdx=rewardIdx(iv);

nTotalAll(shufflePairs,1)=length(rewardIdx);
nAOAll(shufflePairs,1)=length(aoRewardIdx);
nVSAll(shufflePairs,1)=length(vsRewardIdx);
end

save('randomAllDistances.mat', 'nTotalAll','nAOAll','nVSAll','all4BinPairsNR')
%% the above shuffle may not make sense. use fixed distance as the real reward zone
%% only used the zone pairs with the same distances at the two reward zones. but the fake zone pairs can be any where on the track


%create locations of every 4 bins
totalBins=140;
% all4Bins=[0 0];
% for n=1:140/4;
%     all4Bins(n+1,:)=[all4Bins(n,2)+1 n*4];
% end
% all4Bins=all4Bins(2:end,:);

all4Bins=[];
for n=1:140-3;
    all4Bins(n,:)=[n n+3];
end

all4BinPairs=[];%the first 2 columns are the start and end of the first reward zone, the second columns are the start and end of the second reward zone.

for n=1:size(all4Bins,1)-52; %52 is the distance between the two real reward zone (91-39)
    z1=all4Bins(n,:);
        z2=all4Bins(n+52,:);
        all4BinPairs(end+1,:)=[z1 z2];
end
remove1=[39 42 91 94];
[~,i1,~]=intersect(all4BinPairs,remove1,'rows');
i=setdiff([1:1:size(all4BinPairs,1)],i1);
all4BinPairs(i,:);
%reverse
all4BinPairsR=[all4BinPairs(:,[3 4]) all4BinPairs(:,[1 2])];
all4BinPairsAll=[all4BinPairs;all4BinPairsR];

%% only used the zone pairs with the same distances at the two reward zones. but the fake zone pairs can be any where on the track
%get the shuffle calculations

%%%%%%%%%%%%%%%%%%do not re run! takes a long time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nTotalAll2=[];
nAOAll2=[];
nVSAll2=[];

for shufflePairs=1:size(all4BinPairsAll,1)
    clear idx
    clear aoRewardIdx
    clear vsRewardIdx

    disp(shufflePairs)
ra=[all4BinPairsAll(shufflePairs,1):all4BinPairsAll(shufflePairs,2)];
rv=[all4BinPairsAll(shufflePairs,3):all4BinPairsAll(shufflePairs,4)];
r=[ra rv];
%reward scores
d=[];%dfof within cues before the reward compared to other cues
nr=setdiff([1:1:140],r);
rdfof = fields(:, r); rdfof = nanmean(rdfof, 2);
nrdfof = fields(:, nr); nrdfof = nanmean(nrdfof, 2);
d = (rdfof - nrdfof) ./ (rdfof + nrdfof);

allF=cell2mat(shF');
dShuffle=[];
rdfof = allF(:, r); rdfof = nanmean(rdfof, 2);
nrdfof = allF(:, nr); nrdfof = nanmean(nrdfof, 2);
dShuffle = (rdfof - nrdfof) ./ (rdfof + nrdfof);
% 
allFP=allF(abs(dShuffle)<=1,:);
%
dShuffle=dShuffle(abs(dShuffle)<=1);
threshReward=prctile(dShuffle,99);
thresh1=prctile(dShuffle,1);

allFP=allFP(dShuffle>threshReward,:);

idx=find(d>threshReward);
rewardIdx=idx;

dfofR=fields(idx,:);
%calculate their activity preferneces

radfof = dfofR(:, ra); radfof = nanmean(radfof, 2);
rvdfof = dfofR(:, rv); rvdfof = nanmean(rvdfof, 2);
dr = (radfof - rvdfof) ./ (radfof + rvdfof);

nshuffle=1000;
allShuffles={};
gaussianWindow=gausswin(7,1); %smooth the data within reward zone
for n=1:nshuffle;
    a=[];
    for m=1:size(dfofR,1);
        thisdfof=dfofR(m,:);
        thisdfofR=thisdfof(r);
        thisdfofRR=thisdfofR(randperm(length(thisdfofR)));
 thisdfofRR=gaussianSmoothWithNan( thisdfofRR,gaussianWindow);
        thisdfof(ra)=thisdfofRR(1:length(ra));
        thisdfof(rv)=thisdfofRR(length(ra)+1:length(r));  
        

        a(m,:)=thisdfof;
    end
    allShuffles{n,1}=a;
end
all=cell2mat(allShuffles);

radfof = all(:, ra); radfof = nanmean(radfof, 2);
rvdfof = all(:, rv); rvdfof = nanmean(rvdfof, 2);
ds = (radfof - rvdfof) ./ (radfof + rvdfof);

avrewardThresh1=prctile(ds,1);
avrewardThresh99=prctile(ds,99);

ia=find(dr>avrewardThresh99);
iv=find(dr<avrewardThresh1);

%these are idx in reward cells. we need to find idx in all cells
aoRewardIdx=rewardIdx(ia);
vsRewardIdx=rewardIdx(iv);

nTotalAll2(shufflePairs,1)=length(rewardIdx);
nAOAll2(shufflePairs,1)=length(aoRewardIdx);
nVSAll2(shufflePairs,1)=length(vsRewardIdx);
end

save('randomFixDistances.mat', 'nTotalAll2','nAOAll2','nVSAll2','all4BinPairsAll')
%% plot all distances: this does not work well
%remove the zones that overlaps with the two reward zones

AO=nAOAll./nTotalAll;
VS=nVSAll./nTotalAll;
All=AO+VS;

aoZone = [29:38 43:52 63:72 111:120 39:42];
vsZone = [9:18 81:90 95:104 127:136 91:94];


%exact zone +/- 3 bin
aoZone = [26:41 40:55 60:75 108:123 36:45];
vsZone = [6:22 78:94 92:107 124:139 88:97];

overlap=[]; %calculate the overlap between the zones and ao and vs zone

for n=1:size(all4BinPairsNR)
    a=[all4BinPairsNR(n,1):all4BinPairsNR(n,2)];
    b=[all4BinPairsNR(n,3):all4BinPairsNR(n,4)];
    aao=intersect(a,aoZone);
    avs=intersect(a,vsZone);
    % if ~isempty(aao);
    if length(aao)>=2;%need to overlap at lest 2 bins
        overlap(n,1)=1;
    % elseif ~isempty(avs);
    elseif length(avs)>=2;
        overlap(n,1)=2;
    else
        overlap(n,1)=0;
    end

    bao=intersect(b,aoZone);
    bvs=intersect(b,vsZone);
    % if ~isempty(bao);
    if length(bao)>=2;
        overlap(n,2)=1;
    % elseif ~isempty(bvs);
    elseif length(bvs)>=2;
        overlap(n,2)=2;
    else
        overlap(n,2)=0;
    end
end

overlapSum=sum(overlap,2);
i=find(overlapSum==0);
%
figure,
[p,x]=ksdensity(All(i));
plot(x,p)
real=(length(aoRewardIdx)+length(vsRewardIdx))/length(rewardIdx);
hold on
line([real real],[0 max(p)])
per=length(find(All(i)<=real))/length(i);
title(['percent',num2str(per)])

overlapSum=sum(overlap,2);
overlapDiff=overlap(:,1)-overlap(:,2);
% i=find(overlapSum~=3);
i=find(overlapDiff==0); %same cue type or no cue
%
figure,
[p,x]=ksdensity(All(i));
plot(x,p)
real=(length(aoRewardIdx)+length(vsRewardIdx))/length(rewardIdx);
hold on
line([real real],[0 max(p)])
per=length(find(All(i)<=real))/length(i);
title(['percent',num2str(per)])
% 
% color=[1 0 0;1 0 1;0 1 0;0 0 1;0 0 0]; %r, m, g, b, k
% colorSum=[];
% for n=1:length(overlapSum);
%     for m=1:5;
%     if overlapSum(n)==m-1;
%         colorSum(n,:)=color(m,:);
%     end
%     end
% end
% 
% figure,
% for n=1:size(overlapSum,1);
%     hold on
%     plot(n,All(n),'.','Color',colorSum(n,:));
% end

%% plot fixed distances: used this one

%remove the zones that overlaps with the two reward zones

AO=nAOAll2./nTotalAll2;
VS=nVSAll2./nTotalAll2;
All=AO+VS;

%since VS is larger then AO:
AOVS=[AO VS];
VS=[];
AO=[];

for n=1:size(AOVS,1)
    VS(n,1)=max(AOVS(n,:));
    AO(n,1)=min(AOVS(n,:));
end


% %exact zone
% aoZone = [29:38 43:52 63:72 111:120 39:42];
% vsZone = [9:18 81:90 95:104 127:136 91:94];
% 
% %exact zone +/- 1 bin
% aoZone = [28:39 42:53 62:73 110:121 38:43];
% vsZone = [8:19 80:91 94:105 126:137 90:95];
% % 
% %exact zone +/- 2 bin
% aoZone = [27:40 41:54 61:74 109:122 37:44];
% vsZone = [7:20 79:93 93:106 125:138 89:96];

%exact zone +/- 3 bin
% aoZone = [26:41 40:55 60:75 108:123 36:45];
% vsZone = [6:22 78:94 92:107 124:139 88:97];
% emptyZone=setdiff([1:1:140],unique([aoZone vsZone]));
% 
%exact zone +/- 4 bin
aoZone = [25:42 39:56 59:76 107:124 35:46];
vsZone = [5:23 77:95 91:108 123:140 87:98];
emptyZone=setdiff([1:1:140],unique([aoZone vsZone]));
% 
% %exact zone +/- 5 bin
% aoZone = [24:43 38:57 58:77 106:125 34:47];
% vsZone = [4:24 76:96 90:109 122:141 86:99];

% overlap=[]; %calculate the overlap between the zones and ao and vs zone


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% for n=1:size(all4BinPairsAll)
%     a=[all4BinPairsAll(n,1):all4BinPairsAll(n,2)];
%     b=[all4BinPairsAll(n,3):all4BinPairsAll(n,4)];
%     aao=intersect(a,aoZone);
%     avs=intersect(a,vsZone);
%     % if ~isempty(aao);
%     if length(aao)>=2;%need to overlap at lest 2 bins
%         overlap(n,1)=1;
%     % elseif ~isempty(avs);
%     elseif length(avs)>=2;
%         overlap(n,1)=2;
%     else
%         overlap(n,1)=0;
%     end
% 
%     bao=intersect(b,aoZone);
%     bvs=intersect(b,vsZone);
%     % if ~isempty(bao);
%     if length(bao)>=2;
%         overlap(n,2)=1;
%     % elseif ~isempty(bvs);
%     elseif length(bvs)>=2;
%         overlap(n,2)=2;
%     else
%         overlap(n,2)=0;
%     end
% end
% 
% overlapSum=sum(overlap,2);
% overlapDiff=overlap(:,1)-overlap(:,2); %see whether the pairs that have two same cue types or no cues at all
% iEmpty=find(overlapSum==0);
% i=find(overlapDiff==0); %same cue type or no cue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%try a different ways to classify the zones, each zone got three numbers:
%first number: overlap with auditory, second with visual, last with nothing

overlap=[];
for n=1:size(all4BinPairsAll)
    a=[all4BinPairsAll(n,1):all4BinPairsAll(n,2)];
    b=[all4BinPairsAll(n,3):all4BinPairsAll(n,4)];
    aao=intersect(a,aoZone);
    avs=intersect(a,vsZone);
    aep=intersect(a,emptyZone);

    if length(aao)>=2;
        overlap(n,1)=1;
    else
        overlap(n,1)=0;
    end

    if length(avs)>=2;
        overlap(n,2)=1;
    else
        overlap(n,2)=0;
    end

    if length(aep)>=2;
        overlap(n,3)=1;
    else
        overlap(n,3)=0;
    end


    bao=intersect(b,aoZone);
    bvs=intersect(b,vsZone);
    bep=intersect(b,emptyZone);

    if length(bao)>=2;
        overlap(n,4)=1;
    else
        overlap(n,4)=0;
    end

    if length(bvs)>=2;
        overlap(n,5)=1;
    else
        overlap(n,5)=0;
    end

    if length(bep)>=2;
        overlap(n,6)=1;
    else
        overlap(n,6)=0;
    end

end
overlapDiff=abs(overlap(:,[1:3])-overlap(:,[4:6])); 
i=find(sum(overlapDiff,2)==0);

% overlapS=overlap(:,1).*overlap(:,2);%if positive, means it covers two types of cues
% overlapS1=overlap(:,1)+overlap(:,2);
% i1=find(overlapS>0&overlapS1==2);
% 
% i=setdiff(i,i1);


figure,
subplot(141)
[p,x]=ksdensity(All(i));
plot(x,p)
real=(length(aoRewardIdx)+length(vsRewardIdx))/length(rewardIdx);
hold on
line([real real],[0 max(p)])
per=length(find(All(i)<=real))/length(i);
title(['percent ao+vs',num2str(per)])

subplot(142)
[p,x]=ksdensity(AO(i));
plot(x,p)
real=(length(aoRewardIdx))/length(rewardIdx);
hold on
line([real real],[0 max(p)])
per=length(find(AO(i)<=real))/length(i);
title(['percent ao',num2str(per)])

subplot(143)
[p,x]=ksdensity(VS(i));
plot(x,p)
real=(length(vsRewardIdx))/length(rewardIdx);
hold on
line([real real],[0 max(p)])
per=length(find(VS(i)<=real))/length(i);
title(['percent vs',num2str(per)])

subplot(144)
track=zeros(1,140);
a=unique(aoZone);
v=unique(vsZone);
track(a)=1;
track(v)=2;

ca=contiguous(track);
a1=ca{2,2};
v2=ca{3,2};

hold on
trackBaseline=zeros(1,140);
plot(trackBaseline,'k')
for n=1:size(a1,1);
    hold on
    x=[a1(n,1) a1(n,1) a1(n,2) a1(n,2)];
    y=[0 1 1 0];
    patch(x,y,'m');
end

for n=1:size(v2,1);
    hold on
    x=[v2(n,1) v2(n,1) v2(n,2) v2(n,2)];
    y=[0 1 1 0];
    patch(x,y,'g');
end

axis off
saveas(gcf,'Reward_shuffleLocations.fig')
print -painters -depsc Reward_shuffleLocations.eps





% color=[1 0 0;1 0 1;0 1 0;0 0 1;0 0 0]; %r, m, g, b, k
% colorSum=[];
% for n=1:length(i);
%     for m=1:5;
%     if overlapDiff(i(n))==m-1;
%         colorSum(n,:)=color(m,:);
%     end
%     end
% end
% 
% figure,
% for n=1:length(i);
%     hold on
%     plot(n,All(i(n)),'.','Color',colorSum(n,:));
% end