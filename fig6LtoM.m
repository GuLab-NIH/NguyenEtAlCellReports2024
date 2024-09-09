%load data: contact lead author
%%

%plot cells
otherIdx={};
figure,
for ff=1:20;

  AOcells = allROIInfos{ff}(:,:,aoIdx{ff})*3;
    VScells = allROIInfos{ff}(:,:,vsIdx{ff})*2;
otherCellsIdx=setdiff([1:1:size(allROIInfos{ff},3)],[aoIdx{ff};vsIdx{ff}]);
otherCells = allROIInfos{ff}(:,:,otherCellsIdx);
otherIdx{ff}=otherCellsIdx;
    subplot(4,5,ff)
    imagesc(sum(AOcells,3)+sum(VScells,3)+sum(otherCells,3));
   title(['fov',num2str(ff)])
end

saveas(gcf,'AOVScellAnatomicalLocations_withOther.fig');

U=[];
for n=1:20;
U(1,n)=length(aoIdx{n});
U(2,n)=length(vsIdx{n});
U(3,n)=length(otherIdx{n});
end
%%
%plot in dots

figure,
for ff=1:20;

  AOcells = allROIInfos{ff}(:,:,aoIdx{ff})*3;
    VScells = allROIInfos{ff}(:,:,vsIdx{ff})*2;
otherCellsIdx=setdiff([1:1:size(allROIInfos{ff},3)],[aoIdx{ff};vsIdx{ff}]);
otherCells = allROIInfos{ff}(:,:,otherCellsIdx);

 AOcentroids = nan(size(AOcells,3), 2);
    for nn = 1:size(AOcells,3)
        AOcentroids(nn,:) = center_Of_Mass(AOcells(:,:,nn), 0, 0);
    end
    
    VScentroids = nan(size(VScells,3), 2);
    for nn = 1:size(VScells,3)
        VScentroids(nn,:) = center_Of_Mass(VScells(:,:,nn), 0, 0);
    end


    Othercentroids = nan(size(otherCells,3), 2);
    for nn = 1:size(otherCells,3)
        Othercentroids(nn,:) = center_Of_Mass(otherCells(:,:,nn), 0, 0);
    end

    subplot(4,5,ff)

    plot(Othercentroids(:,1),Othercentroids(:,2),'.','MarkerSize',10,'Color',[0.1 0.1 0.1]);
    hold on
    plot(AOcentroids(:,1),AOcentroids(:,2),'m.','MarkerSize',10);
    hold on
    plot(VScentroids(:,1),VScentroids(:,2),'g.','MarkerSize',10);
    axis equal
set(gca,'Color','k')
end

saveas(gcf,'AOVScellAnatomicalLocations_withOther_dots.fig');
%% AO vs VS cells

% p=pwd;
% % AO vs VS cells
% load('idxAOFOV.mat'); aoIdx = idxAOFOV;
% load('idxVSFOV.mat'); vsIdx = idxVSFOV;

% a=[];%find number of cells in each FOV
% for n=1:length(aoIdx);
%     a(n)=length(aoIdx{n});
% end
% i=find(a>10);
i=1:1:length(aoIdx);

%should be at lest 10 cells (three fovs are 1 and 3 and 6 cells)

v=[];%find number of cells in each FOV
for n=1:length(vsIdx);
    v(n)=length(vsIdx{n});
end

allROIInfosUse=allROIInfos(i);
aoIdxUse=aoIdx(i);
vsIdxUse=vsIdx(i);

% % A-rew vs V-rew cells
% load('idxARFOV.mat'); aoIdx = idxARFOV;
% load('idxVRFOV.mat'); vsIdx = idxVRFOV;


ntry=10; %get distances 10 times
nShuffle=10;
ABAll={};
shuffleAll={};%each cell is one fov
for ff = 1:length(allROIInfosUse)
    disp(ff)

    a=length(aoIdxUse{ff});
    b=length(vsIdxUse{ff});

    if a*b==0;
        ABAll{ff}=nan;
        shuffleAll{ff}=nan;
    else
    % get 2 populations
    AOcells = allROIInfosUse{ff}(:,:,aoIdxUse{ff});
    VScells = allROIInfosUse{ff}(:,:,vsIdxUse{ff});
    % 
    % if size(AOcells,3) == 0 || size(VScells,3) == 0
    %     continue
    % end
    
    % centers
    C={};
    AOcentroids = nan(size(AOcells,3), 2);
    for nn = 1:size(AOcells,3)
        AOcentroids(nn,:) = center_Of_Mass(AOcells(:,:,nn), 0, 0);
    end
    
    VScentroids = nan(size(VScells,3), 2);
    for nn = 1:size(VScells,3)
        VScentroids(nn,:) = center_Of_Mass(VScells(:,:,nn), 0, 0);
    end
    
    C{1}=AOcentroids;
    C{2}=VScentroids;

    nCells(1)=size(AOcentroids,1);
    nCells(2)=size(VScentroids,1);
   
    %smallest numbers of cells
    [nUse,i]=min(nCells);

  AB=[];
shuffle=[];
i=i(1);
    for n=1:ntry;
        groupA=C{i};%this can actually be ao or vs
        ii=setdiff([1:1:2],i);
        g=C{ii};%the other group
        g=g(randperm(size(g,1)),:);
        groupB=g(1:nUse,:);

        AtoBTry=[];
        BtoATry=[];
        ABTry=[];

        shuffleTry=[];
        

        A=pdist2(groupA,groupB,'euclidean','Smallest',nUse);%each column is a cell. distances are sorterd from small to large
        B=pdist2(groupB,groupA,'euclidean','Smallest',nUse);
        %take the mean
        for m=1:nUse;
            AtoBTry(1,m)=nanmean(A(1:m,:),'all');
            BtoATry(1,m)=nanmean(B(1:m,:),'all');
        end

        ABTry=nanmean([AtoBTry;BtoATry],1);

        %shuffle
        group=[groupA;groupB];
        ABS=[];
     
        for ii=1:nShuffle;
            
            groupS=group(randperm(size(group,1)),:);
            groupAS=groupS(1:nUse,:);
            groupBS=groupS(nUse+1:end,:);
            AS=pdist2(groupAS,groupBS,'euclidean','Smallest',nUse);%each column is a cell. distances are sorterd from small to large
        BS=pdist2(groupBS,groupAS,'euclidean','Smallest',nUse);
        AtoBTryS=[];
        BtoATryS=[];
        for m=1:nUse;
            AtoBTryS(1,m)=nanmean(AS(1:m,:),'all');
            BtoATryS(1,m)=nanmean(BS(1:m,:),'all');
        end

        ABS(ii,:)=nanmean([AtoBTryS;BtoATryS],1);
        end
        shuffleTry=nanmean(ABS,1);
        %PUT THEM TO AB

        AB(n,:)=ABTry;
        shuffle(n,:)=shuffleTry;
    end

    ABAll{ff}=nanmean(AB,1);
    shuffleAll{ff}=nanmean(shuffle,1);
    end
end
        

save('ABAll.mat','ABAll');
save('shuffleAll.mat','shuffleAll');

%% plot


ABAllM=cell2mat(ABAll);
shuffleAllM=cell2mat(shuffleAll);
%change this to microns
ABAllM=ABAllM*500/512;
shuffleAllM=shuffleAllM*500/512;

ABAllM=ABAllM(~isnan(ABAllM));
shuffleAllM=shuffleAllM(~isnan(shuffleAllM));

figure,

subplot(231)
cdfplot(ABAllM)
hold on
cdfplot(shuffleAllM)

subplot(232)
ksdensity(ABAllM)
hold on
ksdensity(shuffleAllM)
[~,p]=kstest2(ABAllM',shuffleAllM');
title(['kstest2 p=',num2str(p)])

subplot(233)
for n=1:length(ABAllM)
    P=[ABAllM(n) shuffleAllM(n)];
    hold on
    plot([1 2],P);
end
[~,p]=ttest(ABAllM,shuffleAllM);
title(['ttest', num2str(p)])

subplot(234)
M=[ABAllM;shuffleAllM];
bar([1 2],nanmean(M,2));
hold on
errorbar([1 2],nanmean(M,2),nansem(M,2),'.');
[~,p]=ttest2(ABAllM,shuffleAllM);
title(['ttest2 p=',num2str(p)]);
xlim([0.5 2.5]);

%
iii=[];
for n=1:length(ABAll);
    iii(n)=length(ABAll{n});
end

ABAllMM=nan(max(iii),length(ABAll));
for n=1:length(ABAll);
    ABAllMM(1:length(ABAll{n}),n)=ABAll{n};
end

shuffleAllMM=nan(max(iii),length(shuffleAll));
for n=1:length(shuffleAll);
    shuffleAllMM(1:length(shuffleAll{n}),n)=shuffleAll{n};
end

%plot by n cell
subplot(235)

nonan=[];
for n=1:size(ABAllMM,1);
    nonan(n)=length(find(~isnan(ABAllMM(n,:))));
end

iiii=find(nonan>1);
nSample=max(iiii);
errorbar([1:1:nSample]+5,nanmean(ABAllMM(1:nSample,:),2),nansem(ABAllMM(1:nSample,:),2))
hold on
errorbar([1:1:nSample]+5,nanmean(shuffleAllMM(1:nSample,:),2),nansem(shuffleAllMM(1:nSample,:),2))
saveas(gcf,'distances.fig')