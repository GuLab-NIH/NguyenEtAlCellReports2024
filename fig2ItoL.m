
%% load aligment directories
cd ..\
load('anatomicalDistance_2.mat');
cd('NNNearNeighbor_uni_uniqueCells')



%% CELL ANATOMY - AO1 vs VS1

folders = folders_audio_vs_visual;
% 
% distancesInfos = nan( length(folders), 3 );      % AO1 vs VS1
% AOtoAO
% AOtoVS
% VStoAO
% VStoVS

% allROIs = cell( length(folders), 2 );
p=pwd;
ntry=10; %get distances 10 times
nShuffle=10;
ABAll={};
shuffleAll={};%each cell is one fov

for ff = 1:length(folders)
    disp(ff)
    dd = dir([folders{ff} '\1\1_2\cellRegistered*']);
    load([folders{ff} '\1\1_2\' dd(1).name]);
    cellIdx = cell_registered_struct.cell_to_index_map;
    
    % remove non-cell
    load([folders{ff} '\visualVR\pcaica\roiIdxUse.mat']);
    roiIdxUse = find(roiIdxUse);
    cellIdx( ~ismember(cellIdx(:,2), roiIdxUse), 2) = 0;
    
    
    % find unique-cell populations
    cc = find(prod(cellIdx, 2) == 0);
    % cc = find(sum(cellIdx, 2) > 0);
    Uao = cellIdx(cc, 1); Uao = unique(Uao); Uao = Uao(2:end);
    Uvs = cellIdx(cc, 2); Uvs = unique(Uvs); Uvs = Uvs(2:end);
    
    % centers
    AOcentroids = cell_registered_struct.centroid_locations_corrected{1};
    AOcentroids = AOcentroids(Uao,:);
    VScentroids = cell_registered_struct.centroid_locations_corrected{2};
    VScentroids = VScentroids(Uvs,:);
    
    allROIs{ff,1} = AOcentroids;
    allROIs{ff,2} = VScentroids;
    
    C{1} = AOcentroids;
    C{2} = VScentroids;

    nCells(1) = size(AOcentroids,1);
    nCells(2) = size(VScentroids,1);

   %smallest numbers of cells
    [nUse,i]=min(nCells);
    
 AB=[];
shuffle=[];
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
        for m=5:nUse;
            AtoBTry(1,m-4)=nanmean(A(1:m,:),'all');
            BtoATry(1,m-4)=nanmean(B(1:m,:),'all');
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
        for m=5:nUse;
            AtoBTryS(1,m-4)=nanmean(AS(1:m,:),'all');
            BtoATryS(1,m-4)=nanmean(BS(1:m,:),'all');
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
cd(p)
save('avr1vvr1.mat','ABAll','shuffleAll');
%% plot


ABAllM=cell2mat(ABAll);
shuffleAllM=cell2mat(shuffleAll);
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
saveas(gcf,'distances_avr1vvr1.fig')




%% CELL ANATOMY - AO2 vs VS2
% 
% distancesInfos_2 = nan(60, 3);      % AO2 vs VS2
% % AOtoAO
% % AOtoVS
% % VStoAO
% % VStoVS
% allROIs_2 = cell(60, 2);



p=pwd;
ntry=10; %get distances 10 times
nShuffle=10;
ABAll={};
shuffleAll={};%each cell is one fov

cd('Z:\labMembers\DN\_PROJECT_\finalData\allAlignments\audio_NEW_vs_visual_NEW');
dd1 = dir('ID202*');
ff = 1;
for nMice = 1:numel(dd1)
   
    cd( dd1(nMice).name );
    dd2 = dir('loc*');
    for nFov = 1:numel(dd2)
        cd( dd2(nFov).name );
        load('folders.mat');
        
        NvsO = contains( folders(1:end), '\audioVR_NEW\');      % identify which directory correspond to audio and visual environment
        ao = find(NvsO == 1);
        vs = find(NvsO == 0);

        
        for idx1 = 1:numel(ao)
        for idx2 = 1:numel(vs)
            disp(ff)
            env1 = ao(idx1); env2 = vs(idx2);

            % load common cells
            fl = dir([num2str(env1) '\' num2str(env1) '_' num2str(env2) '\cellRegistered*.mat']);
            load([num2str(env1) '\' num2str(env1) '_' num2str(env2) '\' fl(1).name]);
            cellIdx = cell_registered_struct.cell_to_index_map;
            
            % remove non-cell
            load([folders{env2} '\roiIdxUse.mat']);
            roiIdxUse = find(roiIdxUse);
            cellIdx( ~ismember(cellIdx(:,2), roiIdxUse), 2) = 0;
            
            % find unique-cell populations
            cc = find(prod(cellIdx, 2) == 0);
            % cc = find(sum(cellIdx, 2) > 0);
            Uao = cellIdx(cc, 1); Uao = unique(Uao); Uao = Uao(2:end);
            Uvs = cellIdx(cc, 2); Uvs = unique(Uvs); Uvs = Uvs(2:end);

            % centers
            AOcentroids = cell_registered_struct.centroid_locations_corrected{1};
            AOcentroids = AOcentroids(Uao,:);
            VScentroids = cell_registered_struct.centroid_locations_corrected{2};
            VScentroids = VScentroids(Uvs,:);
            
            allROIs_2{ff, 1} = AOcentroids;
            allROIs_2{ff, 2} = VScentroids;
            
    C{1}=AOcentroids;
    C{2}=VScentroids;

    nCells(1)=size(AOcentroids,1);
    nCells(2)=size(VScentroids,1);
   
    %smallest numbers of cells
    [nUse,i]=min(nCells);

  AB=[];
shuffle=[];
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
        for m=5:nUse;
            AtoBTry(1,m-4)=nanmean(A(1:m,:),'all');
            BtoATry(1,m-4)=nanmean(B(1:m,:),'all');
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
        for m=5:nUse;
            AtoBTryS(1,m-4)=nanmean(AS(1:m,:),'all');
            BtoATryS(1,m-4)=nanmean(BS(1:m,:),'all');
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
             
            
            ff = ff + 1;
        end
        end
        cd ..\
    end
    cd ..\
end
cd(p)
save('avr2vvr2.mat','ABAll','shuffleAll');
%% plot


ABAllM=cell2mat(ABAll);
shuffleAllM=cell2mat(shuffleAll);
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
saveas(gcf,'distances_avr2vvr2.fig')



% %% combine both
% 
% allDist = [];
% allDist = [allDist; distancesInfos(:,[1 3 2])];
% allDist = [allDist; distancesInfos_2(:,[1 3 2])];
% 
% figure; hold on
% for n=1:size(allDist,1)
%     plot(allDist(n,:), 'c--', 'LineWidth', 0.2);
% end
% 
% errorbar(1, nanmean(allDist(:,1)), nanstd(allDist(:,1))/sqrt(size(allDist,1)),...
%     'm', 'LineWidth', 3);
% errorbar(2, nanmean(allDist(:,2)), nanstd(allDist(:,2))/sqrt(size(allDist,1)),...
%     'k', 'LineWidth', 3);
% errorbar(3, nanmean(allDist(:,3)), nanstd(allDist(:,3))/sqrt(size(allDist,1)),...
%     'g', 'LineWidth', 3);
% xlim([0.8, 3.2]); ylabel('Average Distance (\mum)');
% set(gca, 'xtick', 1:3, 'xticklabel', {'ao-ao', 'ao-vs', 'vs-vs'});
% set(gcf, 'Position', [222 222 345 333]);
% title('AO vs VS');
% 
% 
% [~, p1] = ttest2(allDist(:,1), allDist(:,2)); disp(p1);
% [~, p2] = ttest2(allDist(:,1), allDist(:,3)); disp(p2);
% [~, p3] = ttest2(allDist(:,2), allDist(:,3)); disp(p3);
% 
% % [isSignificant,adjusted_pvals,~]= bonferroni_holm([p1 p2 p3],0.05) % adjusted_pvals are the corrected p values
% 
% 
% % adjusted_pvals =
% % 
%     0.0000    0.1014    0.0002

%% CELL ANATOMY - AO1 vs AO2
% 
% distancesInfos_3 = nan(111, 3);      % AO2 vs VS2
% allROIs_3 = cell(111,2);

p=pwd;
    
ntry=10; %get distances 10 times
nShuffle=10;
ABAll={};
shuffleAll={};%each cell is one fov


cd('Z:\labMembers\DN\_PROJECT_\finalData\allAlignments\audio_old_vs_new_afterLearned');
dd1 = dir('ID202*');
ff = 1;
for nMice = 1:numel(dd1)
    cd( dd1(nMice).name );
    dd2 = dir('loc*');
    for nFov = 1:numel(dd2)
        if nMice == 2 && nFov == 2
            continue
        end
        cd( dd2(nFov).name );
        load('folders.mat');
        
        NvsO = contains( folders(2:end), '\audioVR_NEW\');      % identify which directory correspond to OLD vs NEW environment
        AOn = find(NvsO, 2, 'last') + 1;
        AOo = find(NvsO == 0) + 1;
        % AOo = find(NvsO == 0, 2, 'last') + 1;

        for idx1 = 1:numel(AOn)
        for idx2 = 1:numel(AOo)
            disp(ff)
            env1 = AOn(idx1); env2 = AOo(idx2);
            
            % load common cells
            fl = dir([num2str(env1) '\' num2str(env1) '_' num2str(env2) '\cellRegistered*.mat']);
            load([num2str(env1) '\' num2str(env1) '_' num2str(env2) '\' fl(1).name]);
            cellIdx = cell_registered_struct.cell_to_index_map;
            
            % find unique-cell populations
            cc = find(prod(cellIdx, 2) == 0);
            % cc = find(sum(cellIdx, 2) > 0);
            Uao2 = cellIdx(cc, 1); Uao2 = unique(Uao2); Uao2 = Uao2(2:end);
            Uao1 = cellIdx(cc, 2); Uao1 = unique(Uao1); Uao1 = Uao1(2:end);

            % centers
            AO2centroids = cell_registered_struct.centroid_locations_corrected{1};
            AO2centroids = AO2centroids(Uao2,:);
            AO1centroids = cell_registered_struct.centroid_locations_corrected{2};
            AO1centroids = AO1centroids(Uao1,:);
            
            allROIs_3{ff, 1} = AO2centroids;
            allROIs_3{ff, 2} = AO1centroids;
            
                C{1}=AO1centroids;
    C{2}=AO2centroids;

    nCells(1)=size(AO1centroids,1);
    nCells(2)=size(AO2centroids,1);
   
    %smallest numbers of cells
    [nUse,i]=min(nCells);

  AB=[];
shuffle=[];
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
        for m=5:nUse;
            AtoBTry(1,m-4)=nanmean(A(1:m,:),'all');
            BtoATry(1,m-4)=nanmean(B(1:m,:),'all');
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
        for m=5:nUse;
            AtoBTryS(1,m-4)=nanmean(AS(1:m,:),'all');
            BtoATryS(1,m-4)=nanmean(BS(1:m,:),'all');
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
             
            
            ff = ff + 1;
        end
        end
        cd ..\
    end
    cd ..\
end
cd(p)
save('avr1avr2.mat','ABAll','shuffleAll');

%% plot


ABAllM=cell2mat(ABAll);
shuffleAllM=cell2mat(shuffleAll);
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

saveas(gcf,'distances_avr1avr2.fig')


%% CELL ANATOMY - VS1 vs VS2

% distancesInfos_4 = nan(111, 4);
% allROIs_4 = cell(111,2);
p=pwd;
    
ntry=10; %get distances 10 times
nShuffle=10;
ABAll={};
shuffleAll={};%each cell is one fov


cd('Z:\labMembers\DN\_PROJECT_\finalData\allAlignments\visual_old_vs_new_afterLearned');
dd1 = dir('ID202*');
ff = 1;
for nMice = 1:numel(dd1)
    cd( dd1(nMice).name );
    dd2 = dir('loc*');
    for nFov = 1:numel(dd2)
        cd( dd2(nFov).name );
        load('folders.mat');
        
        NvsO = contains( folders(1:end), '\visualVR_NEW\');      % identify which directory correspond to OLD vs NEW environment
        VSn = find(NvsO, 2, 'last');
        VSo = find(NvsO == 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%    %yg
%         VSnfolder=folders(VSn);
%         VSofolder=folders(VSo);
%         thisn=[];
%         for i=1:length(VSnfolder);
%             load([VSnfolder{i} '\roiIdxUse.mat']);
%             thisn(end+1)=length(find(roiIdxUse));
%              nCells{nMice}(nFov,i)=length(find(roiIdxUse));
%         end
%         thisn=min(thisn);
%         thiso=[];
%         for i=1:length(VSofolder);
%             load([VSofolder{i} '\roiIdxUse.mat']);
%             thiso(end+1)=length(find(roiIdxUse));
%              nCells{nMice}(nFov,i+2)=length(find(roiIdxUse));
%         end
%  thiso=max(thiso);
%         diffno=thiso-thisn;
%         if diffno>50;
%             cd ..\
%         continue
%         end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for idx1 = 1:numel(VSn)
        for idx2 = 1:numel(VSo)
            disp(ff)
            env1 = VSn(idx1); env2 = VSo(idx2);
            clear fields_ref dfofM_ref SF_ref cell_registered_struct commonCells;
            
            load([folders{env1} '\roiIdxUse.mat']);
            roiIdxUse_ref = find(roiIdxUse); clear roiIdxUse;
            load([folders{env2} '\roiIdxUse.mat']);
            roiIdxUse_target = find(roiIdxUse); clear roiIdxUse;
            
            try
                % load common cells
                fl = dir([num2str(env2) '\' num2str(env2) '_' num2str(env1) '\cellRegistered*.mat']);
                load([num2str(env2) '\' num2str(env2) '_' num2str(env1) '\' fl(1).name]);
                cellIdx = cell_registered_struct.cell_to_index_map;
                cellIdx( ~ismember(cellIdx(:,2), roiIdxUse_ref), 2) = 0;
                cellIdx( ~ismember(cellIdx(:,1), roiIdxUse_target), 1) = 0;
            catch er
                continue
            end
            
            % find unique-cell populations
            cc = find(prod(cellIdx, 2) == 0);
            % cc = find(sum(cellIdx, 2) > 0);
            Uvs1 = cellIdx(cc, 1); Uvs1 = unique(Uvs1); Uvs1 = Uvs1(2:end);
            Uvs2 = cellIdx(cc, 2); Uvs2 = unique(Uvs2); Uvs2 = Uvs2(2:end);

            % centers
            VS1centroids = cell_registered_struct.centroid_locations_corrected{1};
            VS1centroids = VS1centroids(Uvs1,:);
            VS2centroids = cell_registered_struct.centroid_locations_corrected{2};
            VS2centroids = VS2centroids(Uvs2,:);
            
            allROIs_4{ff, 1} = VS1centroids;
            allROIs_4{ff, 2} = VS2centroids;
            
                            C{1}=VS1centroids;
    C{2}=VS2centroids;

    nCells(1)=size(VS1centroids,1);
    nCells(2)=size(VS2centroids,1);
   
    %smallest numbers of cells
    [nUse,i]=min(nCells);

  AB=[];
shuffle=[];
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
        for m=5:nUse;
            AtoBTry(1,m-4)=nanmean(A(1:m,:),'all');
            BtoATry(1,m-4)=nanmean(B(1:m,:),'all');
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
        for m=5:nUse;
            AtoBTryS(1,m-4)=nanmean(AS(1:m,:),'all');
            BtoATryS(1,m-4)=nanmean(BS(1:m,:),'all');
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
             
            
            ff = ff + 1;
        end
        end
        cd ..\
    end
    cd ..\
end
cd(p)
save('vvr1vvr2.mat','ABAll','shuffleAll');

%% plot

ABAllM=cell2mat(ABAll);
shuffleAllM=cell2mat(shuffleAll);
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

saveas(gcf,'distances_vvr1vvr2.fig')


%% CELL ANATOMY - AO1 vs AO2 - only AO1 in AO1-VS1


p=pwd;
    
ntry=10; %get distances 10 times
nShuffle=10;
ABAll={};
shuffleAll={};%each cell is one fov


cd('Z:\labMembers\DN\_PROJECT_\finalData\allAlignments\audio_old_vs_new_afterLearned');
dd1 = dir('ID202*');
ff = 1;

%yg
% nCells={};
useIdx=[];

for nMice = 1:numel(dd1)
    cd( dd1(nMice).name );
    dd2 = dir('loc*');
    for nFov = 1:numel(dd2)
        if nMice == 2 && nFov == 2
            continue
        end
        cd( dd2(nFov).name );
        load('folders.mat');
        
        NvsO = contains( folders(2:end), '\audioVR_NEW\');      % identify which directory correspond to OLD vs NEW environment
        AOn = find(NvsO, 2, 'last') + 1;
        AOo = find(NvsO == 0, 2, 'last') + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %yg
%         AOnfolder=folders(AOn);
%         AOofolder=folders(AOo);
%         for i=1:length(AOnfolder);
%             load([AOnfolder{i} '\allROIs.mat']);
%             nCells{nMice}(nFov,i)=size(roi,3);
%         end
%         thisn=min(nCells{nMice}(nFov,[1 2]));
%         for i=1:length(AOofolder);
%             load([AOofolder{i} '\allROIs.mat']);
%             nCells{nMice}(nFov,i+2)=size(roi,3);
%         end
%  thiso=max(nCells{nMice}(nFov,[3 4]));
%         diffno=thiso-thisn;
%         if diffno>60;
%             cd ..\
%         continue
%         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cd ..\
%     end
%     cd ..\
% end
% 
% diffno=[];
% for i=1:length(nCells);
%     for ii=1:size(nCells{i},1);
%         diffno(end+1)=max(nCells{i}(ii,3:4))-min(nCells{i}(ii,1:2));
%     end
% end
% 
% 
% stdno=[];
% for i=1:length(nCells);
%     for ii=1:size(nCells{i},1);
%         stdno(end+1)=std(nCells{i}(ii,:),0,2);
%     end
% end



       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        for idx1 = 1:numel(AOn)
        for idx2 = 1:numel(AOo)
            disp(ff)
            env1 = AOn(idx1); env2 = AOo(idx2);
            useIdx(end+1,:)=[nMice nFov idx1 idx2];
            % load common cells
            fl = dir([num2str(env1) '\' num2str(env1) '_' num2str(env2) '\cellRegistered*.mat']);
            load([num2str(env1) '\' num2str(env1) '_' num2str(env2) '\' fl(1).name]);
            cellIdx = cell_registered_struct.cell_to_index_map;
            
            % find unique-cell populations
            cc = find(prod(cellIdx, 2) == 0);
            % cc = find(sum(cellIdx, 2) > 0);
            Uao2 = cellIdx(cc, 1); Uao2 = unique(Uao2); Uao2 = Uao2(2:end);
            Uao1 = cellIdx(cc, 2); Uao1 = unique(Uao1); Uao1 = Uao1(2:end);

            % centers
            AO2centroids = cell_registered_struct.centroid_locations_corrected{1};
            AO2centroids = AO2centroids(Uao2,:);
            AO1centroids = cell_registered_struct.centroid_locations_corrected{2};
            AO1centroids = AO1centroids(Uao1,:);

            
            allROIs_3{ff, 1} = AO2centroids;
            allROIs_3{ff, 2} = AO1centroids;
            
                C{1}=AO1centroids;
    C{2}=AO2centroids;

    nCells(1)=size(AO1centroids,1);
    nCells(2)=size(AO2centroids,1);
   
    %smallest numbers of cells
    [nUse,i]=min(nCells);

  AB=[];
shuffle=[];
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
        for m=5:nUse;
            AtoBTry(1,m-4)=nanmean(A(1:m,:),'all');
            BtoATry(1,m-4)=nanmean(B(1:m,:),'all');
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
        for m=5:nUse;
            AtoBTryS(1,m-4)=nanmean(AS(1:m,:),'all');
            BtoATryS(1,m-4)=nanmean(BS(1:m,:),'all');
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
            
            ff = ff + 1;
        end
        end
        
        cd ..\
       
    end
    cd ..\
end
cd(p)

save('avr1avr2_withVVR.mat','ABAll','shuffleAll');

%% plot


ABAllM=cell2mat(ABAll);
shuffleAllM=cell2mat(shuffleAll);
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

saveas(gcf,'distances_avr1avr2WithVVR.fig')

