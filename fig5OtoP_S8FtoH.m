
%% load aligment directories: contact lead author for data



%% CELL ANATOMY - AO1 vs VS1

% folders = folders_audio_vs_visual;
% 
% distancesInfos = nan( length(folders), 3 );      % AO1 vs VS1
% AOtoAO
% AOtoVS
% VStoAO
% VStoVS

% allROIs = cell( length(folders), 2 );

if 0
    ntry = 100;
for ff = 1:length(folders)
    
    dd = dir([folders{ff} '\1\1_2\cellRegistered*']);
    load([folders{ff} '\1\1_2\' dd(1).name]);
    cellIdx = cell_registered_struct.cell_to_index_map;
    
    % remove non-cell
    load([folders{ff} '\visualVR\pcaica\roiIdxUse.mat']);
    roiIdxUse = find(roiIdxUse);
    cellIdx( ~ismember(cellIdx(:,2), roiIdxUse), 2) = 0;
    
    
    % find unique-cell populations
%     cc = find(prod(cellIdx, 2) == 0);
    cc = find(sum(cellIdx, 2) > 0);
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
    [nUse, iidx] = min(nCells);
    
    for nn = 1:length(C)
        
        if nn == iidx   % if the cell group is the one with smallerst cell number
            Dis = pdist2(C{nn},C{nn}, 'euclidean');
            Dis = unique(Dis);
            Dis = Dis(2:end);%remove zero
            distancesInfos(ff,nn) = nanmean(Dis);
            
        else    %if the cell group is the one with a larger cell number
            
            CAll = C{nn};
            disSelfTry = nan(1, ntry);
            disOtherTry = nan(1, ntry);
            
            for t = 1:ntry
                CAllTry=CAll(randperm(size(CAll,1)),:);
                CAllTry=CAllTry(1:nUse,:);%random nUse cells

                %self distances
                Dis = pdist2(CAllTry,CAllTry, 'euclidean');
                Dis = unique(Dis);
                Dis = Dis(2:end);%remove zero
                disSelfTry(t) = nanmean(Dis);

                %distances across groups
                otherGroup = setdiff([1 2],nn);
                Dis = pdist2(CAllTry, C{otherGroup}, 'euclidean');
                disOtherTry(t) = nanmean(Dis,'all');%this one no duplicate values
            end

            distancesInfos(ff,nn) = nanmean(disSelfTry);
            distancesInfos(ff,3) = nanmean(disOtherTry);
        end
    end
end
end



% PLOTTING

allDist=distancesInfos(:,[1 3 2]);%first column is auditory, second is between, the last is visual.

s=sum(allDist,2);
allDist=allDist(~isnan(s),:); %getting rid of the ones with nan

figure; 
for n=1:size(allDist,1)
    hold on
    plot(allDist(n,:), 'c--', 'LineWidth', 0.2);
end

errorbar(1, nanmean(allDist(:,1)), nanstd(allDist(:,1))/sqrt(size(allDist,1)),...
    'm', 'LineWidth', 3);
errorbar(2, nanmean(allDist(:,2)), nanstd(allDist(:,2))/sqrt(size(allDist,1)),...
    'k', 'LineWidth', 3);
errorbar(3, nanmean(allDist(:,3)), nanstd(allDist(:,3))/sqrt(size(allDist,1)),...
    'g', 'LineWidth', 3);
xlim([0.8, 3.2]); ylabel('Average Distance (\mum)');
set(gca, 'xtick', 1:3, 'xticklabel', {'ao-ao', 'ao-vs', 'vs-vs'});
set(gcf, 'Position', [222 222 345 333]);


[~, p1] = ttest(allDist(:,1), allDist(:,2)); disp(p1);
[~, p2] = ttest(allDist(:,1), allDist(:,3)); disp(p2);
[~, p3] = ttest(allDist(:,2), allDist(:,3)); disp(p3);


[isSignificant,adjusted_pvals,~]= bonferroni_holm([p1 p2 p3],0.05); % adjusted_pvals are the corrected p values
disp(adjusted_pvals);

title(['12:',num2str(adjusted_pvals(1)),' 13:',num2str(adjusted_pvals(2)),' 23:',num2str(adjusted_pvals(3))])

figure; 
for ff = 1:size(allROIs,1)
    subplot(4,8,ff); hold on
    xlim([0 512]); ylim([0 512]);
    plot(allROIs{ff,1}(:,1), allROIs{ff,1}(:,2), 'm.');
    plot(allROIs{ff,2}(:,1), allROIs{ff,2}(:,2), 'g.');
end


%% CELL ANATOMY - AO2 vs VS2
% 
% distancesInfos_2 = nan(60, 3);      % AO2 vs VS2
% % AOtoAO
% % AOtoVS
% % VStoAO
% % VStoVS
% allROIs_2 = cell(60, 2);

if 0
tic;

ntry = 100;
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
%             cc = find(prod(cellIdx, 2) == 0);
            cc = find(sum(cellIdx, 2) > 0);
            Uao = cellIdx(cc, 1); Uao = unique(Uao); Uao = Uao(2:end);
            Uvs = cellIdx(cc, 2); Uvs = unique(Uvs); Uvs = Uvs(2:end);

            % centers
            AOcentroids = cell_registered_struct.centroid_locations_corrected{1};
            AOcentroids = AOcentroids(Uao,:);
            VScentroids = cell_registered_struct.centroid_locations_corrected{2};
            VScentroids = VScentroids(Uvs,:);
            
            allROIs_2{ff, 1} = AOcentroids;
            allROIs_2{ff, 2} = VScentroids;
            
            C{1} = AOcentroids;
            C{2} = VScentroids;

            nCells(1) = size(AOcentroids,1);
            nCells(2) = size(VScentroids,1);

            %smallest numbers of cells
            [nUse, iidx] = min(nCells);

            for nn=1:length(C)

                if nn == iidx   % if the cell group is the one with smallerst cell number
                    Dis = pdist2(C{nn},C{nn}, 'euclidean');
                    Dis = unique(Dis);
                    Dis = Dis(2:end);%remove zero
                    distancesInfos_2(ff,nn) = nanmean(Dis);

                else    %if the cell group is the one with a larger cell number

                    CAll = C{nn};
                    disSelfTry = nan(1, ntry);
                    disOtherTry = nan(1, ntry);

                    for t = 1:ntry
                        CAllTry=CAll(randperm(size(CAll,1)),:);
                        CAllTry=CAllTry(1:nUse,:);%random nUse cells

                        %self distances
                        Dis = pdist2(CAllTry,CAllTry, 'euclidean');
                        Dis = unique(Dis);
                        Dis = Dis(2:end);%remove zero
                        disSelfTry(t) = nanmean(Dis);

                        %distances across groups
                        otherGroup = setdiff([1 2],nn);
                        Dis = pdist2(CAllTry, C{otherGroup}, 'euclidean');
                        disOtherTry(t) = nanmean(Dis,'all');%this one no duplicate values
                    end

                    distancesInfos_2(ff,nn) = nanmean(disSelfTry);
                    distancesInfos_2(ff,3) = nanmean(disOtherTry);
                end
            end
            
            ff = ff + 1;
        end
        end
        cd ..\
    end
    cd ..\
end
cd ..\
toc;
end


% PLOTTING


allDist = distancesInfos_2(:,[1 3 2]);%first column is auditory, second is between, the last is visual.

i=[43 44];
i=setdiff([1:1:size(allDist,1)],i);
allDist=allDist(i,:);

s=sum(allDist,2);
allDist=allDist(~isnan(s),:); %getting rid of the ones with nan

figure; hold on
for n=1:size(allDist,1)
    plot(allDist(n,:), 'c--', 'LineWidth', 0.2);
end

errorbar(1, nanmean(allDist(:,1)), nanstd(allDist(:,1))/sqrt(size(allDist,1)),...
    'm', 'LineWidth', 3);
errorbar(2, nanmean(allDist(:,2)), nanstd(allDist(:,2))/sqrt(size(allDist,1)),...
    'k', 'LineWidth', 3);
errorbar(3, nanmean(allDist(:,3)), nanstd(allDist(:,3))/sqrt(size(allDist,1)),...
    'g', 'LineWidth', 3);
xlim([0.8, 3.2]); ylabel('Average Distance (\mum)');
set(gca, 'xtick', 1:3, 'xticklabel', {'ao-ao', 'ao-vs', 'vs-vs'});
set(gcf, 'Position', [222 222 345 333]);


[~, p1] = ttest(allDist(:,1), allDist(:,2)); disp(p1);
[~, p2] = ttest(allDist(:,1), allDist(:,3)); disp(p2);
[~, p3] = ttest(allDist(:,2), allDist(:,3)); disp(p3);


[isSignificant,adjusted_pvals,~]= bonferroni_holm([p1 p2 p3],0.05); % adjusted_pvals are the corrected p values
disp(adjusted_pvals);

title(['12:',num2str(adjusted_pvals(1)),' 13:',num2str(adjusted_pvals(2)),' 23:',num2str(adjusted_pvals(3))])



figure; 
for ff = 1:size(allROIs_2,1)
    subplot(6,10,ff); hold on
    xlim([0 512]); ylim([0 512]);
    plot(allROIs_2{ff,1}(:,1), allROIs_2{ff,1}(:,2), 'm.');
    plot(allROIs_2{ff,2}(:,1), allROIs_2{ff,2}(:,2), 'g.');
end



%% combine both

allDist = [];
allDist = [allDist; distancesInfos(:,[1 3 2])];
allDist = [allDist; distancesInfos_2(:,[1 3 2])];

figure; hold on
for n=1:size(allDist,1)
    plot(allDist(n,:), 'c--', 'LineWidth', 0.2);
end

errorbar(1, nanmean(allDist(:,1)), nanstd(allDist(:,1))/sqrt(size(allDist,1)),...
    'm', 'LineWidth', 3);
errorbar(2, nanmean(allDist(:,2)), nanstd(allDist(:,2))/sqrt(size(allDist,1)),...
    'k', 'LineWidth', 3);
errorbar(3, nanmean(allDist(:,3)), nanstd(allDist(:,3))/sqrt(size(allDist,1)),...
    'g', 'LineWidth', 3);
xlim([0.8, 3.2]); ylabel('Average Distance (\mum)');
set(gca, 'xtick', 1:3, 'xticklabel', {'ao-ao', 'ao-vs', 'vs-vs'});
set(gcf, 'Position', [222 222 345 333]);
title('AO vs VS');


[~, p1] = ttest2(allDist(:,1), allDist(:,2)); disp(p1);
[~, p2] = ttest2(allDist(:,1), allDist(:,3)); disp(p2);
[~, p3] = ttest2(allDist(:,2), allDist(:,3)); disp(p3);

% [isSignificant,adjusted_pvals,~]= bonferroni_holm([p1 p2 p3],0.05) % adjusted_pvals are the corrected p values


% adjusted_pvals =
% 
%     0.0000    0.1014    0.0002

%% CELL ANATOMY - AO1 vs AO2
% 
% distancesInfos_3 = nan(111, 3);      % AO2 vs VS2
% allROIs_3 = cell(111,2);

if 0
    
ntry = 100;
tic;
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
            env1 = AOn(idx1); env2 = AOo(idx2);
            
            % load common cells
            fl = dir([num2str(env1) '\' num2str(env1) '_' num2str(env2) '\cellRegistered*.mat']);
            load([num2str(env1) '\' num2str(env1) '_' num2str(env2) '\' fl(1).name]);
            cellIdx = cell_registered_struct.cell_to_index_map;
            
            % find unique-cell populations
%             cc = find(prod(cellIdx, 2) == 0);
            cc = find(sum(cellIdx, 2) > 0);
            Uao2 = cellIdx(cc, 1); Uao2 = unique(Uao2); Uao2 = Uao2(2:end);
            Uao1 = cellIdx(cc, 2); Uao1 = unique(Uao1); Uao1 = Uao1(2:end);

            % centers
            AO2centroids = cell_registered_struct.centroid_locations_corrected{1};
            AO2centroids = AO2centroids(Uao2,:);
            AO1centroids = cell_registered_struct.centroid_locations_corrected{2};
            AO1centroids = AO1centroids(Uao1,:);
            
            allROIs_3{ff, 1} = AO2centroids;
            allROIs_3{ff, 2} = AO1centroids;
            
            C{1} = AO2centroids;
            C{2} = AO1centroids;

            nCells(1) = size(AO2centroids,1);
            nCells(2) = size(AO1centroids,1);

            %smallest numbers of cells
            [nUse, iidx] = min(nCells);

            for nn=1:length(C);

                if nn == iidx   % if the cell group is the one with smallerst cell number
                    Dis = pdist2(C{nn},C{nn}, 'euclidean');
                    Dis = unique(Dis);
                    Dis = Dis(2:end);%remove zero
                    distancesInfos_3(ff,nn) = nanmean(Dis);

                else    %if the cell group is the one with a larger cell number

                    CAll = C{nn};
                    disSelfTry = nan(1, ntry);
                    disOtherTry = nan(1, ntry);

                    for t = 1:ntry
                        CAllTry=CAll(randperm(size(CAll,1)),:);
                        CAllTry=CAllTry(1:nUse,:);%random nUse cells

                        %self distances
                        Dis = pdist2(CAllTry,CAllTry, 'euclidean');
                        Dis = unique(Dis);
                        Dis = Dis(2:end);%remove zero
                        disSelfTry(t) = nanmean(Dis);

                        %distances across groups
                        otherGroup = setdiff([1 2],nn);
                        Dis = pdist2(CAllTry, C{otherGroup}, 'euclidean');
                        disOtherTry(t) = nanmean(Dis,'all');%this one no duplicate values
                    end

                    distancesInfos_3(ff,nn) = nanmean(disSelfTry);
                    distancesInfos_3(ff,3) = nanmean(disOtherTry);
                end
            end
            
            
            ff = ff + 1;
        end
        end
        cd ..\
    end
    cd ..\
end
cd ..\
toc;
end


% PLOTTING

allDist = distancesInfos_3(:,[1 3 2]);%first column is auditory, second is between, the last is visual.

s=sum(allDist,2);
allDist=allDist(~isnan(s),:); %getting rid of the ones with nan

figure; hold on
for n=1:size(allDist,1)
    plot(allDist(n,:), 'c--', 'LineWidth', 0.2);
end

errorbar(1, nanmean(allDist(:,1)), nanstd(allDist(:,1))/sqrt(size(allDist,1)),...
    'm', 'LineWidth', 3);
errorbar(2, nanmean(allDist(:,2)), nanstd(allDist(:,2))/sqrt(size(allDist,1)),...
    'k', 'LineWidth', 3);
errorbar(3, nanmean(allDist(:,3)), nanstd(allDist(:,3))/sqrt(size(allDist,1)),...
    'g', 'LineWidth', 3);
xlim([0.8, 3.2]); ylabel('Average Distance (\mum)');
set(gca, 'xtick', 1:3, 'xticklabel', {'ao2-ao2', 'ao1-ao2', 'ao1-ao1'});
set(gcf, 'Position', [222 222 345 333]);


[~, p1] = ttest(allDist(:,1), allDist(:,2)); disp(p1);
[~, p2] = ttest(allDist(:,1), allDist(:,3)); disp(p2);
[~, p3] = ttest(allDist(:,2), allDist(:,3)); disp(p3);

% [isSignificant,adjusted_pvals,~]= bonferroni_holm([p1 p2 p3],0.05); % adjusted_pvals are the corrected p values
% disp(adjusted_pvals);
% 
% title(['12:',num2str(adjusted_pvals(1)),' 13:',num2str(adjusted_pvals(2)),' 23:',num2str(adjusted_pvals(3))])


figure;
for ff = 1:size(allROIs_3,1)
    subplot(6,10,ff); hold on
    xlim([0 512]); ylim([0 512]);
    plot(allROIs_3{ff,1}(:,1), allROIs_3{ff,1}(:,2), 'm.');
    plot(allROIs_3{ff,2}(:,1), allROIs_3{ff,2}(:,2), 'c.');
end


%% CELL ANATOMY - VS1 vs VS2

% distancesInfos_4 = nan(111, 4);
% allROIs_4 = cell(111,2);

if 0
ntry = 100;
tic;
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
%             cc = find(prod(cellIdx, 2) == 0);
            cc = find(sum(cellIdx, 2) > 0);
            Uvs1 = cellIdx(cc, 1); Uvs1 = unique(Uvs1); Uvs1 = Uvs1(2:end);
            Uvs2 = cellIdx(cc, 2); Uvs2 = unique(Uvs2); Uvs2 = Uvs2(2:end);

            % centers
            VS1centroids = cell_registered_struct.centroid_locations_corrected{1};
            VS1centroids = VS1centroids(Uvs1,:);
            VS2centroids = cell_registered_struct.centroid_locations_corrected{2};
            VS2centroids = VS2centroids(Uvs2,:);
            
            allROIs_4{ff, 1} = VS1centroids;
            allROIs_4{ff, 2} = VS2centroids;
            
            C{1} = VS1centroids;
            C{2} = VS2centroids;

            nCells(1) = size(VS1centroids,1);
            nCells(2) = size(VS2centroids,1);

            %smallest numbers of cells
            [nUse, iidx] = min(nCells);

            for nn=1:length(C)

                if nn == iidx   % if the cell group is the one with smallerst cell number
                    Dis = pdist2(C{nn},C{nn}, 'euclidean');
                    Dis = unique(Dis);
                    Dis = Dis(2:end);%remove zero
                    distancesInfos_4(ff,nn) = nanmean(Dis);

                else    %if the cell group is the one with a larger cell number

                    CAll = C{nn};
                    disSelfTry = nan(1, ntry);
                    disOtherTry = nan(1, ntry);

                    for t = 1:ntry
                        CAllTry=CAll(randperm(size(CAll,1)),:);
                        CAllTry=CAllTry(1:nUse,:);%random nUse cells

                        %self distances
                        Dis = pdist2(CAllTry,CAllTry, 'euclidean');
                        Dis = unique(Dis);
                        Dis = Dis(2:end);%remove zero
                        disSelfTry(t) = nanmean(Dis);

                        %distances across groups
                        otherGroup = setdiff([1 2],nn);
                        Dis = pdist2(CAllTry, C{otherGroup}, 'euclidean');
                        disOtherTry(t) = nanmean(Dis,'all');%this one no duplicate values
                    end

                    distancesInfos_4(ff,nn) = nanmean(disSelfTry);
                    distancesInfos_4(ff,3) = nanmean(disOtherTry);
                end
            end
            
            ff = ff + 1;
        end
        end
        cd ..\
    end
    cd ..\
end
cd ..\
toc;
end



% PLOTTING

allDist = distancesInfos_4(:,[1 3 2]);%first column is auditory, second is between, the last is visual.

s=sum(allDist,2);
allDist=allDist(~isnan(s),:); %getting rid of the ones with nan

figure; hold on
for n=1:size(allDist,1)
    plot(allDist(n,:), 'c--', 'LineWidth', 0.2);
end

errorbar(1, nanmean(allDist(:,1)), nanstd(allDist(:,1))/sqrt(size(allDist,1)),...
    'm', 'LineWidth', 3);
errorbar(2, nanmean(allDist(:,2)), nanstd(allDist(:,2))/sqrt(size(allDist,1)),...
    'k', 'LineWidth', 3);
errorbar(3, nanmean(allDist(:,3)), nanstd(allDist(:,3))/sqrt(size(allDist,1)),...
    'g', 'LineWidth', 3);
xlim([0.8, 3.2]); ylabel('Average Distance (\mum)');
set(gca, 'xtick', 1:3, 'xticklabel', {'vs1-vs1', 'vs1-vs2', 'vs2-vs2'});
set(gcf, 'Position', [222 222 345 333]);


[~, p1] = ttest(allDist(:,1), allDist(:,2)); disp(p1);
[~, p2] = ttest(allDist(:,1), allDist(:,3)); disp(p2);
[~, p3] = ttest(allDist(:,2), allDist(:,3)); disp(p3);

% [isSignificant,adjusted_pvals,~]= bonferroni_holm([p1 p2 p3],0.05); % adjusted_pvals are the corrected p values
% disp(adjusted_pvals);
% 
% title(['12:',num2str(adjusted_pvals(1)),' 13:',num2str(adjusted_pvals(2)),' 23:',num2str(adjusted_pvals(3))])



figure;
for ff = 1:58
    subplot(6,10,ff); hold on
    xlim([0 512]); ylim([0 512]);
    plot(allROIs_4{ff,1}(:,1), allROIs_4{ff,1}(:,2), 'k.');
    plot(allROIs_4{ff,2}(:,1), allROIs_4{ff,2}(:,2), 'g.');
end














%% CELL ANATOMY - AO1 vs AO2 - only AO1 in AO1-VS1

distancesInfos_5 = cell(222, 4);      % AO1 vs AO2

if 1
tic;
cd('Z:\labMembers\DN\_PROJECT_\finalData\allAlignments\audio_old_vs_new_afterLearned');
dd1 = dir('ID202*');
ff = 1;

%yg
% nCells={};
% useIdx=[];

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
            env1 = AOn(idx1); env2 = AOo(idx2);
            useIdx(end+1,:)=[nMice nFov idx1 idx2];
            % load common cells
            fl = dir([num2str(env1) '\' num2str(env1) '_' num2str(env2) '\cellRegistered*.mat']);
            load([num2str(env1) '\' num2str(env1) '_' num2str(env2) '\' fl(1).name]);
            cellIdx = cell_registered_struct.cell_to_index_map;
            
            % find unique-cell populations
            cc = find(sum(cellIdx, 2) > 0);
            Uao2 = cellIdx(cc, 1); Uao2 = unique(Uao2); Uao2 = Uao2(2:end);
            Uao1 = cellIdx(cc, 2); Uao1 = unique(Uao1); Uao1 = Uao1(2:end);

            % centers
            AO2centroids = cell_registered_struct.centroid_locations_corrected{1};
            AO2centroids = AO2centroids(Uao2,:);
            AO1centroids = cell_registered_struct.centroid_locations_corrected{2};
            AO1centroids = AO1centroids(Uao1,:);

            sz = max( size(AO2centroids,1), size(AO1centroids,1) );
            % distances
            % AOtoAO
            AO22 = pdist2( AO2centroids, AO2centroids, 'euclidean', 'Smallest', sz);

            % AOtoVS
            AO21 = pdist2( AO2centroids, AO1centroids, 'euclidean', 'Smallest', sz);

            % VStoAO
            AO12 = pdist2( AO1centroids, AO2centroids, 'euclidean', 'Smallest', sz);

            % VStoVS
            AO11 = pdist2( AO1centroids, AO1centroids, 'euclidean', 'Smallest', sz);

            for ss = 1:sz
                if ss <= size(AO22,1)-1
                    distancesInfos_5{ ff, 1} = [distancesInfos_5{ ff, 1}; nanmean( AO22(2:ss+1,:), 1)];
                end

                if ss <= size(AO11,1)-1
                    distancesInfos_5{ ff, 4} = [distancesInfos_5{ ff, 4}; nanmean( AO11(2:ss+1,:), 1)];
                end

                if ss <= size(AO21,1)
                    distancesInfos_5{ ff, 2} = [distancesInfos_5{ ff, 2}; nanmean( AO21(1:ss, :), 1)];
                end

                if ss <= size(AO12,1)
                    distancesInfos_5{ ff, 3} = [distancesInfos_5{ ff, 3}; nanmean( AO12(1:ss, :), 1)];
                end

            end
            
            ff = ff + 1;
        end
        end
        
        cd ..\
       
    end
    cd ..\
end
cd ..\
toc;
end


%% PLOTTING

figure; hold on
allDist = [];
for ff = 1:size(distancesInfos_5,1)
    if isempty(distancesInfos_5{ff, 1})
        break
    end
    if ismember(ff, [4 8 10 78 88 89 90 91 92])
        continue
    end
    mDist = nan(1,3);
    mDist(1) = nanmean(distancesInfos_5{ff, 1}, 'all');
    mDist(2) = nanmean(distancesInfos_5{ff, 2}, 'all');
    mDist(3) = nanmean(distancesInfos_5{ff, 4}, 'all');
    allDist = [allDist; mDist];
    plot(mDist, 'c--', 'LineWidth', 0.2);
end
errorbar(1, nanmean(allDist(:,1)), nanstd(allDist(:,1))/sqrt(size(allDist,1)),...
    'm', 'LineWidth', 3);
errorbar(2, nanmean(allDist(:,2)), nanstd(allDist(:,2))/sqrt(size(allDist,1)),...
    'k', 'LineWidth', 3);
errorbar(3, nanmean(allDist(:,3)), nanstd(allDist(:,3))/sqrt(size(allDist,1)),...
    'g', 'LineWidth', 3);
xlim([0.8, 3.2]); ylabel('Average Distance (\mum)');
set(gca, 'xtick', 1:3, 'xticklabel', {'ao2-ao2', 'ao1-ao2', 'ao1-ao1'});
set(gcf, 'Position', [222 222 345 333]);
title('AO1 vs AO2');


[~, p] = ttest2(allDist(:,1), allDist(:,2)); disp(p);
[~, p] = ttest2(allDist(:,1), allDist(:,3)); disp(p);
[~, p] = ttest2(allDist(:,2), allDist(:,3)); disp(p);

% i=setdiff([1:1:size(allDist,1)],[20 21]);
% i=setdiff([1:1:size(allDist,1)],[32 33]);
% i=setdiff([1:1:size(allDist,1)],[36 37]);
% [~, p] = ttest2(allDist(i,1), allDist(i,2)); disp(p);
% [~, p] = ttest2(allDist(i,1), allDist(i,3)); disp(p);
% [~, p] = ttest2(allDist(i,2), allDist(i,3)); disp(p);



%% combine ALL

figure; hold on
allDist = [];
for ff = 1:size(distancesInfos,1)
    if ff == 3
        continue
    end
    mDist = nan(1,3);
    mDist(1) = nanmean(distancesInfos{ff, 1}, 'all');
    mDist(2) = nanmean(distancesInfos{ff, 2}, 'all');
    mDist(3) = nanmean(distancesInfos{ff, 4}, 'all');
    allDist = [allDist; mDist];
    plot(mDist, 'c--', 'LineWidth', 0.2);
end
for ff = 1:size(distancesInfos_2,1)
    if ff == 43 || ff == 44 || ff == 47
        continue
    end
    mDist = nan(1,3);
    mDist(1) = nanmean(distancesInfos_2{ff, 1}, 'all');
    mDist(2) = nanmean(distancesInfos_2{ff, 2}, 'all');
    mDist(3) = nanmean(distancesInfos_2{ff, 4}, 'all');
    allDist = [allDist; mDist];
    plot(mDist, 'c--', 'LineWidth', 0.2);
end
for ff = 1:size(distancesInfos_3,1)
    if isempty(distancesInfos_3{ff,1})
        continue
    end
    mDist = nan(3,3);
    mDist(1,1) = nanmean(distancesInfos_3{ff, 1}, 'all');
    mDist(2,1) = nanmean(distancesInfos_3{ff, 2}, 'all');
    mDist(3,1) = nanmean(distancesInfos_3{ff, 4}, 'all');
    allDist = [allDist; mDist];
    plot(1 - 0.1 + 0.2 * rand(3,1), mDist(:,1), 'k.');
end
for ff = 1:size(distancesInfos_4,1)
    if isempty(distancesInfos_4{ff,1})
        continue
    end
    mDist = nan(3,3);
    mDist(1,3) = nanmean(distancesInfos_4{ff, 1}, 'all');
    mDist(2,3) = nanmean(distancesInfos_4{ff, 2}, 'all');
    mDist(3,3) = nanmean(distancesInfos_4{ff, 4}, 'all');
    allDist = [allDist; mDist];
    plot(3 - 0.1 + 0.2 * rand(3,1), mDist(:,3), 'k.');
end
errorbar(1, nanmean(allDist(:,1)), nanstd(allDist(:,1))/sqrt(size(allDist,1)),...
    'm', 'LineWidth', 3);
errorbar(2, nanmean(allDist(:,2)), nanstd(allDist(:,2))/sqrt(size(allDist,1)),...
    'k', 'LineWidth', 3);
errorbar(3, nanmean(allDist(:,3)), nanstd(allDist(:,3))/sqrt(size(allDist,1)),...
    'g', 'LineWidth', 3);
xlim([0.8, 3.2]); ylabel('Average Distance (\mum)');
set(gca, 'xtick', 1:3, 'xticklabel', {'ao-ao', 'ao-vs', 'vs-vs'});
set(gcf, 'Position', [222 222 345 333]);
title('AO vs VS');


[~, p] = ttest2(allDist(:,1), allDist(:,2)); disp(p);
[~, p] = ttest2(allDist(:,1), allDist(:,3)); disp(p);
[~, p] = ttest2(allDist(:,2), allDist(:,3)); disp(p);



%%
