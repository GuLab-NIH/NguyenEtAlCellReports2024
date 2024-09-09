
% load data
% contact Lead Author for data


% SPATIAL SELECTIVITY ANALYSIS in common vs unique populations
%   1 - number of fields
%   2 - width of fields
%   3 - infield-nonInField



%% width and number of fields

% % for avr1-vvr1
% useData = allFields_1; mice = 1:31; % 4-19 / 1-31;
% trackLen = 2.1; env = 1; ttl = '1'; 

% % for avr2-vvr2
% useData = allFields_2; mice = 1:60; % 5-32; 1-60
% trackLen = 3.3; env = 2; ttl = '2';


figure;

% avr - Unique
xx = cat(2, useData{1,2,mice});
nn1 = zeros(size(xx,2),1);
wd1 = [];
for ii = 1:size(xx,2)
    try
        yy = contiguous(xx(:,ii),1);
    catch
        continue
    end
    yy = yy{2};
    nn1(ii) = size(yy,1) / trackLen;
    wd1 = [wd1; [yy(:,2)-yy(:,1)+1]*5];
end
subplot(1,2,1); hold on
violin2(nn1, 'x', 2, 'facecolor', 'm', 'facealpha', 0.3);
subplot(1,2,2); hold on
violin2(wd1, 'x', 2, 'facecolor', 'm', 'facealpha', 0.3);


% avr - Common
xx = cat(2, useData{1,1,mice});
nn2 = zeros(size(xx,2),1);
wd2 = [];
for ii = 1:size(xx,2)
    try
        yy = contiguous(xx(:,ii),1);
    catch
        continue
    end
    yy = yy{2};
    nn2(ii) = size(yy,1) / trackLen;
    wd2 = [wd2; [yy(:,2)-yy(:,1)+1]*5];
end
subplot(1,2,1); hold on
violin2(nn2, 'x', 1, 'facecolor', 'k', 'facealpha', 0.3);
subplot(1,2,2); hold on
violin2(wd2, 'x', 1, 'facecolor', 'k', 'facealpha', 0.3);


% visual - Common
xx = cat(2, useData{2,1,mice});
nn3 = zeros(size(xx,2),1);
wd3 = [];
for ii = 1:size(xx,2)
    try
        yy = contiguous(xx(:,ii),1);
    catch
        continue
    end
    yy = yy{2};
    nn3(ii) = size(yy,1) / trackLen;
    wd3 = [wd3; [yy(:,2)-yy(:,1)+1]*5];
end
subplot(1,2,1); hold on
violin2(nn3, 'x', 3, 'facecolor', 'k', 'facealpha', 0.3);
subplot(1,2,2); hold on
violin2(wd3, 'x', 3, 'facecolor', 'k', 'facealpha', 0.3);


% visual - Unique
xx = cat(2, useData{2,2,mice});
nn4 = zeros(size(xx,2),1);
wd4 = [];
for ii = 1:size(xx,2)
    try
        yy = contiguous(xx(:,ii),1);
    catch
        continue
    end
    yy = yy{2};
    nn4(ii) = size(yy,1) / trackLen;
    wd4 = [wd4; [yy(:,2)-yy(:,1)+1]*5];
end
subplot(1,2,1); hold on
violin2(nn4, 'x', 4, 'facecolor', 'g', 'facealpha', 0.3);
subplot(1,2,2); hold on
violin2(wd4, 'x', 4, 'facecolor', 'g', 'facealpha', 0.3);

subplot(1,2,1);
xlim([.5 4.5]);
set(gca, 'xtick', 1:4, 'xticklabel', {'avr-c', 'avr-u', 'vvr-c', 'vvr-u'});
title('number of fields');
subplot(1,2,2);
xlim([.5 4.5]);
set(gca, 'xtick', 1:4, 'xticklabel', {'avr-c', 'avr-u', 'vvr-c', 'vvr-u'});
title('average width of fields');

set(gcf, 'position', [222 222 666 333]);



[~,p] = ttest2(nn1, nn2); disp(p);
[~,p] = ttest2(wd1, wd2); disp(p);

[~,p] = ttest2(nn3, nn4); disp(p);
[~,p] = ttest2(wd3, wd4); disp(p);


[~,p1] = ttest2(nn1, nn4);
[~,p2] = ttest2(nn2, nn3);

figure(15); subplot(1,2,env); hold on
violin2(nn1, 'x', 2, 'facecolor', 'm', 'facealpha', 0.3);
violin2(nn2, 'x', 1, 'facecolor', 'k', 'facealpha', 0.3);
violin2(nn3, 'x', 3, 'facecolor', 'k', 'facealpha', 0.3);
violin2(nn4, 'x', 4, 'facecolor', 'g', 'facealpha', 0.3);
xlim([.5 4.5]);
set(gca, 'xtick', 1:4, 'xticklabel', {'avr-c', 'avr-u', 'vvr-c', 'vvr-u'});
ylabel('# field / meter');
title([ttl ' - ' num2str(p1) '; ' num2str(p2)]);

set(gcf, 'position', [333 333 555 234]);






%% PLOT spatial selectivity


figure; hold on
subplot(1,2,1); hold on
for ii = 1:4
    switch ii
        case 1
            clr = 'k';
        case 2
            clr = 'm';
        case 3
            clr = 'k';
        case 4
            clr = 'g';
    end
    xx = cat(1, inFoutF_1{:,ii});
    dd = (xx(:,1) - xx(:,2)) ./ (xx(:,1) + xx(:,2));   
%     kk{ii} = dd(xx(:,1)>0);
    kk{ii} = dd;
    disp(length(find(kk{ii}<0)) / length(kk{ii}));
    violin2(kk{ii}, 'x', ii, 'facecolor', clr, 'facealpha', 0.3);
end
xlim([.5 4.5]);
set(gca, 'xtick', 1:4, 'xticklabel', {'avr-c', 'avr-u', 'vvr-c', 'vvr-u'});
title('avr1 vs vvr1');

[~,p] = ttest2(kk{1}, kk{2}); disp(p);
[~,p] = ttest2(kk{3}, kk{4}); disp(p);



subplot(1,2,2); hold on
for ii = 1:4
    switch ii
        case 1
            clr = 'k';
        case 2
            clr = 'm';
        case 3
            clr = 'k';
        case 4
            clr = 'g';
    end
    xx = cat(1, inFoutF_2{:,ii});
    dd = (xx(:,1) - xx(:,2)) ./ (xx(:,1) + xx(:,2));
%     kk{ii} = dd(xx(:,1)>0);
    kk{ii} = dd;
    disp(length(find(kk{ii}<0)) / length(kk{ii}));
    violin2(kk{ii}, 'x', ii, 'facecolor', clr, 'facealpha', 0.3);
end
xlim([.5 4.5]);
set(gca, 'xtick', 1:4, 'xticklabel', {'avr-c', 'avr-u', 'vvr-c', 'vvr-u'});
title('avr2 vs vvr2');

set(gcf, 'position', [123 123 666 333]);

[~,p] = ttest2(kk{1}, kk{2}); disp(p);
[~,p] = ttest2(kk{3}, kk{4}); disp(p);




%% test

% useData = inFoutF_1;
useData = inFoutF_2;
for ii = 1:4
for jj = 1:4
    if ii >= jj
        continue
    end
    xx = cat(1, inFoutF_2{:,ii});
    dd1 = (xx(:,1) - xx(:,2)) ./ (xx(:,1) + xx(:,2));
    
    xx = cat(1, inFoutF_2{:,jj});
    dd2 = (xx(:,1) - xx(:,2)) ./ (xx(:,1) + xx(:,2));
    
    [~,p] = kstest2(dd1, dd2);
    disp([num2str(ii) ' vs ' num2str(jj) ': ' num2str(p)]);
end
end






%% GRID CELL ANALYSIS
% avr1 vs vvr1


if 1

load('Y:\labMembers\DN\_PROJECT_\finalData\allFolders.mat');
folders = folders_audio_vs_visual;
clearvars -except folders;


overlapPerc_1 = cell(numel(folders), 1);
inFnoninF_1 = cell(numel(folders), 8);

for ff = 1:numel(folders)
    cd(folders{ff});
    clear fields_ref dfofM_ref SF_ref cell_registered_struct commonCells;
   
    
    try
        % load REFERENCE info
        load([folders{ff} '\audioVR\pcaica\corrIncorr_20230701\PValueClassifier_KY2_6_sigC\allCells.mat']);
        field_ref = allCells.inFieldBins;
        df_ref = allCells.dfofaveragesmooth; clear allCells;
        load([folders{ff} '\audioVR\pcaica\corrIncorr_20230701\PValueClassifier_KY2_6_sigC\grid.mat']);
        grid_ref = grid.indices; clear grid;
        
        
        % load TARGET info
        load([folders{ff} '\visualVR\pcaica\corrIncorr_20230701\PValueClassifier_KY2_6_sigC\allCells.mat']);
        field_target = allCells.inFieldBins;
        df_target = allCells.dfofaveragesmooth; clear allCells;
        load([folders{ff} '\visualVR\pcaica\corrIncorr_20230701\PValueClassifier_KY2_6_sigC\grid.mat']);
        grid_target = grid.indices; clear grid;
        load([folders{ff} '\visualVR\pcaica\roiIdxUse.mat']);
        roiIdxUse_target = find(roiIdxUse); clear roiIdxUse;
%         load([folders{ff} '\visualVR\pcaica\corrIncorr_20230701\cueCells.mat']);
%         isCueCell(isnan(isCueCell)) = 0;
%         cue_target = find(~isCueCell); clear isCueCell cueLags cueScores;
        load([folders{ff} '\visualVR\pcaica\corrIncorr_20230701\cueCells_20240619.mat']); cue_target = find(~cueCells.isCueCell); clear cueCells ;
        
        
        % load commonCells
        fl = dir('1\1_2\cellRegistered*.mat');
        load(['1\1_2\' fl(1).name]); clear fl;
        commonCells = cell_registered_struct.cell_to_index_map;
        commonCells( ~ismember(commonCells(:,2), roiIdxUse_target), 2) = 0;
        commonCells( ~ismember(commonCells(:,2), cue_target), 2) = 0;
        
    catch er
        disp(folders{ff});
        disp(er);
        continue
    end
    
    common = prod(commonCells,2) > 0;
    
    cRefCell = unique(commonCells(common,1));
    uRefCell = unique(commonCells(~common,1)); uRefCell = uRefCell(2:end);
    
    cTargetCell = unique(commonCells(common,2));
    uTargetCell = unique(commonCells(~common,2)); uTargetCell = uTargetCell(2:end);
    
    
    % overlap
    nCommon = sum(common);
    nUniqueRef = numel(uRefCell);
    nUniqueTarget = numel(uTargetCell);
    
    nRefGridCommon = numel(intersect( grid_ref, cRefCell ));
    nRefGridUnique = numel(intersect( grid_ref, uRefCell ));
    
    nTargetGridCommon = numel(intersect( grid_target, cTargetCell ));
    nTargetGridUnique = numel(intersect( grid_target, uTargetCell ));
    
    overlapPerc_1{ff} = [overlapPerc_1{ff}; [nRefGridCommon/nCommon ...
        nRefGridUnique/nUniqueRef ...
        nTargetGridCommon/nCommon ...
        nTargetGridUnique/nUniqueTarget ]];
    
    
    % spatial selectivity
    for type = 1:4
        switch type
            case 1
                indices = cRefCell;
                gridIndices = grid_ref;
                fields = field_ref;
                dfofs = df_ref;
            case 2
                indices = uRefCell;
                gridIndices = grid_ref;
                fields = field_ref;
                dfofs = df_ref;
            case 3
                indices = cTargetCell;
                gridIndices = grid_target;
                fields = field_target;
                dfofs = df_target;
            case 4
                indices = uTargetCell;
                gridIndices = grid_target;
                fields = field_target;
                dfofs = df_target;
        end
        for nn = 1:numel(indices)
            
            % get cell indice
            cInd = indices(nn);
            
            % check grid or not
            if ismember(cInd, gridIndices)
                count = type;
            else
                count = type + 4;
            end
            
            % get dfof infield-nonInField
            if ~isnan(fields{cInd})
                xx = nanmean(dfofs(fields{cInd}, cInd));
                yy = nanmean(dfofs(setdiff(1:42, fields{cInd}), cInd));
            else
                xx = 0;
                yy = nanmean(dfofs(:, cInd));
            end
            inFnoninF_1{ff,count} = [inFnoninF_1{ff,count}; [xx yy]];
        end
    end
    
end

end




%% GRID CELL ANALYSIS
% avr2 vs vvr2


if 1

cd('Z:\labMembers\DN\_PROJECT_\finalData\allAlignments\audio_NEW_vs_visual_NEW');
count2 = 1;

overlapPerc_2 = cell(60, 1);
inFnoninF_2 = cell(60, 8);

dd1 = dir('ID202*');
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
            clear fields_ref dfofM_ref SF_ref cell_registered_struct commonCells;
            try
                % load REFERENCE info
                load([folders{env1} '\corrIncorr_20230701\PValueClassifier_KY2_6_sigC\allCells.mat']);
                field_ref = allCells.inFieldBins;
                df_ref = allCells.dfofaveragesmooth; clear allCells;
                load([folders{env1} '\corrIncorr_20230701\PValueClassifier_KY2_6_sigC\grid.mat']);
                grid_ref = grid.indices; clear grid;
                

                % load TARGET info
                load([folders{env2} '\corrIncorr_20230701\PValueClassifier_KY2_6_sigC\allCells.mat']);
                field_target  = allCells.inFieldBins;
                df_target = allCells.dfofaveragesmooth; clear allCells;
                load([folders{env2} '\corrIncorr_20230701\PValueClassifier_KY2_6_sigC\grid.mat']);
                grid_target = grid.indices; clear grid;
                load([folders{env2} '\roiIdxUse.mat']);
                roiIdxUse_target = find(roiIdxUse); clear roiIdxUse;
%                 load([folders{env2} '\corrIncorr_20230701\cueCells.mat']);
%                 cue_target = find(~isCueCell); clear isCueCell cueLags cueScores;
                load([folders{env2} '\corrIncorr_20230701\cueCells_20240619.mat']); cue_ref = find(~cueCells.isCueCell); clear cueCells ;

                % load common cells
                fl = dir([num2str(env1) '\' num2str(env1) '_' num2str(env2) '\cellRegistered*.mat']);
                load([num2str(env1) '\' num2str(env1) '_' num2str(env2) '\' fl(1).name]);
                commonCells = cell_registered_struct.cell_to_index_map;
                commonCells( ~ismember(commonCells(:,2), roiIdxUse_target), 2) = 0;
                commonCells( ~ismember(commonCells(:,2), cue_target), 2) = 0;
                
            catch er
                disp(er);
                continue
            end
            disp(pwd());
            
            common = prod(commonCells,2) > 0;

            cRefCell = unique(commonCells(common,1));
            uRefCell = unique(commonCells(~common,1)); uRefCell = uRefCell(2:end);

            cTargetCell = unique(commonCells(common,2));
            uTargetCell = unique(commonCells(~common,2)); uTargetCell = uTargetCell(2:end);


            % overlap
            nCommon = sum(common);
            nUniqueRef = numel(uRefCell);
            nUniqueTarget = numel(uTargetCell);

            nRefGridCommon = numel(intersect( grid_ref, cRefCell ));
            nRefGridUnique = numel(intersect( grid_ref, uRefCell ));

            nTargetGridCommon = numel(intersect( grid_target, cTargetCell ));
            nTargetGridUnique = numel(intersect( grid_target, uTargetCell ));

            overlapPerc_2{count2} = [overlapPerc_2{count2}; [nRefGridCommon/nCommon ...
                nRefGridUnique/nUniqueRef ...
                nTargetGridCommon/nCommon ...
                nTargetGridUnique/nUniqueTarget ]];


            % spatial selectivity
            for type = 1:4
                switch type
                    case 1
                        indices = cRefCell;
                        gridIndices = grid_ref;
                        fields = field_ref;
                        dfofs = df_ref;
                    case 2
                        indices = uRefCell;
                        gridIndices = grid_ref;
                        fields = field_ref;
                        dfofs = df_ref;
                    case 3
                        indices = cTargetCell;
                        gridIndices = grid_target;
                        fields = field_target;
                        dfofs = df_target;
                    case 4
                        indices = uTargetCell;
                        gridIndices = grid_target;
                        fields = field_target;
                        dfofs = df_target;
                end
                for nn = 1:numel(indices)

                    % get cell indice
                    cInd = indices(nn);

                    % check grid or not
                    if ismember(cInd, gridIndices)
                        count = type;
                    else
                        count = type + 4;
                    end

                    % get dfof infield-nonInField
                    if ~isnan(fields{cInd})
                        xx = nanmean(dfofs(fields{cInd}, cInd));
                        yy = nanmean(dfofs(setdiff(1:66, fields{cInd}), cInd));
                    else
                        xx = 0;
                        yy = nanmean(dfofs(:, cInd));
                    end
                    inFnoninF_2{count2,count} = [inFnoninF_2{count2,count}; [xx yy]];
                end
            end
            
            count2 = count2 + 1;
        end
        end
        cd ..\
    end
    cd ..\
end
cd ..\

end



%% plot for grid cell analysis
% overlap

% load('extra_21e.mat');

figure; hold on
for type = 1:2
    switch type
        case 1
            xx = cat(1, overlapPerc_1{:});
            ttl = 'avr1 vs. vvr1';
        case 2
            xx = cat(1, overlapPerc_2{:});
            ttl = 'avr2 vs. vvr2';
    end
    
    subplot(1,2,type); hold on
    for kk = 1:4
        switch kk
            case {1,3}
%                 yy = xx(:,2);
                clr = 'k';
            case 2
%                 yy = nanmean(xx(:,[1 3]),2);
                clr = 'm';
            case 4
%                 yy = xx(:,4);
                clr = 'g';
        end
        yy = xx(:,kk);
        violin2(yy, 'x', kk, 'facecolor', clr, 'facealpha', 0.3);
    end
    
    
    [~,p] = ttest( xx(:,1), xx(:,2) ); disp(p);
    [~,p] = ttest( xx(:,3), xx(:,4) ); disp(p);
    
    xlim([.5 4.5]); title(ttl);
%     set(gca, 'xtick', 1:3, 'xticklabel', {'%grid in AVR-u', '%grid in Common', '% grid in VVR-u'}, 'xticklabelrotation', 45);
    set(gca, 'xtick', 1:4, 'xticklabel', {'%grid in AVR-c', '%grid in AVR-u', '%grid in VVR-c', '% grid in VVR-u'}, 'xticklabelrotation', 45);
end
set(gcf, 'position', [123 123 777 333]);


%%

xx = cat(1, overlapPerc_1{:});
yy = cat(1, overlapPerc_2{:});

zz = [xx; yy];
for kk = 1:4
    switch kk
        case {1,3}
%                 yy = xx(:,2);
            clr = 'k';
        case 2
%                 yy = nanmean(xx(:,[1 3]),2);
            clr = 'm';
        case 4
%                 yy = xx(:,4);
            clr = 'g';
    end
    yy = zz(:,kk);
    violin2(yy, 'x', kk, 'facecolor', clr, 'facealpha', 0.3);
end
[~,p] = ttest( zz(:,1), zz(:,2) ); disp(p);
[~,p] = ttest( zz(:,3), zz(:,4) ); disp(p);



%% plot for grid cell analysis
% spatial selectivity

figure; hold on
for VR = 1:2
    switch VR
        case 1
            useData = inFnoninF_1;
            ttl = 'avr1 vs. vvr1';
        case 2
            useData = inFnoninF_2;
            ttl = 'avr2 vs. vvr2';
    end
    subplot(1,2,VR); hold on
    for ii = 1:4
        xx = cat(1, useData{:,ii});
        yy = cat(1, useData{:,ii+4});
        switch ii
            case {1, 3}
                clr = 'k';
            case 2
                clr = 'm';
            case 4
                clr = 'g';
        end
        dd1 = (xx(:,1) - xx(:,2)) ./ (xx(:,1) + xx(:,2));
        dd1 = dd1(xx(:,1)>0);
        violin2(dd1, 'x', 2*ii-1, 'facecolor', clr, 'facealpha', 0.3);
        
        dd2 = (yy(:,1) - yy(:,2)) ./ (yy(:,1) + yy(:,2));
        dd2 = dd2(yy(:,1)>0);
        violin2(dd2, 'x', 2*ii, 'facecolor', clr, 'facealpha', 0.3);
    end
    xlim([.5 8.5]);
    set(gca, 'xtick', 1.5:2:9, 'xticklabel', {'AVR-C', 'AVR-U', 'VVR-C', 'VVR-U'}, 'xticklabelrotation', 45 );
    title(ttl);
end
set(gcf, 'Position', [200 200 999 333]);







%%
