%load data: contact lead author

%%
folders=folders_combo;
allFieldsA = cell(1, numel(folders_combo));  
allDfofsA = cell(1, numel(folders_combo));
allAverageDfofsA= cell(1, numel(folders_combo));
rbrConsistenciesA = cell(1, numel(folders_combo));
cueConsistenciesA = cell(3, numel(folders_combo));   % VS/AO/AO-6-bin - toOthers

allFieldsV = cell(1, numel(folders_combo));  
allDfofsV = cell(1, numel(folders_combo));
allAverageDfofsV= cell(1, numel(folders_combo));
rbrConsistenciesV = cell(1, numel(folders_combo));
cueConsistenciesV = cell(3, numel(folders_combo));   % VS/AO/AO-6-bin - toOthers

allFieldsN = cell(1, numel(folders_combo));  
allDfofsN = cell(1, numel(folders_combo));
allAverageDfofsN= cell(1, numel(folders_combo));
rbrConsistenciesN = cell(1, numel(folders_combo));
cueConsistenciesN = cell(3, numel(folders_combo));   % VS/AO/AO-6-bin - toOthers

% rbrConsistencies = cell(1, numel(folders_combo));
% cueConsistencies = cell(3, numel(folders_combo));   % VS/AO/AO-6-bin - toOthers
% nRuns = nan(numel(folders_combo), 2);               % number of success runs / number of runs
% cueScores = cell(1, numel(folders_combo));          % for each FOV - audio-10-bin-full; audio-6-bin-full; audio-10-bin-2cue; audio-6-bin-2cue; visual
% cueLags = cell(1, numel(folders_combo));            % for each FOV - audio-10-bin-full; audio-6-bin-full; audio-10-bin-2cue; audio-6-bin-2cue; visual
% allDfofMs = [];
nBin=140;

zones{1,1} = 9:18;          % VS - 10-bin
zones{1,2} = 81:90;         % VS - 10-bin
zones{1,3} = 95:104;        % VS - 10-bin
zones{1,4} = 127:136;       % VS - 10-bin

zones{2,1} = 29:38;         % AO - 10-bin
zones{2,2} = 43:52;         % AO - 10-bin
zones{2,3} = 63:72;         % AO - 10-bin
zones{2,4} = 111:120;       % AO - 10-bin

zones{3,1} = 31:36;         % AO - 6-bin
zones{3,2} = 45:50;         % AO - 6-bin
zones{3,3} = 65:70;         % AO - 6-bin
zones{3,4} = 113:118;       % AO - 6-bin

%

nRunAVN=[]; %number of runs for auditory reward only, visual only, and no reward at all
p=pwd;


for ff = 1:numel(folders)
  % progress
    disp([num2str(ff) '/' num2str(numel(folders))]);
cd(folders{ff})
cd('visualAuditoryOnly_20230819')
load('isVRNoFirstLast.mat');
load([folders{ff} '\corrIncorr_20230701\dfofMInterpM.mat']);
    dfofM = dfofMInterpM; clear dfofMInterpM;
    load([folders{ff} '\roiIdxUse.mat']);
    
    load([folders{ff} '\RunByRun_sig\dfofMInterpM_sig.mat']);
    dfofM2 = dfofMInterpM_sig; clear dfofMInterpM_sig;

 if exist('PValueClassifier_KY2_6_sigA','dir')
     nRunAVN(ff,1)=length(find(isVRNoFirstLast==0));
    % load all data
    load('PValueClassifier_KY2_6_sigA\allCells.mat');
    fields = allCells.inFieldBins;
    dfof_sig = allCells.dfofaveragesmooth;
    dfof_moments = allCells.dfofUse;
    clear allCells;
    load('dfofInterpM_sigA.mat')
    

    for nn = 1:numel(fields)

        if ~roiIdxUse(nn)
            continue
        end
        if any(isnan(dfofM{nn}), 'all')
            continue
        end

         % FIELDS
        tempBINs = fields{nn};
        tempZeros = zeros(1, nBin);
        if ~isnan(tempBINs)
            tempZeros(tempBINs) = 1;
        end
        allFieldsA{ff} = [allFieldsA{ff}; tempZeros];

                   % DFOF

                allDfofsA{ff} = [allDfofsA{ff}; dfof_sig(:, nn)'];
     allAverageDfofsA{ff} = [allAverageDfofsA{ff} nanmean( dfof_moments(:, nn), 'all')];
        rbr = dfofInterpM_sigA{nn};
        if any(~isnan(rbr), 'all')

            % RBR-CONSISTENCY
            corrInfo = calculateRBR_singleCell( rbr );
            rbrConsistenciesA{ff} = [rbrConsistenciesA{ff}; [corrInfo.meantoMean corrInfo.meantoNext corrInfo.meantoOthers]];

            % CUE-CONSISTENCY
            for cueType = 1:3
                toOthers = nan(1, nBin);
                for cues = 1:4
                    tempRBR = rbr( :, zones{cueType, cues} );
                    corrInfo = calculateRBR_singleCellDUC( tempRBR );
                    if isnan(corrInfo.meantoOthers); corrInfo.meantoOthers = 0; end
                    toOthers(zones{cueType, cues}) = corrInfo.meantoOthers;
                end
                cueConsistenciesA{cueType, ff} = [cueConsistenciesA{cueType, ff}; toOthers]; clear toOthers;
            end
            
        else
            % RBR-CONSISTENCY
            rbrConsistenciesA{ff} = [rbrConsistenciesA{ff}; nan(1,3)];
            
            % CUE-CONSISTENCY
            for cueType = 1:3
                cueConsistenciesA{cueType, ff} = [cueConsistenciesA{cueType, ff}; nan(1, nBin)];
            end
        end
     
    end
 else
     nRunAVN(ff,1)=0;
       allFieldsA{ff} = [];

                   % DFOF

                allDfofsA{ff} = [];
     allAverageDfofsA{ff} = [];
     rbrConsistenciesA{ff}=[];
     for cueType = 1:3
                cueConsistenciesA{cueType, ff} = [];
            end

 end


  if exist('PValueClassifier_KY2_6_sigV','dir')
     nRunAVN(ff,2)=length(find(isVRNoFirstLast==1));
    % load all data
    load('PValueClassifier_KY2_6_sigV\allCells.mat');
    fields = allCells.inFieldBins;
    dfof_sig = allCells.dfofaveragesmooth;
    dfof_moments = allCells.dfofUse;
    clear allCells;
      load('dfofInterpM_sigV.mat')


    for nn = 1:numel(fields)

        if ~roiIdxUse(nn)
            continue
        end
        if any(isnan(dfofM{nn}), 'all')
            continue
        end

         % FIELDS
        tempBINs = fields{nn};
        tempZeros = zeros(1, nBin);
        if ~isnan(tempBINs)
            tempZeros(tempBINs) = 1;
        end
        allFieldsV{ff} = [allFieldsV{ff}; tempZeros];

                   % DFOF

                allDfofsV{ff} = [allDfofsV{ff}; dfof_sig(:, nn)'];
     allAverageDfofsV{ff} = [allAverageDfofsV{ff} nanmean( dfof_moments(:, nn), 'all')];
 rbr = dfofInterpM_sigV{nn};
        if any(~isnan(rbr), 'all')

            % RBR-CONSISTENCY
            corrInfo = calculateRBR_singleCell( rbr );
            rbrConsistenciesV{ff} = [rbrConsistenciesV{ff}; [corrInfo.meantoMean corrInfo.meantoNext corrInfo.meantoOthers]];

            % CUE-CONSISTENCY
            for cueType = 1:3
                toOthers = nan(1, nBin);
                for cues = 1:4
                    tempRBR = rbr( :, zones{cueType, cues} );
                    corrInfo = calculateRBR_singleCellDUC( tempRBR );
                    if isnan(corrInfo.meantoOthers); corrInfo.meantoOthers = 0; end
                    toOthers(zones{cueType, cues}) = corrInfo.meantoOthers;
                end
                cueConsistenciesV{cueType, ff} = [cueConsistenciesV{cueType, ff}; toOthers]; clear toOthers;
            end
            
        else
            % RBR-CONSISTENCY
            rbrConsistenciesV{ff} = [rbrConsistenciesV{ff}; nan(1,3)];
            
            % CUE-CONSISTENCY
            for cueType = 1:3
                cueConsistenciesV{cueType, ff} = [cueConsistenciesV{cueType, ff}; nan(1, nBin)];
            end
        end
     
    end
 else
     nRunAVN(ff,2)=0;
       allFieldsV{ff} = [];

                   % DFOF

                allDfofsV{ff} = [];
     allAverageDfofsV{ff} = [];
          rbrConsistenciesV{ff}=[];
     for cueType = 1:3
                cueConsistenciesV{cueType, ff} = [];
            end
 end

 cd ..\
 cd('noRewardOnly_20230820')
load('isCorrectNoFirstLast.mat')
if exist('PValueClassifier_KY2_6_sigIC','dir')
    nRunAVN(ff,3)=length(find(isCorrectNoFirstLast==2));

  % load all data
    load('PValueClassifier_KY2_6_sigIC\allCells.mat');
    fields = allCells.inFieldBins;
    dfof_sig = allCells.dfofaveragesmooth;
    dfof_moments = allCells.dfofUse;
    clear allCells;
  load('dfofInterpM_sigIC.mat')

    for nn = 1:numel(fields)

        if ~roiIdxUse(nn)
            continue
        end
        if any(isnan(dfofM{nn}), 'all')
            continue
        end

         % FIELDS
        tempBINs = fields{nn};
        tempZeros = zeros(1, nBin);
        if ~isnan(tempBINs)
            tempZeros(tempBINs) = 1;
        end
        allFieldsN{ff} = [allFieldsN{ff}; tempZeros];

                   % DFOF

                allDfofsN{ff} = [allDfofsN{ff}; dfof_sig(:, nn)'];
     allAverageDfofsN{ff} = [allAverageDfofsN{ff} nanmean( dfof_moments(:, nn), 'all')];
 rbr = dfofInterpM_sigIC{nn};
        if any(~isnan(rbr), 'all')

            % RBR-CONSISTENCY
            corrInfo = calculateRBR_singleCell( rbr );
            rbrConsistenciesN{ff} = [rbrConsistenciesN{ff}; [corrInfo.meantoMean corrInfo.meantoNext corrInfo.meantoOthers]];

            % CUE-CONSISTENCY
            for cueType = 1:3
                toOthers = nan(1, nBin);
                for cues = 1:4
                    tempRBR = rbr( :, zones{cueType, cues} );
                    corrInfo = calculateRBR_singleCellDUC( tempRBR );
                    if isnan(corrInfo.meantoOthers); corrInfo.meantoOthers = 0; end
                    toOthers(zones{cueType, cues}) = corrInfo.meantoOthers;
                end
                cueConsistenciesN{cueType, ff} = [cueConsistenciesN{cueType, ff}; toOthers]; clear toOthers;
            end
            
        else
            % RBR-CONSISTENCY
            rbrConsistenciesN{ff} = [rbrConsistenciesN{ff}; nan(1,3)];
            
            % CUE-CONSISTENCY
            for cueType = 1:3
                cueConsistenciesN{cueType, ff} = [cueConsistenciesN{cueType, ff}; nan(1, nBin)];
            end
        end
     
    end
 else
     nRunAVN(ff,3)=0;
       allFieldsN{ff} = [];

                   % DFOF

                allDfofsN{ff} = [];
     allAverageDfofsN{ff} = [];
               rbrConsistenciesN{ff}=[];
     for cueType = 1:3
                cueConsistenciesN{cueType, ff} = [];
            end
 end

end

cd(p)

save('infoAVN.mat','nRunAVN','allFieldsV','allDfofsV','allAverageDfofsV','rbrConsistenciesV','cueConsistenciesV','allFieldsA','allDfofsA','allAverageDfofsA','rbrConsistenciesA','cueConsistenciesA','allFieldsN','allDfofsN','allAverageDfofsN','rbrConsistenciesN','cueConsistenciesN')

%% also load Z:\labMembers\DN\_PROJECT_\finalData\figures\20230705_v10\fig7c.mat


i3=nRunAVN(:,1).*nRunAVN(:,2).*nRunAVN(:,3);
i3=find(i3>0);%all categories have runs
i=sum(nRunAVN,2);
usei=find(i>=3);
usei=intersect(usei,i3);
i=find(nRuns(:,1)>1);
i=intersect(usei,i);

%only use the fov with idential number of cells
numberdfof=[];
for n=1:length(allDfofs);
numberdfof(n,1)=size(allDfofs{n},1);
end
for n=1:length(allDfofsA);
    if isempty(allDfofsA{n});
       numberdfof(n,2)=nan; 
    else
numberdfof(n,2)=size(allDfofsA{n},1);
    end
end
for n=1:length(allDfofsV);
    if isempty(allDfofsV{n});
       numberdfof(n,3)=nan; 
    else
numberdfof(n,3)=size(allDfofsV{n},1);
    end
end
di=[];
di(:,1)=numberdfof(:,1)-numberdfof(:,2);
di(:,2)=numberdfof(:,1)-numberdfof(:,3);
s=sum(di,2);
ii=find(s==0);
ifinal=intersect(ii,i);

save('ifinal.mat','ifinal')

%% two rewards
fields = cat(1, allDfofs{ifinal});
for n=1:size(fields,1);
    fields(n,:)=naninterp(fields(n,:));
end

shF={};
nBin=140;
nShuffle=100;%kept individual shuffles
for n=1:nShuffle;
    s=[];
for m=1:size(fields,1);
    f=fields(m,:);
    s(m,:)=f(randperm(length(f)));
end
shF{n}=s;
end

% 10-bin
aoIdx = [29:38 43:52 63:72 111:120];
vsIdx = [9:18 81:90 95:104 127:136];
% a=zeros(1,140);
% a(aoIdx)=1;
% a(vsIdx)=1;
% figure,plot(a)
d=[];

avr = fields(:, aoIdx); avr = nanmean(avr, 2);
vvr = fields(:, vsIdx); vvr = nanmean(vvr, 2);
d = (avr - vvr) ./ (avr + vvr);
d2=[];

for n=1:nShuffle;
shAvr = shF{n}(:, aoIdx); 
shAvr = nanmean(shAvr, 2);
shVvr = shF{n}(:, vsIdx); 
shVvr = nanmean(shVvr, 2);

d2 (:,n)= (shAvr - shVvr) ./ (shAvr + shVvr);
end

edges = linspace(-1,1, 22 );
[N, ~] = histcounts(d, edges);
N = N / sum(N); N = N*100;
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'k', 'FaceAlpha', 0.5 );

NShuffle=[];
for n=1:nShuffle;
[NS, ~] = histcounts(d2(:,n), edges);
NS = NS / sum(NS); NS = NS*100;
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'r', 'FaceAlpha', 0.5 );
NShuffle(n,:)=NS;
end
%plot above: actually used 1 st and 99 percentile
thresh5=prctile(reshape(d2,[size(d2,1)*size(d2,2) 1]),1);
thresh95=prctile(reshape(d2,[size(d2,1)*size(d2,2) 1]),99);

figure
subplot(141)
for n=1:size(NShuffle,1);
    hold on
    plot(edges(1:end-1) + diff(edges)/2,NShuffle(n,:),'Color',[0.8 0.8 0.8])
end


hold on
semshade(nanmean(NShuffle,1),0.2,'k',edges(1:end-1) + diff(edges)/2,1);
hold on
semshade(N,0.2,'r',edges(1:end-1) + diff(edges)/2,1);
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'k', 'FaceAlpha', 0.5 );

hold on
line([thresh5 thresh5],[0 max(N)]);
hold on
line([thresh95 thresh95],[0 max(N)]);
title('2 rewards only')

save('d.mat','d')
save('thresh5.mat','thresh5')
save('thresh95.mat','thresh95')

%% auditory cue only

fieldsA = cat(1, allDfofsA{ifinal});
for n=1:size(fieldsA,1);
    fieldsA(n,:)=naninterp(fieldsA(n,:));
end

shFA={};
nBin=140;
nShuffle=100;%kept individual shuffles
for n=1:nShuffle;
    s=[];
for m=1:size(fieldsA,1);
    f=fieldsA(m,:);
    s(m,:)=f(randperm(length(f)));
end
shFA{n}=s;
end

% 10-bin
aoIdx = [29:38 43:52 63:72 111:120];
vsIdx = [9:18 81:90 95:104 127:136];
% a=zeros(1,140);
% a(aoIdx)=1;
% a(vsIdx)=1;
% figure,plot(a)

da=[];
avr = fieldsA(:, aoIdx); avr = nanmean(avr, 2);
vvr = fieldsA(:, vsIdx); vvr = nanmean(vvr, 2);
da = (avr - vvr) ./ (avr + vvr);
d2=[];

for n=1:nShuffle;
shAvr = shFA{n}(:, aoIdx); 
shAvr = nanmean(shAvr, 2);
shVvr = shFA{n}(:, vsIdx); 
shVvr = nanmean(shVvr, 2);

d2 (:,n)= (shAvr - shVvr) ./ (shAvr + shVvr);
end

edges = linspace(-1,1, 22 );
[N, ~] = histcounts(da, edges);
N = N / sum(N); N = N*100;
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'k', 'FaceAlpha', 0.5 );

NShuffle=[];
for n=1:nShuffle;
[NS, ~] = histcounts(d2(:,n), edges);
NS = NS / sum(NS); NS = NS*100;
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'r', 'FaceAlpha', 0.5 );
NShuffle(n,:)=NS;
end
%plot above
thresh5A=prctile(reshape(d2,[size(d2,1)*size(d2,2) 1]),1);
thresh95A=prctile(reshape(d2,[size(d2,1)*size(d2,2) 1]),99);


subplot(142)
for n=1:size(NShuffle,1);
    hold on
    plot(edges(1:end-1) + diff(edges)/2,NShuffle(n,:),'Color',[0.8 0.8 0.8])
end


hold on
semshade(nanmean(NShuffle,1),0.2,'k',edges(1:end-1) + diff(edges)/2,1);
hold on
semshade(N,0.2,'r',edges(1:end-1) + diff(edges)/2,1);
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'k', 'FaceAlpha', 0.5 );

hold on
line([thresh5A thresh5A],[0 max(N)]);
hold on
line([thresh95A thresh95A],[0 max(N)]);
title('auditory only')

save('da.mat','da')

%% visual cue only

fieldsV = cat(1, allDfofsV{ifinal});
for n=1:size(fieldsV,1);
    fieldsV(n,:)=naninterp(fieldsV(n,:));
end

shFV={};
nBin=140;
nShuffle=100;%kept individual shuffles
for n=1:nShuffle;
    s=[];
for m=1:size(fieldsV,1);
    f=fieldsV(m,:);
    s(m,:)=f(randperm(length(f)));
end
shFV{n}=s;
end

% 10-bin
aoIdx = [29:38 43:52 63:72 111:120];
vsIdx = [9:18 81:90 95:104 127:136];
% a=zeros(1,140);
% a(aoIdx)=1;
% a(vsIdx)=1;
% figure,plot(a)

dv=[];
avr = fieldsV(:, aoIdx); avr = nanmean(avr, 2);
vvr = fieldsV(:, vsIdx); vvr = nanmean(vvr, 2);
dv = (avr - vvr) ./ (avr + vvr);
d2=[];

for n=1:nShuffle;
shAvr = shFV{n}(:, aoIdx); 
shAvr = nanmean(shAvr, 2);
shVvr = shFV{n}(:, vsIdx); 
shVvr = nanmean(shVvr, 2);

d2 (:,n)= (shAvr - shVvr) ./ (shAvr + shVvr);
end

edges = linspace(-1,1, 22 );
[N, ~] = histcounts(dv, edges);
N = N / sum(N); N = N*100;
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'k', 'FaceAlpha', 0.5 );

NShuffle=[];
for n=1:nShuffle;
[NS, ~] = histcounts(d2(:,n), edges);
NS = NS / sum(NS); NS = NS*100;
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'r', 'FaceAlpha', 0.5 );
NShuffle(n,:)=NS;
end
%plot above
thresh5V=prctile(reshape(d2,[size(d2,1)*size(d2,2) 1]),1);
thresh95V=prctile(reshape(d2,[size(d2,1)*size(d2,2) 1]),99);


subplot(143)
for n=1:size(NShuffle,1);
    hold on
    plot(edges(1:end-1) + diff(edges)/2,NShuffle(n,:),'Color',[0.8 0.8 0.8])
end


hold on
semshade(nanmean(NShuffle,1),0.2,'k',edges(1:end-1) + diff(edges)/2,1);
hold on
semshade(N,0.2,'r',edges(1:end-1) + diff(edges)/2,1);
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'k', 'FaceAlpha', 0.5 );

hold on
line([thresh5V thresh5V],[0 max(N)]);
hold on
line([thresh95V thresh95V],[0 max(N)]);
title('visual only')

save('dv.mat','dv')


%% no reward

fieldsN = cat(1, allDfofsN{ifinal});
for n=1:size(fieldsN,1);
    fieldsN(n,:)=naninterp(fieldsN(n,:));
end

shFN={};
nBin=140;
nShuffle=100;%kept individual shuffles
for n=1:nShuffle;
    s=[];
for m=1:size(fieldsN,1);
    f=fieldsN(m,:);
    s(m,:)=f(randperm(length(f)));
end
shFN{n}=s;
end

% 10-bin
aoIdx = [29:38 43:52 63:72 111:120];
vsIdx = [9:18 81:90 95:104 127:136];
save('aoIdx.mat','aoIdx');
save('vsIdx.mat','vsIdx');

% a=zeros(1,140);
% a(aoIdx)=1;
% a(vsIdx)=1;
% figure,plot(a)

dn=[];
avr = fieldsN(:, aoIdx); avr = nanmean(avr, 2);
vvr = fieldsN(:, vsIdx); vvr = nanmean(vvr, 2);
dn = (avr - vvr) ./ (avr + vvr);
dn(isnan(dn))=0;
d2=[];
% dn(isnan(dn))=0;
for n=1:nShuffle;
shAvr = shFN{n}(:, aoIdx); 
shAvr = nanmean(shAvr, 2);
shVvr = shFN{n}(:, vsIdx); 
shVvr = nanmean(shVvr, 2);

d2 (:,n)= (shAvr - shVvr) ./ (shAvr + shVvr);
end

edges = linspace(-1,1, 22 );
[N, ~] = histcounts(dn, edges);
N = N / sum(N); N = N*100;
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'k', 'FaceAlpha', 0.5 );

NShuffle=[];
for n=1:nShuffle;
[NS, ~] = histcounts(d2(:,n), edges);
NS = NS / sum(NS); NS = NS*100;
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'r', 'FaceAlpha', 0.5 );
NShuffle(n,:)=NS;
end
%plot above
thresh5N=prctile(reshape(d2,[size(d2,1)*size(d2,2) 1]),1);
thresh95N=prctile(reshape(d2,[size(d2,1)*size(d2,2) 1]),99);


subplot(144)
for n=1:size(NShuffle,1);
    hold on
    plot(edges(1:end-1) + diff(edges)/2,NShuffle(n,:),'Color',[0.8 0.8 0.8])
end


hold on
semshade(nanmean(NShuffle,1),0.2,'k',edges(1:end-1) + diff(edges)/2,1);
hold on
semshade(N,0.2,'r',edges(1:end-1) + diff(edges)/2,1);
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'k', 'FaceAlpha', 0.5 );

hold on
line([thresh5N thresh5N],[0 max(N)]);
hold on
line([thresh95N thresh95N],[0 max(N)]);
title('no reward')
save('dn.mat','dn')

saveas(gcf,'allConsitionSensoryScores.fig')
%

%% use the thresh used in previous data
load('Z:\labMembers\DN\_PROJECT_\finalData\figures\20230705_v10\fig_YG_DFOF_completeRandom_99\thresh.mat')
%again here thresh 95 and 5 are actually 99 and 1

aoCellIdx=find(d>thresh95);
vsCellIdx=find(d<thresh5);
aovsCellIdx=find(d>=thresh5&d<=thresh95);

aoCellIdxa=find(da>thresh95);
vsCellIdxa=find(da<thresh5);
aovsCellIdxa=find(da>=thresh5&da<=thresh95);

aoCellIdxv=find(dv>thresh95);
vsCellIdxv=find(dv<thresh5);
aovsCellIdxv=find(dv>=thresh5&dv<=thresh95);

aoCellIdxn=find(dn>thresh95);
vsCellIdxn=find(dn<thresh5);
aovsCellIdxn=find(dn>=thresh5&dn<=thresh95);
save('thresh5.mat','thresh5')
save('thresh95.mat','thresh95')
save('allIdx.mat','aoCellIdx','vsCellIdx','aovsCellIdx','aoCellIdxa','vsCellIdxa','aovsCellIdxa','aoCellIdxv','vsCellIdxv','aovsCellIdxv','aoCellIdxn','vsCellIdxn','aovsCellIdxn')
figure
nTotal = length(d);
subplot(141)
nAO = length(aoCellIdx);
nVS = length(vsCellIdx);
ax = gca();
% pie(ax, [(nTotal-nAO-nVS) nAO nVS], [1 1 1], {'other', 'audio', 'visual'}');
pie(ax, [(nTotal-nAO-nVS) nAO nVS], {'other', 'audio', 'visual'}');

ax.Colormap = [...
    0.8 0.8 0.8;
    1 0 1;
    0 1 0];
set(gcf, 'position', [222 222 333 333]);
title('both reward')

subplot(142)
nAO = length(aoCellIdxa);
nVS = length(vsCellIdxa);
ax = gca();
% pie(ax, [(nTotal-nAO-nVS) nAO nVS], [1 1 1], {'other', 'audio', 'visual'}');
pie(ax, [(nTotal-nAO-nVS) nAO nVS], {'other', 'audio', 'visual'}');

ax.Colormap = [...
    0.8 0.8 0.8;
    1 0 1;
    0 1 0];
set(gcf, 'position', [222 222 333 333]);
title('auditory')

subplot(143)
nAO = length(aoCellIdxv);
nVS = length(vsCellIdxv);
ax = gca();
% pie(ax, [(nTotal-nAO-nVS) nAO nVS], [1 1 1], {'other', 'audio', 'visual'}');
pie(ax, [(nTotal-nAO-nVS) nAO nVS], {'other', 'audio', 'visual'}');

ax.Colormap = [...
    0.8 0.8 0.8;
    1 0 1;
    0 1 0];
set(gcf, 'position', [222 222 333 333]);
title('visual')

subplot(144)
nAO = length(aoCellIdxn);
nVS = length(vsCellIdxn);
ax = gca();
% pie(ax, [(nTotal-nAO-nVS) nAO nVS], [1 1 1], {'other', 'audio', 'visual'}');
pie(ax, [(nTotal-nAO-nVS) nAO nVS], {'other', 'audio', 'visual'}');

ax.Colormap = [...
    0.8 0.8 0.8;
    1 0 1;
    0 1 0];
set(gcf, 'position', [222 222 333 333]);
title('no reward')

saveas(gcf,'cellFractions4Conditions.fig')

%%
%cell overlap WITH two reward condition
figure
subplot(131);
a1=aoCellIdx;
b1=aoCellIdxa;
c1=aoCellIdxv;
d1=aoCellIdxn;

iab=length(intersect(a1,b1));
iac=length(intersect(a1,c1));
iad=length(intersect(a1,d1));

bar([1 2 3],[iab/(length(a1)+length(b1)-iab) iac/(length(a1)+length(c1)-iac) iac/(length(a1)+length(d1)-iad)])
names={'2r-a';'2r-v';'2r-n'};
set(gca,'xtick',[1:3],'xticklabel',names)

title('ao cells')

subplot(132);
a1=vsCellIdx;
b1=vsCellIdxa;
c1=vsCellIdxv;
d1=vsCellIdxn;

iab=length(intersect(a1,b1));
iac=length(intersect(a1,c1));
iad=length(intersect(a1,d1));

bar([1 2 3],[iab/(length(a1)+length(b1)-iab) iac/(length(a1)+length(c1)-iac) iac/(length(a1)+length(d1)-iad)])
names={'2r-a';'2r-v';'2r-n'};
set(gca,'xtick',[1:3],'xticklabel',names)

title('vs cells')

subplot(133);
a1=aovsCellIdx;
b1=aovsCellIdxa;
c1=aovsCellIdxv;
d1=aovsCellIdxn;

iab=length(intersect(a1,b1));
iac=length(intersect(a1,c1));
iad=length(intersect(a1,d1));

bar([1 2 3],[iab/(length(a1)+length(b1)-iab) iac/(length(a1)+length(c1)-iac) iac/(length(a1)+length(d1)-iad)])
names={'2r-a';'2r-v';'2r-n'};
set(gca,'xtick',[1:3],'xticklabel',names)

title('aovs cells')
saveas(gcf,'overlapWith2Rewards.fig')

%% random chance
randomOverlap={};
for n=1:3;
    randomOverlap{n}=[];
end
%three categories: ao, vs, aovs cells overlaps

nTotal=length(d);

for n=1:100;
    dr=randperm(nTotal);
    aoCellIdxs=dr(1:length(aoCellIdx));
    vsCellIdxs=dr(length(aoCellIdx)+1:length(vsCellIdx)+length(aoCellIdx));
    aovsCellIdxs=dr(length(vsCellIdx)+length(aoCellIdx)+1:end);

        da=randperm(nTotal);
    aoCellIdxas=da(1:length(aoCellIdxa));
    vsCellIdxas=da(length(aoCellIdxa)+1:length(vsCellIdxa)+length(aoCellIdxa));
    aovsCellIdxas=da(length(vsCellIdxa)+length(aoCellIdxa)+1:end);

      dv=randperm(nTotal);
    aoCellIdxvs=dv(1:length(aoCellIdxv));
    vsCellIdxvs=dv(length(aoCellIdxv)+1:length(vsCellIdxv)+length(aoCellIdxv));
    aovsCellIdxvs=dv(length(vsCellIdxv)+length(aoCellIdxv)+1:end);

     dn=randperm(nTotal);
    aoCellIdxns=dn(1:length(aoCellIdxn));
    vsCellIdxns=dn(length(aoCellIdxn)+1:length(vsCellIdxn)+length(aoCellIdxn));
    aovsCellIdxns=dn(length(vsCellIdxn)+length(aoCellIdxn)+1:end);

    a1=aoCellIdxs;
b1=aoCellIdxas;
c1=aoCellIdxvs;
d1=aoCellIdxns;

iab=length(intersect(a1,b1));
iac=length(intersect(a1,c1));
iad=length(intersect(a1,d1));


    randomOverlap{1}(n,:)=[iab/(length(a1)+length(b1)-iab) iac/(length(a1)+length(c1)-iac) iac/(length(a1)+length(d1)-iad)];

 a1=vsCellIdxs;
b1=vsCellIdxas;
c1=vsCellIdxvs;
d1=vsCellIdxns;

iab=length(intersect(a1,b1));
iac=length(intersect(a1,c1));
iad=length(intersect(a1,d1));


    randomOverlap{2}(n,:)=[iab/(length(a1)+length(b1)-iab) iac/(length(a1)+length(c1)-iac) iac/(length(a1)+length(d1)-iad)];
    
 a1=aovsCellIdxs;
b1=aovsCellIdxas;
c1=aovsCellIdxvs;
d1=aovsCellIdxns;

iab=length(intersect(a1,b1));
iac=length(intersect(a1,c1));
iad=length(intersect(a1,d1));


    randomOverlap{3}(n,:)=[iab/(length(a1)+length(b1)-iab) iac/(length(a1)+length(c1)-iac) iac/(length(a1)+length(d1)-iad)];
end

figure,
for n=1:3;
    subplot(1,3,n);
    bar([1 2 3],nanmean(randomOverlap{n},1))
    hold on
    errorbar([1 2 3],nanmean(randomOverlap{n},1),nansem(randomOverlap{n},1),'.');
    names={'2r-a';'2r-v';'2r-n'};
set(gca,'xtick',[1:3],'xticklabel',names)
if n==1;
    title('ao cells')
elseif n==2;
    title('vs cells')
else
    title('aovs cells')
end
end

saveas(gcf,'shuffledOverlaps.fig')
%%
%calcualte fraction of cells per fov
f=allDfofs(ifinal);
for n=1:length(f);
    for m=1:size(f{n},1)
    f{n}(m,:)=naninterp(f{n}(m,:));
    end
end

dfindv={};
for n=1:length(f)
    fThis=f{n};
avr =  fThis(:, aoIdx); avr = nanmean(avr, 2);
vvr =  fThis(:, vsIdx); vvr = nanmean(vvr, 2);
dd = (avr - vvr) ./ (avr + vvr);
dd(isnan(dd))=0;
dfindv{n}=dd;
end

fa=allDfofsA(ifinal);
for n=1:length(fa);
    for m=1:size(fa{n},1)
    fa{n}(m,:)=naninterp(fa{n}(m,:));
    end
end

daindv={};
for n=1:length(fa)
    fThis=fa{n};
avr =  fThis(:, aoIdx); avr = nanmean(avr, 2);
vvr =  fThis(:, vsIdx); vvr = nanmean(vvr, 2);
dd = (avr - vvr) ./ (avr + vvr);
dd(isnan(dd))=0;
daindv{n}=dd;
end

fv=allDfofsV(ifinal);
for n=1:length(fv);
    for m=1:size(fv{n},1)
    fv{n}(m,:)=naninterp(fv{n}(m,:));
    end
end

dvindv={};
for n=1:length(fv)
    fThis=fv{n};
avr =  fThis(:, aoIdx); avr = nanmean(avr, 2);
vvr =  fThis(:, vsIdx); vvr = nanmean(vvr, 2);
dd = (avr - vvr) ./ (avr + vvr);
dd(isnan(dd))=0;
dvindv{n}=dd;
end

fn=allDfofsN(ifinal);
for n=1:length(fn);
    for m=1:size(fn{n},1)
    fn{n}(m,:)=naninterp(fn{n}(m,:));
    end
end

dnindv={};
for n=1:length(fn)
    fThis=fn{n};
avr =  fThis(:, aoIdx); avr = nanmean(avr, 2);
vvr =  fThis(:, vsIdx); vvr = nanmean(vvr, 2);
dd = (avr - vvr) ./ (avr + vvr);
dd(isnan(dd))=0;
dnindv{n}=dd;
end

%cells classified

allCellTypes=[]; %percentage. first cell is 2 rewards, second is auditory reward, third is visua, last is no reward
aThis=[];
dThis=dfindv;
for n=1:length(dfindv);
    ao=length(find(dThis{n}>thresh95));
    vs=length(find(dThis{n}<thresh5));
    aovs=length(find(dThis{n}>=thresh5&dThis{n}<=thresh95));
    aThis(n,1)=ao/length(dThis{n});
    aThis(n,2)=vs/length(dThis{n});
    aThis(n,3)= aovs/length(dThis{n});
end
allCellTypes{1}=aThis;

aThis=[];
dThis=daindv;
for n=1:length(dfindv);
    ao=length(find(dThis{n}>thresh95));
    vs=length(find(dThis{n}<thresh5));
    aovs=length(find(dThis{n}>=thresh5&dThis{n}<=thresh95));
    aThis(n,1)=ao/length(dThis{n});
    aThis(n,2)=vs/length(dThis{n});
    aThis(n,3)= aovs/length(dThis{n});
end
allCellTypes{2}=aThis;

aThis=[];
dThis=dvindv;
for n=1:length(dfindv);
    ao=length(find(dThis{n}>thresh95));
    vs=length(find(dThis{n}<thresh5));
    aovs=length(find(dThis{n}>=thresh5&dThis{n}<=thresh95));
    aThis(n,1)=ao/length(dThis{n});
    aThis(n,2)=vs/length(dThis{n});
    aThis(n,3)= aovs/length(dThis{n});
end
allCellTypes{3}=aThis;

aThis=[];
dThis=dnindv;
for n=1:length(dfindv);
    ao=length(find(dThis{n}>thresh95));
    vs=length(find(dThis{n}<thresh5));
    aovs=length(find(dThis{n}>=thresh5&dThis{n}<=thresh95));
    aThis(n,1)=ao/length(dThis{n});
    aThis(n,2)=vs/length(dThis{n});
    aThis(n,3)= aovs/length(dThis{n});
end
allCellTypes{4}=aThis;

allCellTypes2={}; %organize by cell type: A, V and AV, each column is each behavior category: two reward, auditory only, visual only, no reward
for n=1:3
    allCellTypes2{n}=[];
    for m=1:4;
        allCellTypes2{n}(:,m)=allCellTypes{m}(:,n);
    end
end
save('allCellTypes.mat','allCellTypes')
save('allCellTypes2.mat','allCellTypes2')

aoPFOV=[];
for n=1:3;
    [~,aoPFOV(n)]=ttest(allCellTypes2{1}(:,1),allCellTypes2{1}(:,n+1));
end

vsPFOV=[];
for n=1:3;
    [~,vsPFOV(n)]=ttest(allCellTypes2{2}(:,1),allCellTypes2{2}(:,n+1));
end

aovsPFOV=[];
for n=1:3;
    [~,aovsPFOV(n)]=ttest(allCellTypes2{3}(:,1),allCellTypes2{3}(:,n+1));
end

figure,
for n=1:length(allCellTypes2)
subplot(1,4,n)
plotThis=allCellTypes2{n};
bar([1:1:4],nanmean(plotThis,1))
hold on
errorbar([1:1:4],nanmean(plotThis,1),nansem(plotThis,1),'.');
names={'2r';'a';'v';'n'};

set(gca,'xtick',[1:4],'xticklabel',names)
if n==1;
    title('auditory cells')
elseif n==2;
    title('visual cells')
else
    title('visual auditory cells')
end
end

%ration of visual and auditory
ratioAV=[];
for n=1:length(allCellTypes);
    ratioAV(:,n)=allCellTypes{n}(:,1)./allCellTypes{n}(:,2);
end

ratioP=[];
for n=1:3;
    [~,ratioP(n)]=ttest(ratioAV(:,1),ratioAV(:,n+1));
end

subplot(144);
bar([1:1:4],nanmean(ratioAV,1))
hold on
errorbar([1:1:4],nanmean(ratioAV,1),nansem(ratioAV,1),'.');
names={'2r';'a';'v';'n'};
set(gca,'xtick',[1:4],'xticklabel',names)

title('ratio aoCells/vsCells');
saveas(gcf,'aoVSCellsPerFOV.fig')
%% the field distributions of all cells

all4Fields={};
figure
subplot(311)
fieldsThis=allFields(ifinal);
fieldsAll1=[];
for n=1:length(fieldsThis)
    fieldsAll1(n,:)=sum(fieldsThis{n},1)/size(fieldsThis{n},1);
end
semshade(fieldsAll1,0.2,'b',[1:1:140],0);
all4Fields{1}=fieldsAll1;
hold on
fieldsThis=allFieldsA(ifinal);
fieldsAll2=[];
for n=1:length(fieldsThis)
    fieldsAll2(n,:)=sum(fieldsThis{n},1)/size(fieldsThis{n},1);
end
semshade(fieldsAll2,0.2,'m',[1:1:140],0);
ylim([0 0.6])
all4Fields{2}=fieldsAll2;

diffp=[];
for n=1:size(fieldsAll1,2)
    [~,diffp(n)]=ttest(fieldsAll1(:,n),fieldsAll2(:,n),'tail','right');
end

for n=1:length(diffp);
    if diffp(n)<0.05;
        hold on
        plot(n,0.65,'b*')
    end
end

ylim([0 0.7])
title('compare two reward and auditory only')


subplot(312)
fieldsThis=allFields(ifinal);
fieldsAll1=[];
for n=1:length(fieldsThis)
    fieldsAll1(n,:)=sum(fieldsThis{n},1)/size(fieldsThis{n},1);
end
semshade(fieldsAll1,0.2,'b',[1:1:140],0);

hold on
fieldsThis=allFieldsV(ifinal);
fieldsAll2=[];
for n=1:length(fieldsThis)
    fieldsAll2(n,:)=sum(fieldsThis{n},1)/size(fieldsThis{n},1);
end
semshade(fieldsAll2,0.2,'m',[1:1:140],0);
ylim([0 0.6])
all4Fields{3}=fieldsAll2;
diffp=[];
for n=1:size(fieldsAll1,2)
    [~,diffp(n)]=ttest(fieldsAll1(:,n),fieldsAll2(:,n),'tail','right');
end

for n=1:length(diffp);
    if diffp(n)<0.05;
        hold on
        plot(n,0.65,'b*')
    end
end

ylim([0 0.7])
title('compare two reward and visual only')

subplot(313)
fieldsThis=allFields(ifinal);
fieldsAll1=[];
for n=1:length(fieldsThis)
    fieldsAll1(n,:)=sum(fieldsThis{n},1)/size(fieldsThis{n},1);
end
semshade(fieldsAll1,0.2,'b',[1:1:140],0);

hold on
fieldsThis=allFieldsN(ifinal);
fieldsAll2=[];
for n=1:length(fieldsThis)
    fieldsAll2(n,:)=sum(fieldsThis{n},1)/size(fieldsThis{n},1);
end
semshade(fieldsAll2,0.2,'m',[1:1:140],0);
ylim([0 0.6])
all4Fields{4}=fieldsAll2;
diffp=[];
for n=1:size(fieldsAll1,2)
    [~,diffp(n)]=ttest(fieldsAll1(:,n),fieldsAll2(:,n),'tail','right');
end

for n=1:length(diffp);
    if diffp(n)<0.05;
        hold on
        plot(n,0.65,'b*')
    end
end

ylim([0 0.7])
title('compare two reward and no reward')
saveas(gcf,'compareFields.fig')
save('all4Fields.mat','all4Fields')
%% calculate the fields in auditory cue zones and visual cue zones
aoIdx = [29:38 43:52 63:72 111:120];
vsIdx = [9:18 81:90 95:104 127:136];
aovsidx{1}=aoIdx;
aovsidx{2}=vsIdx;

fieldsAOVS={};
for n=1:2;
    fieldsAOVS{n}=[];
    for m=1:length(all4Fields);
        fThis=all4Fields{m};
        fieldsAOVS{n}(:,m)=mean(fThis(:,aovsidx{n}),2);
    end
end

figure,
for n=1:2
subplot(1,2,n);
bar([1:1:4],mean(fieldsAOVS{n},1))
hold on
errorbar([1:1:4],mean(fieldsAOVS{n},1),nansem(fieldsAOVS{n},1),'.')
names={'2r';'a';'v';'n'};
set(gca,'xtick',[1:4],'xticklabel',names)
if n==1;
    title('fields at auditory cues')
else
    title('fields at visual cues')
end
end
saveas(gcf,'compareFields_AO_VS.fig')


%% dfof

all4Dfof={};
fields = cat(1, allDfofs{ifinal});
for n=1:size(fields,1);
    fields(n,:)=naninterp(fields(n,:));
end
all4Dfof{1}=fields;

fields = cat(1, allDfofsA{ifinal});
for n=1:size(fields,1);
    fields(n,:)=naninterp(fields(n,:));
end
all4Dfof{2}=fields;

fields = cat(1, allDfofsV{ifinal});
for n=1:size(fields,1);
    fields(n,:)=naninterp(fields(n,:));
end
all4Dfof{3}=fields;

fields = cat(1, allDfofsN{ifinal});
for n=1:size(fields,1);
    fields(n,:)=naninterp(fields(n,:));
end
all4Dfof{4}=fields;

aoIdx = [29:38 43:52 63:72 111:120];
vsIdx = [9:18 81:90 95:104 127:136];
aovsidx{1}=aoIdx;
aovsidx{2}=vsIdx;

dfofAOVS={};
for n=1:2;
    dfofAOVS{n}=[];
    for m=1:length(all4Dfof);
        fThis=all4Dfof{m};
        dfofAOVS{n}(:,m)=mean(fThis(:,aovsidx{n}),2);
    end
end

figure,
for n=1:2
subplot(1,2,n);
bar([1:1:4],mean(dfofAOVS{n},1))
hold on
errorbar([1:1:4],mean(dfofAOVS{n},1),nansem(dfofAOVS{n},1),'.')
names={'2r';'a';'v';'n'};
set(gca,'xtick',[1:4],'xticklabel',names)
if n==1;
    title('DFOF at auditory cues')
else
    title('DFOF at visual cues')
end
end

saveas(gcf,'compareDFOF.fig')

%% AUDITORY CELLS, CUE CONSISTENCY AT VISUAL AND AUDITORY CUES

%AO cells at ao cues
cueCAOAOAll={};
%ao cells in two reward
cueCAO=cueConsistencies(2,:);
cueCAO = cat(1, cueCAO{ifinal});
A=nanmean(cueCAO(aoCellIdx,:),2);
A(A==0)=nan;
cueCAOAOAll{1}=A(~isnan(A));

%ao cells in auditory only 
cueCAO=cueConsistenciesA(2,:);
cueCAO = cat(1, cueCAO{ifinal});
A=nanmean(cueCAO(aoCellIdxa,:),2);
A(A==0)=nan;
cueCAOAOAll{2}=A(~isnan(A));

%ao cells in visual only 
cueCAO=cueConsistenciesV(2,:);
cueCAO = cat(1, cueCAO{ifinal});
A=nanmean(cueCAO(aoCellIdxv,:),2);
A(A==0)=nan;
cueCAOAOAll{3}=A(~isnan(A));

%ao cells in no reward
cueCAO=cueConsistenciesN(2,:);
cueCAO = cat(1, cueCAO{ifinal});
A=nanmean(cueCAO(aoCellIdxn,:),2);
A(A==0)=nan;
cueCAOAOAll{4}=A(~isnan(A));

cueCAOAllMean=[];
cueCAOAllSEM=[];
for n=1:length(cueCAOAOAll);
    cueCAOAllMean(n)=mean(cueCAOAOAll{n});
    cueCAOAllSEM(n)=nansem(cueCAOAOAll{n},1);
end

figure
subplot(221)
bar([1:1:4],cueCAOAllMean)
hold on
errorbar([1:1:4],cueCAOAllMean,cueCAOAllSEM,'.')
title('AUDITORY CELLS at auditory cue')

names={'2r';'a';'v';'n'};
set(gca,'xtick',[1:4],'xticklabel',names)

%AO cells at VS cues
cueCAOVSAll={};
%ao cells in two reward
cueCVS=cueConsistencies(1,:);
cueCVS = cat(1, cueCVS{ifinal});
A=nanmean(cueCVS(aoCellIdx,:),2);
A(A==0)=nan;
cueCAOVSAll{1}=A(~isnan(A));

%ao cells in auditory only 
cueCVS=cueConsistenciesA(1,:);
cueCVS = cat(1, cueCVS{ifinal});
A=nanmean(cueCVS(aoCellIdxa,:),2);
A(A==0)=nan;
cueCAOVSAll{2}=A(~isnan(A));

%ao cells in visual only 
cueCVS=cueConsistenciesV(1,:);
cueCVS = cat(1, cueCVS{ifinal});
A=nanmean(cueCVS(aoCellIdxv,:),2);
A(A==0)=nan;
cueCAOVSAll{3}=A(~isnan(A));

%ao cells in no reward
cueCVS=cueConsistenciesN(1,:);
cueCVS = cat(1, cueCVS{ifinal});
A=nanmean(cueCVS(aoCellIdxn,:),2);
A(A==0)=nan;
cueCAOVSAll{4}=A(~isnan(A));

cueCVSAllMean=[];
cueCVSAllSEM=[];
for n=1:length(cueCAOVSAll);
    cueCVSAllMean(n)=mean(cueCAOVSAll{n});
    cueCVSAllSEM(n)=nansem(cueCAOVSAll{n},1);
end

subplot(222)

bar([1:1:4],cueCVSAllMean)
hold on
errorbar([1:1:4],cueCVSAllMean,cueCVSAllSEM,'.')
title('AUDITORY CELLS at visual cues')
names={'2r';'a';'v';'n'};
set(gca,'xtick',[1:4],'xticklabel',names)

%AO cells at ao cues
cueCVSAOAll={};
%ao cells in two reward
cueCAO=cueConsistencies(2,:);
cueCAO = cat(1, cueCAO{ifinal});
A=nanmean(cueCAO(vsCellIdx,:),2);
A(A==0)=nan;
cueCVSAOAll{1}=A(~isnan(A));

%ao cells in auditory only 
cueCAO=cueConsistenciesA(2,:);
cueCAO = cat(1, cueCAO{ifinal});
A=nanmean(cueCAO(vsCellIdxa,:),2);
A(A==0)=nan;
cueCVSAOAll{2}=A(~isnan(A));

%ao cells in visual only 
cueCAO=cueConsistenciesV(2,:);
cueCAO = cat(1, cueCAO{ifinal});
A=nanmean(cueCAO(vsCellIdxv,:),2);
A(A==0)=nan;
cueCVSAOAll{3}=A(~isnan(A));

%ao cells in no reward
cueCAO=cueConsistenciesN(2,:);
cueCAO = cat(1, cueCAO{ifinal});
A=nanmean(cueCAO(vsCellIdxn,:),2);
A(A==0)=nan;
cueCVSAOAll{4}=A(~isnan(A));

cueCAOAllMean=[];
cueCAOAllSEM=[];
for n=1:length(cueCVSAOAll);
    cueCAOAllMean(n)=mean(cueCVSAOAll{n});
    cueCAOAllSEM(n)=nansem(cueCVSAOAll{n},1);
end


subplot(223)
bar([1:1:4],cueCAOAllMean)
hold on
errorbar([1:1:4],cueCAOAllMean,cueCAOAllSEM,'.')
title('VISUAL CELLS at auditory cue')

names={'2r';'a';'v';'n'};
set(gca,'xtick',[1:4],'xticklabel',names)

%AO cells at VS cues
cueCVSVSAll={};
%ao cells in two reward
cueCVS=cueConsistencies(1,:);
cueCVS = cat(1, cueCVS{ifinal});
A=nanmean(cueCVS(vsCellIdx,:),2);
A(A==0)=nan;
cueCVSVSAll{1}=A(~isnan(A));

%ao cells in auditory only 
cueCVS=cueConsistenciesA(1,:);
cueCVS = cat(1, cueCVS{ifinal});
A=nanmean(cueCVS(vsCellIdxa,:),2);
A(A==0)=nan;
cueCVSVSAll{2}=A(~isnan(A));

%ao cells in visual only 
cueCVS=cueConsistenciesV(1,:);
cueCVS = cat(1, cueCVS{ifinal});
A=nanmean(cueCVS(vsCellIdxv,:),2);
A(A==0)=nan;
cueCVSVSAll{3}=A(~isnan(A));

%ao cells in no reward
cueCVS=cueConsistenciesN(1,:);
cueCVS = cat(1, cueCVS{ifinal});
A=nanmean(cueCVS(vsCellIdxn,:),2);
A(A==0)=nan;
cueCVSVSAll{4}=A(~isnan(A));

cueCVSAllMean=[];
cueCVSAllSEM=[];
for n=1:length(cueCVSVSAll);
    cueCVSAllMean(n)=mean(cueCVSVSAll{n});
    cueCVSAllSEM(n)=nansem(cueCVSVSAll{n},1);
end

subplot(224)

bar([1:1:4],cueCVSAllMean)
hold on
errorbar([1:1:4],cueCVSAllMean,cueCVSAllSEM,'.')
title('VISUAL CELLS at visual cues')
names={'2r';'a';'v';'n'};
set(gca,'xtick',[1:4],'xticklabel',names)

saveas(gcf,'cueConsistencyAllCells')
%% ao and vs cells: all consistency
aovsRunC={};
%the first one is auditory cells, the second one is visual cells
aovsRunC{1}={};
%ao cells in 2 reward
aoRunC=cat(1, rbrConsistencies{ifinal});
aovsRunC{1}{1}=aoRunC(aoCellIdx,3);
%ao cells in auditory only reward
aoRunC=cat(1, rbrConsistenciesA{ifinal});
aovsRunC{1}{2}=aoRunC(aoCellIdxa,3);
%ao cells in visual only reward
aoRunC=cat(1, rbrConsistenciesV{ifinal});
aovsRunC{1}{3}=aoRunC(aoCellIdxv,3);
%ao cells in no reward
aoRunC=cat(1, rbrConsistenciesN{ifinal});
aovsRunC{1}{4}=aoRunC(aoCellIdxn,3);

%the first one is auditory cells, the second one is visual cells
aovsRunC{2}={};
%ao cells in 2 reward
vsRunC=cat(1, rbrConsistencies{ifinal});
aovsRunC{2}{1}=vsRunC(vsCellIdx,3);
%ao cells in auditory only reward
vsRunC=cat(1, rbrConsistenciesA{ifinal});
aovsRunC{2}{2}=vsRunC(vsCellIdxa,3);
%ao cells in visual only reward
vsRunC=cat(1, rbrConsistenciesV{ifinal});
aovsRunC{2}{3}=vsRunC(vsCellIdxv,3);
%ao cells in no reward
vsRunC=cat(1, rbrConsistenciesN{ifinal});
aovsRunC{2}{4}=vsRunC(vsCellIdxn,3);

aovsRunCMean=[];
aovsRunCSEM=[];
for n=1:length(aovsRunC)
    for m=1:length(aovsRunC{n})
    aovsRunCMean(n,m)=nanmean(aovsRunC{n}{m});
    aovsRunCSEM(n,m)=nansem(aovsRunC{n}{m},1);
    end
end

figure
for n=1:2;
    subplot(1,2,n);
    bar([1:1:4],aovsRunCMean(n,:));
    hold on
    errorbar([1:1:4],aovsRunCMean(n,:),aovsRunCSEM(n,:),'.');
    if n==1;
    title('AUDITORY CELLS RBR consistency')
    else
      title('VISUAL CELLS RBR consistency')  
    end
names={'2r';'a';'v';'n'};
set(gca,'xtick',[1:4],'xticklabel',names)
end

aoP=[];
for n=1:3;
[~,aoP(n)]=ttest2(aovsRunC{1}{1},aovsRunC{1}{n+1});
end

vsP=[];
for n=1:3;
[~,vsP(n)]=ttest2(aovsRunC{2}{1},aovsRunC{2}{n+1});
end

saveas(gcf,'runConsistency.fig')

%% plot for paper figure

load("thresh5.mat")
load("thresh95.mat")

figure
subplot(141);
load('d.mat');

edges = linspace(-1,1, 22 );
[N, ~] = histcounts(d, edges);
N = N / sum(N); N = N*100;
semshade(N,0.2,'r',edges(1:end-1) + diff(edges)/2,1);
hold on
line([thresh5 thresh5],[0 max(N)]);
hold on
line([thresh95 thresh95],[0 max(N)]);
title('2 rewards only')
ylim tight

subplot(142);
load('da.mat');

edges = linspace(-1,1, 22 );
[N, ~] = histcounts(da, edges);
N = N / sum(N); N = N*100;
semshade(N,0.2,'r',edges(1:end-1) + diff(edges)/2,1);
hold on
line([thresh5 thresh5],[0 max(N)]);
hold on
line([thresh95 thresh95],[0 max(N)]);
title('auditory only')
ylim tight

subplot(143);
load('dv.mat');

edges = linspace(-1,1, 22 );
[N, ~] = histcounts(dv, edges);
N = N / sum(N); N = N*100;
semshade(N,0.2,'r',edges(1:end-1) + diff(edges)/2,1);
hold on
line([thresh5 thresh5],[0 max(N)]);
hold on
line([thresh95 thresh95],[0 max(N)]);
title('visual only')
ylim tight

subplot(144);
load('dn.mat');

edges = linspace(-1,1, 22 );
[N, ~] = histcounts(dn, edges);
N = N / sum(N); N = N*100;
semshade(N,0.2,'r',edges(1:end-1) + diff(edges)/2,1);
hold on
line([thresh5 thresh5],[0 max(N)]);
hold on
line([thresh95 thresh95],[0 max(N)]);
title('no reward')
ylim tight

saveas(gcf,'sensoryScores.fig')
%% better plot the shuffles and real overlaps

%only look at auditory and visual cells: 6 numbers:

overlapReal=[];

a1=aoCellIdx;
b1=aoCellIdxa;
c1=aoCellIdxv;
d1=aoCellIdxn;

iab=length(intersect(a1,b1));
iac=length(intersect(a1,c1));
iad=length(intersect(a1,d1));

overlapReal(1:3)=[iab/(length(a1)+length(b1)-iab) iac/(length(a1)+length(c1)-iac) iac/(length(a1)+length(d1)-iad)];

a1=vsCellIdx;
b1=vsCellIdxa;
c1=vsCellIdxv;
d1=vsCellIdxn;

iab=length(intersect(a1,b1));
iac=length(intersect(a1,c1));
iad=length(intersect(a1,d1));

overlapReal(4:6)=[iab/(length(a1)+length(b1)-iab) iac/(length(a1)+length(c1)-iac) iac/(length(a1)+length(d1)-iad)];

a1=aovsCellIdx;
b1=aovsCellIdxa;
c1=aovsCellIdxv;
d1=aovsCellIdxn;

iab=length(intersect(a1,b1));
iac=length(intersect(a1,c1));
iad=length(intersect(a1,d1));

overlapReal(7:9)=[iab/(length(a1)+length(b1)-iab) iac/(length(a1)+length(c1)-iac) iac/(length(a1)+length(d1)-iad)];


%shuffle
randomOverlap={};
for n=1:3;
    randomOverlap{n}=[];
end
%three categories: ao, vs, aovs cells overlaps

nTotal=length(d);

for n=1:100;
    dr=randperm(nTotal);
    aoCellIdxs=dr(1:length(aoCellIdx));
    vsCellIdxs=dr(length(aoCellIdx)+1:length(vsCellIdx)+length(aoCellIdx));
    aovsCellIdxs=dr(length(vsCellIdx)+length(aoCellIdx)+1:end);

        da=randperm(nTotal);
    aoCellIdxas=da(1:length(aoCellIdxa));
    vsCellIdxas=da(length(aoCellIdxa)+1:length(vsCellIdxa)+length(aoCellIdxa));
    aovsCellIdxas=da(length(vsCellIdxa)+length(aoCellIdxa)+1:end);

      dv=randperm(nTotal);
    aoCellIdxvs=dv(1:length(aoCellIdxv));
    vsCellIdxvs=dv(length(aoCellIdxv)+1:length(vsCellIdxv)+length(aoCellIdxv));
    aovsCellIdxvs=dv(length(vsCellIdxv)+length(aoCellIdxv)+1:end);

     dn=randperm(nTotal);
    aoCellIdxns=dn(1:length(aoCellIdxn));
    vsCellIdxns=dn(length(aoCellIdxn)+1:length(vsCellIdxn)+length(aoCellIdxn));
    aovsCellIdxns=dn(length(vsCellIdxn)+length(aoCellIdxn)+1:end);

    a1=aoCellIdxs;
b1=aoCellIdxas;
c1=aoCellIdxvs;
d1=aoCellIdxns;

iab=length(intersect(a1,b1));
iac=length(intersect(a1,c1));
iad=length(intersect(a1,d1));


    randomOverlap{1}(n,:)=[iab/(length(a1)+length(b1)-iab) iac/(length(a1)+length(c1)-iac) iac/(length(a1)+length(d1)-iad)];

 a1=vsCellIdxs;
b1=vsCellIdxas;
c1=vsCellIdxvs;
d1=vsCellIdxns;

iab=length(intersect(a1,b1));
iac=length(intersect(a1,c1));
iad=length(intersect(a1,d1));


    randomOverlap{2}(n,:)=[iab/(length(a1)+length(b1)-iab) iac/(length(a1)+length(c1)-iac) iac/(length(a1)+length(d1)-iad)];
    
 a1=aovsCellIdxs;
b1=aovsCellIdxas;
c1=aovsCellIdxvs;
d1=aovsCellIdxns;

iab=length(intersect(a1,b1));
iac=length(intersect(a1,c1));
iad=length(intersect(a1,d1));


    randomOverlap{3}(n,:)=[iab/(length(a1)+length(b1)-iab) iac/(length(a1)+length(c1)-iac) iac/(length(a1)+length(d1)-iad)];
end

%replot
overlapRandom=[];
overlapRandom(:,1:3)=randomOverlap{1};
overlapRandom(:,4:6)=randomOverlap{2};
overlapRandom(:,7:9)=randomOverlap{3};

figure,
range=0.2;
for n=1:length(overlapReal);
    plot(n-0.5,overlapReal(n),'r.','MarkerSize',10);
    hold on
    x=n-0.5+(rand(size(overlapRandom,1),1)-0.5)*range;
    plot(x,overlapRandom(:,n),'.','Color',[0.8 0.8 0.8],'MarkerSize',3);
    hold on
    errorbar(n-0.5,mean(overlapRandom(:,n)),nansem(overlapRandom(:,n),1),'k')
end

saveas(gcf,'overlap_forPaper.fig')

