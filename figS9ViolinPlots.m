%load data: contact lead author

%%
figure
for n=1:length(allCellTypes2)

    subplot(1,3,n)
a=allCellTypes2{n};

c={};
c{1}='both';%common
c{2}='ao';%unique
c{3}='vs';%unique
c{4}='none';%unique

if n==1;
    color=[1 0 1];
elseif n==2;
    color=[0 1 0];
else
    color=[0.5 0.5 0.5];
end

color=repmat(color,4,1);
grouporder={'both','ao','vs','none'};
%data
All=[];
for n=1:size(a,2);
    All=[All;a(:,n)];
end

%category
AllC={};
for n=1:size(a,2),
    AllC=[AllC;repmat(c(n),size(a,1),1)];
end

violinplot(All, AllC,'ShowData',false,'ViolinColor',color,'ShowMean',true,'ShowWhiskers',false,'ShowBox',false,'GroupOrder',grouporder);
xlim tight;
end
saveas(gcf,'aoVSCellsPerFOV_violin.fig')

%% ao and vs cells: all consistency

load('Z:\labMembers\DN\_PROJECT_\finalData\figures\20230705_v10\fig7.mat')
load('ifinal.mat');
load("infoAVN.mat");
load('allIdx.mat');

%%

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

figure
for n=1:length(aovsRunC);
subplot(1,2,n);
%data
All=[];
for m=1:length(aovsRunC{n});
    All=[All;aovsRunC{n}{m}];
end

c={};
c{1}='both';%common
c{2}='ao';%unique
c{3}='vs';%unique
c{4}='none';%unique

%category
AllC={};
for m=1:length(c),
    AllC=[AllC;repmat(c(m),length(aovsRunC{n}{m}),1)];
end

if n==1;
    color=[1 0 1];
else 
    color=[0 1 0];

end

color=repmat(color,4,1);
grouporder={'both','ao','vs','none'};
violinplot(All, AllC,'ShowData',false,'ViolinColor',color,'ShowMean',true,'ShowWhiskers',false,'ShowBox',false,'GroupOrder',grouporder);
xlim tight;
end

saveas(gcf,'runConsistency_violin.fig')

