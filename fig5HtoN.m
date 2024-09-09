
%% dfof range

%need to load data: contact lead author for data
%%

%dfof
% randChance=40/700;%0.571, 10%
% rate=nRuns(:,1)./nRuns(:,2)
%remove the one that the sccess rate is 
i=find(nRuns(:,1)>1);
fields = cat(1, allDfofs{i});
for n=1:size(fields,1);
    fields(n,:)=naninterp(fields(n,:));
end

% iUse=find(sum(fields,2)>0); %removing cells that dfof are zeros
% fields=fields(iUse,:);
dfofUse=fields;
save('dfofUse.mat','dfofUse');
% save('iUse.mat','iUse');

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

save('shF.mat','shF');
%%
% 10-bin
aoIdx = [29:38 43:52 63:72 111:120];
vsIdx = [9:18 81:90 95:104 127:136];
% a=zeros(1,140);
% a(aoIdx)=1;
% a(vsIdx)=1;
% figure,plot(a)


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
%plot above
thresh5=prctile(reshape(d2,[size(d2,1)*size(d2,2) 1]),1);
thresh95=prctile(reshape(d2,[size(d2,1)*size(d2,2) 1]),99);

save('thresh.mat','thresh95','thresh5');

figure,
subplot(131)
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'k', 'FaceAlpha', 0.5 );
semshade(N,0.2,'r',edges(1:end-1) + diff(edges)/2,1);
hold on
semshade(NShuffle,0.2,'k',edges(1:end-1) + diff(edges),1);
hold on
line([thresh5 thresh5],[0 max(N)]);
hold on
line([thresh95 thresh95],[0 max(N)]);

ylabel('% of data');
legend('real', 'shuffled');
set(gcf, 'Position', [222 222 234*3 222]);
title('histogram distribution')

subplot(132)

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
subplot(133)
w=0.05;
[p,x]=ksdensity(d,'width',w);
hold on
x2=[];
p2=[];
for n=1:nShuffle;
    [p2(n,:),x2(n,:)]=ksdensity(d2(:,n),'width',w);
end

semshade(p,0.2,'r',x,1);
hold on
semshade(p2,0.2,'k',x,1);

hold on
line([thresh5 thresh5],[0 1]);
hold on
line([thresh95 thresh95],[0 1]);


ylabel('probability');
legend('real', 'shuffled');
set(gcf, 'Position', [222 222 234*3 222]);
title('ksdensity distribution')

saveas(gcf,'Dfof_distributionOfRealShuffle.fig') %this figure actually has sem but since n shuffle =100, the sem looks very small


%%
% 
% thresh5=prctile(reshape(d2,[size(d2,1)*size(d2,2) 1]),5);
% thresh95=prctile(reshape(d2,[size(d2,1)*size(d2,2) 1]),95);

%percentage of data, not absolute data
totald=length(find(~isnan(d)));
real5=length(find(d<thresh5))/totald;
real95=length(find(d>thresh95))/totald;
aoCellIdx=find(d>thresh95);
vsCellIdx=find(d<thresh5);

save('aoCellIdx.mat','aoCellIdx');
save('vsCellIdx.mat','vsCellIdx');
otherIdx=setdiff([1:1:size(fields,1)],[vsCellIdx;aoCellIdx]);
save('otherIdx.mat','otherIdx');
%percentage of each cagetories
shuffle5=[];
shuffle95=[];
shuffle595=[];
shuffleout595=[];
for n=1:nShuffle;
    totald2=length(find(~isnan(d2(:,n))));
    shuffle5(n)=length(find(d2(:,n)<thresh5))/totald2;
    shuffle95(n)=length(find(d2(:,n)>thresh95))/totald2;
      shuffle595(n)=length(find(d2(:,n)>thresh5 & d2(:,n)<thresh95))/totald2;
    shuffleout595(n)=length(find(d2(:,n)>thresh95 | d2(:,n)<thresh5))/totald2;
end

real595=length(find(d>thresh5 & d<thresh95))/totald;
   realout595=length(find(d>thresh95 | d<thresh5))/totald;


%determine 95th percentile and 5th percentile in shuffle, and calcualte how
%many cells fall outside this zone: >95th and <5th percentile.
%and see whether real numbers are above this range compared to shuffle. If
%so, this means that the real data indeed has cells for visual and auditory
figure

subplot(221)
[x,p]=ksdensity(shuffle5);
plot(p,x)

hold on
line([real5 real5],[0 max(x)])
p5=length(find(shuffle5>real5))/nShuffle;
[~,p5]=ttest(shuffle5,real5);
title(['5th p=',num2str(p5)])
xlim([0 0.33])

subplot(222)
[x,p]=ksdensity(shuffle95);
plot(p,x)

hold on
line([real95 real95],[0 max(x)])
p95=length(find(shuffle95>real95))/nShuffle;
[~,p95]=ttest(shuffle95,real95);
title(['95th p=',num2str(p95)])
xlim([0 0.12])

subplot(223)
[x,p]=ksdensity(shuffle595);
plot(p,x)

hold on
line([real595 real595],[0 max(x)])
p5=length(find(shuffle595>real595))/nShuffle;
[~,p5]=ttest(shuffle595,real595);

title(['5-95th p=',num2str(p5)])
xlim([0.55 1])

subplot(224)
[x,p]=ksdensity(shuffleout595);
plot(p,x)

hold on
line([realout595 realout595],[0 max(x)])
p95=length(find(shuffleout595>realout595))/nShuffle;
[~,p95]=ttest(shuffleout595,realout595);

title(['out 5 95th p=',num2str(p95)])
xlim([0 0.45])

saveas(gcf,'Dfof_distributionOfRealShuffle_compareToShuffle.fig')
%% plotting fields


i=find(nRuns(:,1)>1);
fieldsReal = cat(1, allFields{i});
% for n=1:size(fields,1);
%     fields(n,:)=naninterp(fields(n,:));
% end


% fieldsReal=fieldsReal(iUse,:);
fieldDistriUse=fieldsReal;
save('fieldDistriUse.mat','fieldDistriUse');

fieldsReal=fields.*fieldsReal;


a=fieldsReal(aoCellIdx,:);
maxa=[];
for n=1:size(a,1);
    [~,imax]=max(a(n,:));
    maxa(n)=imax;
    a(n,:)=a(n,:)-min(a(n,:));
    a(n,:)=a(n,:)/max(a(n,:));
end

[~,ia]=sort(maxa);
a=a(ia,:);

v=fieldsReal(vsCellIdx,:);
maxv=[];
for n=1:size(v,1);
    [~,imax]=max(v(n,:));
    maxv(n)=imax;
    v(n,:)=v(n,:)-min(v(n,:));
    v(n,:)=v(n,:)/max(v(n,:));
end

[~,iv]=sort(maxv);
v=v(iv,:);

figure,
subplot(121)
imagesc(a);
title('ao')
subplot(122)
imagesc(v);
title('vs')

saveas(gcf,'aovsCellFields.fig')

%% try different combinations of cue selection: visua and audio combined

%generating idx

i=[1:1:8]; %the first four is audio, the second four is visual
aoIdx = [29:38 43:52 63:72 111:120];
vsIdx = [9:18 81:90 95:104 127:136];
zones={};
zones{1}=[29:38];
zones{2}=[43:52];
zones{3}=[63:72];
zones{4}=[111:120];
zones{5}=[9:18];
zones{6}=[81:90];
zones{7}=[95:104];
zones{8}=[127:136];

% fields = cat(1, allDfofs{:});
% for n=1:size(fields,1);
%     fields(n,:)=naninterp(fields(n,:));
% end

NC=1000;
r=[];
for n=1:NC;
    r(n,:)=randperm(8);
    r(n,1:4)=sort(r(n,1:4));
    r(n,5:8)=sort(r(n,5:8));
end

ru=unique(r,'rows');
%removing the rows that are like this:[1 2 3 4 5 6 7 8]
ru=ru(1:end,:); %the first one is the indices we use. The last one isi the complete opposite one
% 
ru=ru(2:35,:); %the first one is the indices we use. The last one isi the complete opposite one

%half of it is symetry: first 35 and the last 35 are reversed.
% ru=ru(2:35,:);%remove the mirrored
% 
%only pick at least two auditory cues (2 or 3)
% a=[1 2 3 4];
% i=[];
% for n=1:size(ru,1);
%     ii=length(intersect(a,ru(n,[1:4])));
%     if ii>=2;
%         i(n)=1;
%     else
%         i(n)=0;
%     end
% end
% 
% ru=ru(find(i),:); 

p5r=[];
p95r=[];
real5r=[];
real95r=[];
real595r=[];
realout595r=[];

dr=[];
for n=1:size(ru,1);
aoIdxr=cell2mat(zones(ru(n,1:4)));
vsIdxr=cell2mat(zones(ru(n,5:8)));
avr = fields(:, aoIdxr); avr = nanmean(avr, 2);
vvr = fields(:, vsIdxr); vvr = nanmean(vvr, 2);
dThis = (avr - vvr) ./ (avr + vvr);% the ratio in this combination
dr(:,n)=dThis; 
end
% first compare to shuffle

for n=1:size(ru,1)
% p5r(n)=length(find(shuffle5>real5r(n)))/nShuffle; %compared to shuffle distritbution: how different the abvoe number of cells is from shuffle
% p95r(n)=length(find(shuffle95>real95r(n)))/nShuffle;%compared to shuffle distritbution: how different the abvoe number of cells is from shuffle
% real5r(n)=length(find(dr(:,n)<thresh5)); %number of cells below 5th percentile of combination
% real95r(n)=length(find(dr(:,n)>thresh95));%number of cells above 95th percentile of combination
% 
real5r(n)=length(find(dr(:,n)<thresh5))/length(dr(:,n)); %number of cells below 5th percentile of combination
real95r(n)=length(find(dr(:,n)>thresh95))/length(dr(:,n));%number of cells above 95th percentile of combination
real595r(n)=length(find(dr(:,n)>thresh5&dr(:,n)<thresh95))/length(dr(:,n)); 
realout595r(n)=length(find(dr(:,n)<thresh5 | dr(:,n)>thresh95))/length(dr(:,n)); 

end
real5=length(find(d<thresh5))/totald;
real95=length(find(d>thresh95))/totald;

real595=length(find(d>thresh5&d<thresh95))/totald;
realout595=length(find(d<thresh5 | d>thresh95))/totald;

%take the smallest peak: compare the smallest peak is lower than what we have 

rtogether=[];
rtogether(:,1)=real5r;
rtogether(:,2)=real95r;
[real5r,i5]=max(rtogether,[],2);
[real95r,i95]=min(rtogether,[],2);

%always plot the short peak on the right side
drr=[];
for n=1:size(ru,1);
    if i5(n)==1;
        drr(:,n)=dr(:,n);
    else
        drr(:,n)=-dr(:,n);
    end
end

%

figure,
edges = linspace(-1,1, 22 );
[N, ~] = histcounts(d, edges);
N = N / sum(N); N = N*100;
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'k', 'FaceAlpha', 0.5 );

Nr=[];
for n=1:size(ru,1);
[NS, ~] = histcounts(drr(:,n), edges);
NS = NS / sum(NS); NS = NS*100;
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'r', 'FaceAlpha', 0.5 );
Nr(n,:)=NS;
end
for n=1:size(Nr,1);
    hold on
    plot(edges(1:end-1) + diff(edges)/2,Nr(n,:),'Color',[0.8 0.8 0.8])
end
hold on
semshade(Nr,0.2,'k',edges(1:end-1) + diff(edges)/2,1);
semshade(N,0.2,'r',edges(1:end-1) + diff(edges)/2,1);

hold on
line([thresh5 thresh5],[0 11])
hold on
line([thresh95 thresh95],[0 11])

title('compareWithOtherCombination')

saveas(gcf,'Dfof_distributionOfRealCombinations.fig') %this figure actually has sem but since n shuffle =100, the sem looks very small


figure

subplot(221)
[x,p]=ksdensity(real5r);
plot(p,x)

hold on
% line([real5 real5],[0 max(x)])
% p5=length(find(real5r>real5))/size(ru,1);
line([real5 real5],[0 max(x)])
% p=length(find(real5r>real5))/size(ru,1);
[r,p]=ttest(real5r,real5);
title(['5 p=',num2str(p)])
xlim([0 0.35])

subplot(222)
[x,p]=ksdensity(real95r);
plot(p,x)

hold on
% line([real5 real5],[0 max(x)])
% p5=length(find(real5r>real5))/size(ru,1);
line([real95 real95],[0 max(x)])
% p=length(find(real95r>real95))/size(ru,1);
[r,p]=ttest(real95r,real95);

title(['95 p=',num2str(p)])


subplot(223)
[x,p]=ksdensity(real595r);
plot(p,x)

hold on
% line([real5 real5],[0 max(x)])
% p5=length(find(real5r>real5))/size(ru,1);
line([real595 real595],[0 max(x)])
p=length(find(real595r>real595))/size(ru,1);
[~,p]=ttest(real595r,real595);

title(['5-95 p=',num2str(p)])


subplot(224)
[x,p]=ksdensity(realout595r);
plot(p,x)

hold on
% line([real5 real5],[0 max(x)])
% p5=length(find(real5r>real5))/size(ru,1);
line([realout595 realout595],[0 max(x)])
p=length(find(realout595r>realout595))/size(ru,1);
[~,p]=ttest(realout595r,realout595);

title(['out5-95 p=',num2str(p)])


saveas(gcf,'Dfof_distributionOfRealShuffle_compareToCombinations.fig')

%%

% % %use the thredhold of these combinations
thresh5r=prctile(reshape(dr,[size(dr,1)*size(dr,2) 1]),5);
thresh95r=prctile(reshape(dr,[size(dr,1)*size(dr,2) 1]),95);

for n=1:size(ru,1)
real5r(n)=length(find(dr(:,n)<thresh5r)); %number of cells below 5th percentile of combination
real95r(n)=length(find(dr(:,n)>thresh95r));%number of cells above 95th percentile of combination
real595r(n)=length(find(dr(:,n)>thresh5r&dr(:,n)<thresh95r));
realout595r(n)=length(find(dr(:,n)<thresh5r | dr(:,n)>thresh95r));

end

real5rreal=length(find(d<thresh5r));
real95rreal=length(find(d>thresh95r));
real595rreal=length(find(d>thresh5r&d<thresh95r));
realout595rreal=length(find(d<thresh5r | d>thresh95r));

rtogether=[];
rtogether(:,1)=real5r;
rtogether(:,2)=real95r;
[real5r,i5]=max(rtogether,[],2);
[real95r,i95]=min(rtogether,[],2);

%always plot the short peak on the right side
drr=[];
for n=1:size(ru,1);
    if i5(n)==1;
        drr(:,n)=dr(:,n);
    else
        drr(:,n)=-dr(:,n);
    end
end


%

figure,
edges = linspace(-1,1, 22 );
[N, ~] = histcounts(d, edges);
N = N / sum(N); N = N*100;
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'k', 'FaceAlpha', 0.5 );

Nr=[];
for n=1:size(ru,1);
[NS, ~] = histcounts(drr(:,n), edges);
NS = NS / sum(NS); NS = NS*100;
% bar( edges(1:end-1) + diff(edges), N, 'FaceColor', 'r', 'FaceAlpha', 0.5 );
Nr(n,:)=NS;
end
semshade(N,0.2,'r',edges(1:end-1) + diff(edges),1);
for n=1:size(Nr,1);
    hold on
    plot(edges(1:end-1) + diff(edges),Nr(n,:),'Color',[0.8 0.8 0.8])
end
hold on
semshade(Nr,0.2,'k',edges(1:end-1) + diff(edges),1);
title('compareWithOtherCombination')

saveas(gcf,'Dfof_distributionOfRealCombinations_comShuffle.fig') 

figure

subplot(221)
[x,p]=ksdensity(real5r);
plot(p,x)

hold on
% line([real5 real5],[0 max(x)])
% p5=length(find(real5r>real5))/size(ru,1);
line([real5rreal real5rreal],[0 max(x)])
p=length(find(real5r>real5rreal))/size(ru,1);

title(['5 p=',num2str(p)])

subplot(222)
[x,p]=ksdensity(real95r);
plot(p,x)

hold on
% line([real5 real5],[0 max(x)])
% p5=length(find(real5r>real5))/size(ru,1);
line([real95rreal real95rreal],[0 max(x)])
p=length(find(real95r>real95rreal))/size(ru,1);

title(['95 p=',num2str(p)])

subplot(223)
[x,p]=ksdensity(real595r);
plot(p,x)

hold on
% line([real5 real5],[0 max(x)])
% p5=length(find(real5r>real5))/size(ru,1);
line([real595rreal real595rreal],[0 max(x)])
p=length(find(real595r>real595rreal))/size(ru,1);

title(['5-95 p=',num2str(p)])

subplot(224)
[x,p]=ksdensity(realout595r);
plot(p,x)

hold on
% line([real5 real5],[0 max(x)])
% p5=length(find(real5r>real5))/size(ru,1);
line([realout595rreal realout595rreal],[0 max(x)])
p=length(find(realout595r>realout595rreal))/size(ru,1);

title(['out5-95 p=',num2str(p)])

saveas(gcf,'Dfof_distributionOfRealShuffle_compareToCombinations_comShuffle.fig')


%% negative correlations between dfof in auditory and visual

% 10-bin, same as above
aoIdx = [29:38 43:52 63:72 111:120];
vsIdx = [9:18 81:90 95:104 127:136];
% a=zeros(1,140);
% a(aoIdx)=1;
% a(vsIdx)=1;
% figure,plot(a)


avr = fields(:, aoIdx); avr = nanmean(avr, 2);
vvr = fields(:, vsIdx); vvr = nanmean(vvr, 2);

%negative correlation between avr and vvr
% all=fields(:, [aoIdx vsIdx]);
% da=normalize(avr);
% dv=normalize(vvr);
all=nanmean(fields,2);
da=avr./all;
dv=vvr./all;
ia=find(~isnan(da));
iv=find(~isnan(dv));
iav=intersect(ia,iv);
figure
subplot(121)

plot(da(iav),dv(iav),'r.')
[r,p]=corr(da(iav),dv(iav));
title(['avr./all, p=',num2str(p),' r=',num2str(r)]);
%this correlation is positive


%ignore this one below because this is basically the above one - 1
% da=(avr-all)./all;
% dv=(vvr-all)./all;
% ia=find(~isnan(da));
% iv=find(~isnan(dv));
% iav=intersect(ia,iv);
% subplot(132)
% 
% plot(da(iav),dv(iav),'r.')
% [r,p]=corr(da(iav),dv(iav));
% title(['(avr-all)./all, p=',num2str(p),' r=',num2str(r)]);

da=(avr-all);
dv=(vvr-all);
ia=find(~isnan(da));
iv=find(~isnan(dv));
iav=intersect(ia,iv);
subplot(122)

plot(da(iav),dv(iav),'r.')
[rm,pm]=corr(da(iav),dv(iav));
title(['avr-all, p=',num2str(p),' r=',num2str(r)]);

saveas(gcf,'auditoryVisualDFOFCORR.fig')

%% other combinations

dar=[]; %avr./all;
dvr=[]; %vvr./all;

damr=[]; %(avr-all);
dvmr=[]; %(vvr-all);
for n=1:size(ru,1);
aoIdxr=cell2mat(zones(ru(n,1:4)));
vsIdxr=cell2mat(zones(ru(n,5:8)));
avr = fields(:, aoIdxr); avr = nanmean(avr, 2);
vvr = fields(:, vsIdxr); vvr = nanmean(vvr, 2);
dar(:,n)=avr./all;
dvr(:,n)=vvr./all;
damr(:,n)=avr./all;
dvmr(:,n)=vvr./all;

end

rr=[];
rp=[];
rmr=[];
rmp=[];
for n=1:size(ru,1);
    a=dar(:,n);
    v=dvr(:,n);
    ia=find(~isnan(a));
iv=find(~isnan(v));
iav=intersect(ia,iv);
    [rr(n,1),rp(n,1)]=corr(a(iav),v(iav));

    a=damr(:,n);
    v=dvmr(:,n);
    ia=find(~isnan(a));
iv=find(~isnan(v));
iav=intersect(ia,iv);
     [rmr(n,1),rmp(n,1)]=corr(a(iav),v(iav));
end

figure
subplot(221)
plot(rr,'r.')
hold on
line([0 size(ru,1)],[r r])
title('correlation,avr./all')

subplot(222)
plot(rp,'r.')
hold on
line([0 size(ru,1)],[p p])
title('p value,avr./all')

subplot(223)
plot(rmr,'r.')
hold on
line([0 size(ru,1)],[rm rm])
title('correlation,avr-all')

subplot(224)
plot(rmp,'r.')
hold on
line([0 size(ru,1)],[pm pm])
title('p value,avr-all')

saveas(gcf,'auditoryVisualDFOFCORR_compareToShuffle')

%% Absolute differnece 

figure,
ksdensity(abs(d))
hold on
ksdensity(abs(reshape(d2,[3333*100 1])))
hold on
ksdensity(abs(reshape(dr,[3333*34 1])))

title('absolute difference in real, complete shuffle, and other cue combinations')
saveas(gcf,'absoluteDifferences.fig')
