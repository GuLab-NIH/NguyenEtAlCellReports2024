
% load data
% contact Lead Author for data


%% PLOT %
% S11A-B

indices = [14 9];
ttls = {'1-1', '2-2'};
figure; hold on
for ii = 1:numel(indices)
    switch indices(ii)
        case {1, 5, 6}
            clr = 'm';
            ttl = 'avr - ';
        case {11, 15, 16}
            clr = 'g';
            ttl = 'vvr - ';
        case {9, 14}
            clr = 'k';
            ttl = 'avr v.s. vvr - ';
    end
    xx = finalOverlap{indices(ii)}*100;
    
    subplot(1,2,ii); hold on
    for kk = 1:6
        switch kk
            case 1
                offset = -0.25;
                clr = 'k';
            case 2
                offset = -0.15;
            case 3
                offset = -0.05;
                clr = 'm';
            case 4
                offset = 0.05;
            case 5
                offset = 0.15;
                clr = 'g';
            case 6
                offset = 0.25;
        end
        yy = xx(:,kk);
        yy = yy( yy < mean(yy)+3*std(yy) & yy > mean(yy)-3*std(yy) );
        violin2(yy, 'x', kk, 'facecolor', clr, 'facealpha', 0.6-0.3*mod(kk,2));
        
        set(gca, 'xtick', 1:6, 'xticklabel', {'%P in C', '%S in C', '%P in unique1', '%S in unique1', '%P in unique2', '%S in unique2'}, ...
            'xticklabelrotation', 45);
        xlim([0.3 6.7]);
    end
    title([ttl ttls{ii}]);
    
end
set(gcf, 'Position', [111 111 1111 456]);







%% test everything in avr-vvr



if 1
indices = [14 9];
ttls = {'avr1 - vvr1', 'avr2 - vvr2'};
figure; hold on
for ii = 1:numel(indices)
    xx = finalOverlap{indices(ii)}*100;
    
    subplot(1,2,ii); hold on
    for kk = 1:6
        switch kk
            case {1,2}
                clr = 'k';
            case {3,4}
                clr = 'm';
            case {5,6}
                clr = 'g';
        end
        yy = xx(:,kk);
        violin2(yy, 'x', kk, 'facecolor', clr, 'facealpha', 0.3);        
        set(gca, 'xtick', 1:6, 'xticklabel', {'%P in C', '%S in C', '%P in avr', '%S in avr', '%P in vvr', '%S in vvr'}, ...
            'xticklabelrotation', 45);
    end
    title(ttls{ii});
    
end
set(gcf, 'Position', [111 111 777 345]);

end




% avr1-vvr1
disp('AVR1-VVR1');
xx = finalOverlap{14}*100;
for aa = 1:6
for bb = 1:6
    if aa >= bb
        continue
    end
    yy = xx(:,aa);
    zz = xx(:,bb);
    
    [~,p] = ttest(yy,zz);
    disp([num2str(aa) ' vs. ' num2str(bb) ': ' num2str(p)]);
end
end



% avr2-vvr2
disp('AVR2-VVR2');
xx = finalOverlap{9}*100;
for aa = 1:6
for bb = 1:6
    if aa >= bb
        continue
    end
    yy = xx(:,aa);
    zz = xx(:,bb);
    
    [~,p] = ttest(yy,zz);
    disp([num2str(aa) ' vs. ' num2str(bb) ': ' num2str(p)]);
end
end






%%
