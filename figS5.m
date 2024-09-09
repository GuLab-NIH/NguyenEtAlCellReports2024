
%% cue scores
% figS5

load('allFolders.mat');
thresholds = nan(1,5);
for VR = 1:5
    clear folders;
    switch VR
        case 1  % VAVR
            folders = folders_combo;
            clr = 'c';
            
        case 2  % AVR1
            folders = folders_audio;
            clr = 'm';
            
        case 3  % VVR1
            folders = folders_visual;
            clr = 'g';
            
        case 4  % AVR2
            folders = folders_audio_NEW_learned;
            clr = 'm';
            
        case 5  % VVR2
            folders = folders_visual_NEW_learned;
            clr = 'g';
    end
    
    realScores = [];
    shuffledScores = [];
    for ff = 1:numel(folders)
        disp(['VR: ' num2str(VR) ' - ' num2str(ff) '/' num2str(numel(folders))]);
        
        cd(folders{ff});
        load('corrIncorr_20230701\cueCells_20240619.mat');
        
        realScores = [realScores cueCells.realScores];
        shuffledScores = [shuffledScores reshape(cueCells.shuffleScores, 1, [])];
        
    end
    
    thresholds(VR) = prctile(shuffledScores, 95);
    
    figure; hold on
    
    [f, xi] = ksdensity(realScores, 'bandwidth', .1);
    plot(xi, f, clr, 'linestyle', '-');
    [f, xi] = ksdensity(shuffledScores, 'bandwidth', .1);
    plot(xi, f, clr, 'linestyle', '--');
    
    [~,p] = kstest2(realScores, shuffledScores);
    title(p);
    set(gcf, 'position', [234*VR-111 111 234 222]);
    
end


%%
