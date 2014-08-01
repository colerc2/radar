% Fork of mainLateral.m to show classification performance for all ranges
% and a single delX

clear all; close all; clc;
addpath scripts


c = 299792458;                      % Speed of light, m/s


fs_scope = 5e9;
ts_scope = 1/fs_scope;
fs_tx = 1e9;
Nsub = 16;


delX = 1;
maxDoor = 9;



plotname = 'redblade_6in_step';
scens = {'redblade_r12in_6in', 'redblade_r26in_6in', 'redblade_r38in_6in'};

% plotname = 'redblade_12in_step';
% scens = {'redblade_r12in_12in', 'redblade_r26in_12in', 'redblade_r38in_12in'};

% plotname = 'redblade_18in_step';
% scens = {'redblade_r12in_18in', 'redblade_r26in_18in', 'redblade_r38in_18in'};

% plotname = 'redblade_12in_o6in_step';
% scens = {'redblade_r12in_12in_o6in', 'redblade_r26in_12in_o6in', 'redblade_r38in_12in_o6in'};


for scennum = 1:length(scens)
    scenario = scens{scennum};
    
    numpos = calcNumPos(scenario);
    
    
    % Read in the ambient data and downsample it to the spectral magnitudes at
    % the specific sub-carrier frequencies
    ambdat = ambientRead;
    [ambSpec resampfn] = downsamp_to_subcarriers(ambdat,fs_scope,fs_tx,Nsub);
    
    
    tx = load(['data\lateral\' scenario '\tx_signal.txt']);
    
    avgSpec = zeros(Nsub,numpos);
    subPower = zeros(Nsub,numpos);
    ind = zeros(1,numpos);
    for pos = 1:numpos
        fprintf('Calculating Range-to-Target for Position %d of %d\n', ...
            pos, numpos);
        
        % Load data from pre-parsed tekdat files
        [dat truth] = lateralRead(scenario,pos);
        
        
        % Get spectral contribution of each sub-carrier
        [datSpec resampfn] = downsamp_to_subcarriers(dat,fs_scope,fs_tx,Nsub);
        
        
        % Mean spectral magnitude for each frequency/position
        avgSpec(:,pos) = mean(datSpec,2);
        
        
        % Calculate the average sub-carrier power
        subPower(:,pos) = (mean(datSpec,2)').^2;
        
        % Calculate range-to-environment
        corrdat = tekMF2(dat,tx,fs_scope,fs_tx,Nsub);
        [m, allind] = max(corrdat,[],2);
        ind(pos) = mean(allind);
        
        clc
    end
    
    ranges = c*(ts_scope*ind)/2;
    
    
    
    % Method of Moments -> Doorway ID and Classification.
    noiseVar = var(ambSpec,0,2)';
    
    A = zeros(numpos,Nsub);
    for pos = 1:numpos
        A(pos,:) = subPower(:,pos)' - noiseVar;
    end
    
    
    
    % Using MoM results (A), estimate the location and position:
    %       Ported from IEEE USA paper -> radarcon entrance (I think)
    Adat = mean(A,2);
    Adat = Adat - mean(Adat);
    
    numRef = floor(maxDoor / delX);         % The maximum number of steps we will be in the "door"
    bestWidth = 0; currentMax = 0; currentPos = 0;
    
    
    for i = 1:numRef
        refSig = [zeros(1,length(Adat)-i), -1*ones(1,i)];
        door_corr = xcorr(Adat,refSig);
        
        
        if (max(door_corr) > currentMax)
            bestWidth = (i+1)*delX;         % i+1 because putting a -1 at 3 indices means a square pulse that covers 2 actual delta-X's
            currentMax = max(door_corr);
            
            maxInd = find(door_corr == currentMax);
            currentPos = (maxInd - ((i-1)/2)) * delX;       % ((i-1)/2) to get the # of del-X's that the ref covers, divided by 2... to get the center
        end
    end
    
    fprintf('Door is %d steps wide centered at position %.1f\n', bestWidth, currentPos);
    
    Aout(scennum,:) = Adat';
    widthOut(scennum) = bestWidth;
    posOut(scennum) = currentPos;
end

figure('Name','Doorway Classification Results')
hold on
grid
minA = min(min(Aout));
maxA = max(max(Aout));

% True door edges
plot_truth = line([truth(1) truth(1)],[1.1*minA 1.1*maxA], 'linewidth',5,'color','k');
line([truth(2) truth(2)],[1.1*minA 1.1*maxA], 'linewidth',5,'color','k');

% Experimental Door Edges
plot_exp3 = line([posOut(3)-widthOut(3)/2 posOut(3)-widthOut(3)/2],[minA maxA], 'linewidth',3,'color','g','linestyle','-');
line([posOut(3)+widthOut(3)/2 posOut(3)+widthOut(3)/2],[minA maxA], 'linewidth',3,'color','g','linestyle','-');

plot_exp2 = line([posOut(2)-widthOut(2)/2 posOut(2)-widthOut(2)/2],[.9*minA .9*maxA], 'linewidth',3,'color','r','linestyle','-');
line([posOut(2)+widthOut(2)/2 posOut(2)+widthOut(2)/2],[.9*minA .9*maxA], 'linewidth',3,'color','r','linestyle','-');

plot_exp1 = line([posOut(1)-widthOut(1)/2 posOut(1)-widthOut(1)/2],[.8*minA .8*maxA], 'linewidth',3,'color','b','linestyle','-');
line([posOut(1)+widthOut(1)/2 posOut(1)+widthOut(1)/2],[.8*minA .8*maxA], 'linewidth',3,'color','b','linestyle','-');




% Plots of A
plot_A1 = plot(1:numpos,Aout(1,:),'-b*','linewidth',2,'markersize',10);
plot_A2 = plot(1:numpos,Aout(2,:),'-rs','linewidth',2,'markersize',10);
plot_A3 = plot(1:numpos,Aout(3,:),'-gx','linewidth',2,'markersize',10);





legend([plot_A1,plot_A2,plot_A3,plot_truth,plot_exp1,plot_exp2,plot_exp3], ...
    {'$$A: R=12in$$','$$A: R=26in$$','$$A: R=38in$$','True Doorway Edges','$$Edges: R=12in$$','$$Edges: R=26in$$','$$Edges: R=38in$$'},'location','EastOutside','interpreter','latex')
xlabel('Position Index, $$\Delta X$$','interpreter','latex')
ylabel('A')
% title(regexprep(scenario,'_','\\_'))

if exist('images\lateral','dir') ~= 7
    mkdir('images\lateral');
end
set(gcf,'PaperPositionMode','manual')
set(gcf,'paperposition',[1 1 7 4])
print(gcf,'-depsc',['images\lateral\' plotname])
print(gcf,'-dpng',['images\lateral\' plotname])



rmpath scripts


