% Breakout of the old main.m to just be the lateral case. I got sick of
% checking inputs, etc, so I split the lateral (door scan) apart from the
% rotate (detection/SAR) process.

clear all; close all; clc;
addpath scripts


c = 299792458;                      % Speed of light, m/s


fs_scope = 5e9;
ts_scope = 1/fs_scope;
fs_tx = 1e9;
Nsub = 16;


delX = 1;
maxDoor = 9;

% MOVE ACROSS THE CLUTTERED GRAD OFFICE DOOR
% scenario = '2013-05-23-1100';           % 12in steps
% scenario = '2013-05-23-1130';           % 18in steps
% scenario = '2013-05-23-1200';           % 6in steps
% scenario = '2013-05-24-1040';           % 12in steps, 6in start offset



% REDBLADE LAIR 2013-05-27
scenario = 'redblade_r38in_6in';
% scenario = 'redblade_r38in_12in';
% scenario = 'redblade_r38in_18in';
% scenario = 'redblade_r38in_12in_o6in';


% scenario = 'redblade_r26in_6in';
% scenario = 'redblade_r26in_12in';
% scenario = 'redblade_r26in_18in';
% scenario = 'redblade_r26in_12in_o6in';


% scenario = 'redblade_r12in_6in';
% scenario = 'redblade_r12in_12in';
% scenario = 'redblade_r12in_18in';
% scenario = 'redblade_r12in_12in_o6in';



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


% figure('Name','Range vs. Position')
% plot(ranges)
% xlabel('Position Index, \DeltaX')
% ylabel('Range, m')


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



figure('Name','Doorway Classification Results')
hold on
plot_truth = line([truth(1) truth(1)],[min(Adat) max(Adat)], ...
    'linewidth',4,'color','r');
line([truth(2) truth(2)],[min(Adat) max(Adat)], ...
    'linewidth',4,'color','r');
plot_exp = line([currentPos-bestWidth/2 currentPos-bestWidth/2],[min(Adat) max(Adat)], ...
    'linewidth',3,'color','b','linestyle','--');
line([currentPos+bestWidth/2 currentPos+bestWidth/2],[min(Adat) max(Adat)], ...
    'linewidth',3,'color','b','linestyle','--');
plot_A = plot(1:numpos,Adat,'-k*','linewidth',2,'markersize',10);

legend([plot_A,plot_truth,plot_exp], ...
    {'A','True Doorway Edges','Experimental Doorway Edges'},'location','best')
xlabel('Position Index, $$\Delta X$$','interpreter','latex')
ylabel('A')
% title(regexprep(scenario,'_','\\_'))

if exist('images\lateral','dir') ~= 7
    mkdir('images\lateral');
end
set(gcf,'PaperPositionMode','manual')
set(gcf,'paperposition',[1 1 7 3])
print(gcf,'-depsc',['images\lateral\' regexprep(scenario,'\W','')])
print(gcf,'-dpng',['images\lateral\' regexprep(scenario,'\W','')])



rmpath scripts


