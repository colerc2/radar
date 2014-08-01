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
%scenario = 'redblade_r38in_6in';
scenario = 'r12in_6in';

%numpos = calcNumPos(scenario);
numpos = 21; %should be 21???



% Read in the ambient data and downsample it to the spectral magnitudes at
% the specific sub-carrier frequencies
%ambdat = ambientRead;
%[ambSpec resampfn] = downsamp_to_subcarriers(ambdat,fs_scope,fs_tx,Nsub);

%tx = load(['data\lateral\' scenario '\tx_signal.txt']);
tx = load(['data\lateral\redblade_r38in_6in\tx_signal.txt']);

avgSpec = zeros(Nsub,numpos);
subPower = zeros(Nsub,numpos);
ind = zeros(1,numpos);
fpm_true = zeros(1,numpos);
for pos = 1:numpos
    fprintf('Calculating Range-to-Target for Position %d of %d\n', ...
        pos, numpos);
 
    %FPM stuff
    testdat_fpm = modelRead('12in',pos-1);
    [numtest_fpm, ~] = size(testdat_fpm);
    [testspec_fpm, ~] = downsamp_to_subcarriers(testdat_fpm, fs_scope, fs_tx, Nsub);
    subsToUse_fpm = 8:2:Nsub;
    scenario_fpm = 'data\door_12in_v2\i_data_12in_drywall';
    scenario_fpm_2 = 'data\door_12in_v2\q_data_12in_drywall';
    h0_fpm = modelReadBob(scenario_fpm, scenario_fpm_2);
    scenario_fpm = 'data\door_12in_v2\i_data_12in_drywall_corner';
    scenario_fpm_2 = 'data\door_12in_v2\q_data_12in_drywall_corner';
    h1_fpm = modelReadBob(scenario_fpm, scenario_fpm_2);
    
    [h0spec, ~] = downsamp_to_subcarriers(h0_fpm, fs_scope, fs_tx, Nsub);
    [h1spec, ~] = downsamp_to_subcarriers(h1_fpm, fs_scope, fs_tx, Nsub);
     h0FP = mean(h0spec,2);
     h1FP = mean(h1spec,2);
     
%      figure;
%      subplot(3,1,1);
%      plot(h0FP);
%      subplot(3,1,2);
%      plot(h1FP);
%      subplot(3,1,3);
%      plot(testspec_fpm);
     %pause;
     h0error = zeros(1,numtest_fpm);
        h1error = zeros(1,numtest_fpm);
        for i = 1:numtest_fpm     
            %sub = testspec_fpm(:,1)-max(testspec_fpm(:,1));
            %h0error(i) = sum((sub-(h0FP-max(h0FP))).^2);
            %h1error(i) = sum((sub-(h1FP-max(h1FP))).^2);
            h0error(i) = sum((testspec_fpm(:,i)-h0FP).^2);
            h1error(i) = sum((testspec_fpm(:,i)-h1FP).^2);
        end; clear i;
        %h0error
        %h1error
        FPM = h1error < h0error;
        sum(FPM)
        sum(h0error)
        sum(h1error)

        
    fpm_at_pos(pos) = (sum(FPM)>80);
    
    % Load data from pre-parsed tekdat files
    %[dat truth] = lateralRead(scenario,pos);
    ifile = ['data\door_12in_v2\i_data_12in_' num2str(pos-1) '.txt'];
    qfile = ['data\door_12in_v2\q_data_12in_' num2str(pos-1) '.txt'];
    idat = load(ifile);
    qdat = load(qfile);
    [numpulse datlen] = size(idat);
    
    % Combine I and Q
    dat = sqrt(idat.^2 + qdat.^2);
    
    % Zero-mean each signal... NOT SURE IF NEED
    for i = 1:numpulse
        dat(i,:) = dat(i,:) - mean(dat(i,:));
    end
    
    % MAKE IT 5001 samples long. NOT SURE IF NEED OR WANT
    if datlen < 5001
        filler = zeros(numpulse,5001-datlen);
        dat = [dat filler];
    elseif datlen > 5001
        dat = dat(:,1:5001);
    end
    
    
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

truth = [8 14];

ranges = c*(ts_scope*ind)/2;


% figure('Name','Range vs. Position')
% plot(ranges)
% xlabel('Position Index, \DeltaX')
% ylabel('Range, m')


% Method of Moments -> Doorway ID and Classification.
%noiseVar = var(ambSpec,0,2)';

A = zeros(numpos,Nsub);
for pos = 1:numpos
    A(pos,:) = subPower(:,pos)';% - noiseVar;
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

target = 18;

figure('Name','Doorway Classification Results')
hold on
plot_truth = line([truth(1) truth(1)],[min(Adat) max(Adat)], ...
    'linewidth',6,'color','r');
line([truth(2) truth(2)],[min(Adat) max(Adat)], ...
    'linewidth',6,'color','r');
plot_exp = line([currentPos-bestWidth/2 currentPos-bestWidth/2],[min(Adat) max(Adat)], ...
    'linewidth',3,'color','g','linestyle','-.');
line([currentPos+bestWidth/2 currentPos+bestWidth/2],[min(Adat) max(Adat)], ...
    'linewidth',3,'color','g','linestyle','-.');
plot_true_target = line([target target],[min(Adat) max(Adat)], ...
    'linewidth',6,'color','b');
for ii = 1:length(fpm_at_pos)
   if(fpm_at_pos(ii) == 1)
       plot_ex_target = line([ii ii],[min(Adat) max(Adat)], ...
    'linewidth',3,'color','c','linestyle','-.');
   end
end
plot_A = plot(1:numpos,Adat,'-k*','linewidth',2,'markersize',10);
set(gca,'fontsize',16);

legend([plot_A,plot_truth,plot_exp,plot_true_target,plot_ex_target], ...
    {'A','True Doorway Edges','Experimental Doorway Edges','True Target Position',...
    'Experimental Target Position'},'location','best','fontsize',16)
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

% figure;
% plot(fpm_at_pos);

rmpath scripts


