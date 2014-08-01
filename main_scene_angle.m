clear all; close all; clc;

addpath scripts

%hardcoded magnitudes for different subcarriers
%tx_sub_mag = [ 0.6294    0.8116   -0.7460    0.8268    0.2647   -0.8049   -0.4430    0.0938    0.9150    0.9298   -0.6848    0.9412    0.9143   -0.0292    0.6006 -0.7162];
%tx_sub_mag = [ 0.4964    0.6576    0.9421    0.7052    0.1084    1.0000    0.6185    0.1013    0.8010    0.8299    0.8719    0.8472    0.8426    0.2402    0.5338    0.6619];
tx_sub_mag = [0.121474674605681,0.327584561486124,0.0989478168282373,0.575789031655161,0.704429800943910,0.651381655427185,0.121981011878003,0.252700756842630,1,0.157742384332021,0.638787682419542,0.769284878645375,0.230583604352612,0.655406435184950,0.0594966708657510,0.174360489198279];

% System parameters
fs_scope = 5e9;                     % Tek sampling frequency
fs_tx = 1e9;                        % GageTx sampling frequency
Nsub = 16;                          % Number of sub-carriers in tx
CFAR = 1e-3;                        % Constant false alarm rate (for GLRT)

% Brute force parameters
numruns = 100;
%Mvec = 25:25:1000;
Mvec = 1:100;
subsToUse = 8:2:Nsub;

brick_profile = [];
brick_corner_profile = [];
drywall_profile = [];
drywall_corner_profile = [];
stand_profile = [];
corner_profile = [];
scenario = '2ft_angles_v4';
angles = 20:10:150;
%read in data for each scenario
for angle = 1:length(angles)
    %brick
    %read in I, Q, and combined data
    [data, i_data, q_data] = model_read_angle(scenario,'brick',angles(angle));
    %downsample each one, in the end, i only used the combined
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub, tx_sub_mag);
    [testspec_i, ~] = downsamp_to_subcarriers(i_data, fs_scope, fs_tx, Nsub, tx_sub_mag);
    [testspec_q, ~] = downsamp_to_subcarriers(q_data, fs_scope, fs_tx, Nsub, tx_sub_mag);
    %set the profile
    brick_profile(angle,:,:) = testspec;
    %plot a bunch of stuff i used to debug
    if(angles(angle) == 90)
        fre_prof = mean(testspec, 2);
        fre_std = std(testspec') ./ max(abs(fre_prof));
        fre_prof = fre_prof ./ max(abs(fre_prof));
        figure; hold on;
        mean_ = plot(fre_prof, '-bs', 'LineWidth',2);
        sd_ = plot(fre_prof-fre_std', 'r', 'LineWidth',2);
        plot(fre_prof+fre_std', 'r', 'LineWidth',2);
        sub_mag_ = plot(tx_sub_mag, '-gs', 'LineWidth', 2);
        title('Scenario: Brick');
        xlabel('Subcarrier No.');
        ylabel('Spectral Magnitude');
        %legend([mean_ sd_ sub_mag_ norm_ i_ q_], 'RX FP Mean', 'RX FP SD', 'TX FP',...
        %   'RX FP Normalized', 'RX I FP Mean', 'RX Q FP Mean');
        legend([mean_ sd_ sub_mag_], 'RX FP Mean', 'RX FP SD', 'TX FP');
        grid on;
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\brick_fp_w_sd')
        
        figure; hold on;
        plot(mean(i_data(:,100:500)), 'c', 'LineWidth', 2);
        plot(mean(q_data(:,100:500)), 'g', 'LineWidth', 2);
        plot(mean(data(:,100:500)), 'b', 'LineWidth', 2);
        grid on;
        xlabel('Sample No.');
        ylabel('Magnitude');
        title('Scenario: Brick');
        legend('I','Q','Combined');
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\brick_time')

        
    end
    %brick_profile = [brick_profile ; fre_prof'];
    
    %brick w/ reflector
    [data, i_data, q_data] = model_read_angle(scenario,'brick_corner',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub, tx_sub_mag);
    brick_corner_profile(angle,:,:) = testspec;
     if(angles(angle) == 90)
        fre_prof = mean(testspec, 2);
        fre_std = std(testspec') ./ max(abs(fre_prof));
        fre_prof = fre_prof ./ max(abs(fre_prof));
        figure; hold on;
        mean_ = plot(fre_prof, '-bs', 'LineWidth',2);
        sd_ = plot(fre_prof-fre_std', 'r', 'LineWidth',2);
        plot(fre_prof+fre_std', 'r', 'LineWidth',2);
        sub_mag_ = plot(tx_sub_mag, '-gs', 'LineWidth', 2);
        title('Scenario: Brick Corner');
        xlabel('Subcarrier No.');
        ylabel('Spectral Magnitude');
        legend([mean_ sd_ sub_mag_], 'RX FP Mean', 'RX FP SD', 'TX FP');
        grid on;
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\brick_corner_fp_w_sd')
        
        figure; hold on;
        plot(mean(i_data(:,100:500)), 'c', 'LineWidth', 2);
        plot(mean(q_data(:,100:500)), 'g', 'LineWidth', 2);
        plot(mean(data(:,100:500)), 'b', 'LineWidth', 2);
        grid on;
        xlabel('Sample No.');
        ylabel('Magnitude');
        title('Scenario: Brick Corner');
        legend('I','Q','Combined');
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\brick_corner_time')

     end
    
    %drywall
    [data, i_data, q_data] = model_read_angle(scenario,'drywall',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub, tx_sub_mag);
    drywall_profile(angle,:,:) = testspec;
    if(angles(angle) == 90)
        fre_prof = mean(testspec, 2);
        fre_std = std(testspec') ./ max(abs(fre_prof));
        fre_prof = fre_prof ./ max(abs(fre_prof));
        figure; hold on;
        mean_ = plot(fre_prof, '-bs', 'LineWidth',2);
        sd_ = plot(fre_prof-fre_std', 'r', 'LineWidth',2);
        plot(fre_prof+fre_std', 'r', 'LineWidth',2);
        sub_mag_ = plot(tx_sub_mag, '-gs', 'LineWidth', 2);
        title('Scenario: Drywall');
        xlabel('Subcarrier No.');
        ylabel('Spectral Magnitude');
        legend([mean_ sd_ sub_mag_], 'RX FP Mean', 'RX FP SD', 'TX FP');
        grid on;
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\drywall_fp_w_sd')
        
        figure; hold on;
        plot(mean(i_data(:,100:500)), 'c', 'LineWidth', 2);
        plot(mean(q_data(:,100:500)), 'g', 'LineWidth', 2);
        plot(mean(data(:,100:500)), 'b', 'LineWidth', 2);
        grid on;
        xlabel('Sample No.');
        ylabel('Magnitude');
        title('Scenario: Drywall');
        legend('I','Q','Combined');
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\drywall_time')

    end

    
    %drywall w/ reflector
    [data, i_data, q_data] = model_read_angle(scenario,'drywall_corner',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub, tx_sub_mag);
    drywall_corner_profile(angle,:,:) = testspec;
    if(angles(angle) == 90)
        fre_prof = mean(testspec, 2);
        fre_std = std(testspec') ./ max(abs(fre_prof));
        fre_prof = fre_prof ./ max(abs(fre_prof));
        figure; hold on;
        mean_ = plot(fre_prof, '-bs', 'LineWidth',2);
        sd_ = plot(fre_prof-fre_std', 'r', 'LineWidth',2);
        plot(fre_prof+fre_std', 'r', 'LineWidth',2);
        sub_mag_ = plot(tx_sub_mag, '-gs', 'LineWidth', 2);
        title('Scenario: Drywall Corner');
        xlabel('Subcarrier No.');
        ylabel('Spectral Magnitude');
        legend([mean_ sd_ sub_mag_], 'RX FP Mean', 'RX FP SD', 'TX FP');
        grid on;
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\drywall_corner_fp_w_sd')
        
        figure; hold on;
        plot(mean(i_data(:,100:500)), 'c', 'LineWidth', 2);
        plot(mean(q_data(:,100:500)), 'g', 'LineWidth', 2);
        plot(mean(data(:,100:500)), 'b', 'LineWidth', 2);
        grid on;
        xlabel('Sample No.');
        ylabel('Magnitude');
        title('Scenario: Drywall Corner');
        legend('I','Q','Combined');
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\drywall_corner_time')
    end
    
    %stand
    [data, i_data, q_data] = model_read_angle(scenario,'stand',0);
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub, tx_sub_mag);
    stand_profile(angle,:,:) = testspec;
    if(angles(angle) == 90)
        fre_prof = mean(testspec, 2);
        fre_std = std(testspec') ./ max(abs(fre_prof));
        fre_prof = fre_prof ./ max(abs(fre_prof));
        figure; hold on;
        mean_ = plot(fre_prof, '-bs', 'LineWidth',2);
        sd_ = plot(fre_prof-fre_std', 'r', 'LineWidth',2);
        plot(fre_prof+fre_std', 'r', 'LineWidth',2);
        sub_mag_ = plot(tx_sub_mag, '-gs', 'LineWidth', 2);
        title('Scenario: Stand');
        xlabel('Subcarrier No.');
        ylabel('Spectral Magnitude');
        legend([mean_ sd_ sub_mag_], 'RX FP Mean', 'RX FP SD', 'TX FP');
        grid on;
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\stand_fp_w_sd')
        
        figure; hold on;
        plot(mean(i_data(:,100:500)), 'c', 'LineWidth', 2);
        plot(mean(q_data(:,100:500)), 'g', 'LineWidth', 2);
        plot(mean(data(:,100:500)), 'b', 'LineWidth', 2);
        grid on;
        xlabel('Sample No.');
        ylabel('Magnitude');
        title('Scenario: Stand');
        legend('I','Q','Combined');
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\stand_time')
    end

    
    %stand w/ reflector
    [data, i_data, q_data] = model_read_angle(scenario,'corner',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub, tx_sub_mag);
    corner_profile(angle,:,:) = testspec;
    if(angles(angle) == 90)
        fre_prof = mean(testspec, 2);
        fre_std = std(testspec') ./ max(abs(fre_prof));
        fre_prof = fre_prof ./ max(abs(fre_prof));
        figure; hold on;
        mean_ = plot(fre_prof, '-bs', 'LineWidth',2);
        sd_ = plot(fre_prof-fre_std', 'r', 'LineWidth',2);
        plot(fre_prof+fre_std', 'r', 'LineWidth',2);
        sub_mag_ = plot(tx_sub_mag, '-gs', 'LineWidth', 2);
        title('Scenario: Corner');
        xlabel('Subcarrier No.');
        ylabel('Spectral Magnitude');
        legend([mean_ sd_ sub_mag_], 'RX FP Mean', 'RX FP SD', 'TX FP');
        grid on;
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\corner_fp_w_sd')
        
        figure; hold on;
        plot(mean(i_data(:,100:500)), 'c', 'LineWidth', 2);
        plot(mean(q_data(:,100:500)), 'g', 'LineWidth', 2);
        plot(mean(data(:,100:500)), 'b', 'LineWidth', 2);
        grid on;
        xlabel('Sample No.');
        ylabel('Magnitude');
        title('Scenario: Corner');
        legend('I','Q','Combined');
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\corner_time')
    end
    
    %drywall w/ plane
    [data, i_data, q_data] = model_read_angle(scenario,'drywall_plane',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub, tx_sub_mag);
    drywall_plane_profile(angle,:,:) = testspec;
    if(angles(angle) == 90)
        fre_prof = mean(testspec, 2);
        fre_std = std(testspec') ./ max(abs(fre_prof));
        fre_prof = fre_prof ./ max(abs(fre_prof));
        figure; hold on;
        mean_ = plot(fre_prof, '-bs', 'LineWidth',2);
        sd_ = plot(fre_prof-fre_std', 'r', 'LineWidth',2);
        plot(fre_prof+fre_std', 'r', 'LineWidth',2);
        sub_mag_ = plot(tx_sub_mag, '-gs', 'LineWidth', 2);
        title('Scenario: Drywall Plane');
        xlabel('Subcarrier No.');
        ylabel('Spectral Magnitude');
        legend([mean_ sd_ sub_mag_], 'RX FP Mean', 'RX FP SD', 'TX FP');
        grid on;
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\drywall_plane_fp_w_sd')
        
        figure; hold on;
        plot(mean(i_data(:,100:500)), 'c', 'LineWidth', 2);
        plot(mean(q_data(:,100:500)), 'g', 'LineWidth', 2);
        plot(mean(data(:,100:500)), 'b', 'LineWidth', 2);
        grid on;
        xlabel('Sample No.');
        ylabel('Magnitude');
        title('Scenario: Drywall Plane');
        legend('I','Q','Combined');
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\drywall_plane_time')
    end

    
    %drywall w/ bob
    [data, i_data, q_data] = model_read_angle(scenario,'drywall_human',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub, tx_sub_mag);
    drywall_human_profile(angle,:,:) = testspec;
    if(angles(angle) == 90)
        fre_prof = mean(testspec, 2);
        fre_std = std(testspec') ./ max(abs(fre_prof));
        fre_prof = fre_prof ./ max(abs(fre_prof));
        figure; hold on;
        mean_ = plot(fre_prof, '-bs', 'LineWidth',2);
        sd_ = plot(fre_prof-fre_std', 'r', 'LineWidth',2);
        plot(fre_prof+fre_std', 'r', 'LineWidth',2);
        sub_mag_ = plot(tx_sub_mag, '-gs', 'LineWidth', 2);
        title('Scenario: Drywall Human');
        xlabel('Subcarrier No.');
        ylabel('Spectral Magnitude');
        legend([mean_ sd_ sub_mag_], 'RX FP Mean', 'RX FP SD', 'TX FP');
        grid on;
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\drywall_human_fp_w_sd')
        
        figure; hold on;
        plot(mean(i_data(:,100:500)), 'c', 'LineWidth', 2);
        plot(mean(q_data(:,100:500)), 'g', 'LineWidth', 2);
        plot(mean(data(:,100:500)), 'b', 'LineWidth', 2);
        grid on;
        xlabel('Sample No.');
        ylabel('Magnitude');
        title('Scenario: Drywall Human');
        legend('I','Q','Combined');
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\drywall_human_time')
    end

    
    %drywall w/ cyl
    [data, i_data, q_data] = model_read_angle(scenario,'drywall_cyl',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub, tx_sub_mag);
    drywall_cyl_profile(angle,:,:) = testspec;
    if(angles(angle) == 90)
        fre_prof = mean(testspec, 2);
        fre_std = std(testspec') ./ max(abs(fre_prof));
        fre_prof = fre_prof ./ max(abs(fre_prof));
        figure; hold on;
        mean_ = plot(fre_prof, '-bs', 'LineWidth',2);
        sd_ = plot(fre_prof-fre_std', 'r', 'LineWidth',2);
        plot(fre_prof+fre_std', 'r', 'LineWidth',2);
        sub_mag_ = plot(tx_sub_mag, '-gs', 'LineWidth', 2);
        title('Scenario: Drywall Cyl');
        xlabel('Subcarrier No.');
        ylabel('Spectral Magnitude');
        legend([mean_ sd_ sub_mag_], 'RX FP Mean', 'RX FP SD', 'TX FP');
        grid on;
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\drywall_cyl_fp_w_sd')
        
        figure; hold on;
        plot(mean(i_data(:,100:500)), 'c', 'LineWidth', 2);
        plot(mean(q_data(:,100:500)), 'g', 'LineWidth', 2);
        plot(mean(data(:,100:500)), 'b', 'LineWidth', 2);
        grid on;
        xlabel('Sample No.');
        ylabel('Magnitude');
        title('Scenario: Drywall Cyl');
        legend('I','Q','Combined');
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\drywall_cyl_time')
    end

    
    %brick w/ cyl
    [data, i_data, q_data] = model_read_angle(scenario,'brick_cyl',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub, tx_sub_mag);
    brick_cyl_profile(angle,:,:) = testspec;
    if(angles(angle) == 90)
        fre_prof = mean(testspec, 2);
        fre_std = std(testspec') ./ max(abs(fre_prof));
        fre_prof = fre_prof ./ max(abs(fre_prof));
        figure; hold on;
        mean_ = plot(fre_prof, '-bs', 'LineWidth',2);
        sd_ = plot(fre_prof-fre_std', 'r', 'LineWidth',2);
        plot(fre_prof+fre_std', 'r', 'LineWidth',2);
        sub_mag_ = plot(tx_sub_mag, '-gs', 'LineWidth', 2);
        title('Scenario: Brick Cyl');
        xlabel('Subcarrier No.');
        ylabel('Spectral Magnitude');
        legend([mean_ sd_ sub_mag_], 'RX FP Mean', 'RX FP SD', 'TX FP');
        grid on;
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\brick_cyl_fp_w_sd')
        
        figure; hold on;
        plot(mean(i_data(:,100:500)), 'c', 'LineWidth', 2);
        plot(mean(q_data(:,100:500)), 'g', 'LineWidth', 2);
        plot(mean(data(:,100:500)), 'b', 'LineWidth', 2);
        grid on;
        xlabel('Sample No.');
        ylabel('Magnitude');
        title('Scenario: Brick Cyl');
        legend('I','Q','Combined');
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\brick_cyl_time')
    end

    
    %brick w/ plane
    [data, i_data, q_data] = model_read_angle(scenario,'brick_plane',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub, tx_sub_mag);
    brick_plane_profile(angle,:,:) = testspec;
    if(angles(angle) == 90)
        fre_prof = mean(testspec, 2);
        fre_std = std(testspec') ./ max(abs(fre_prof));
        fre_prof = fre_prof ./ max(abs(fre_prof));
        figure; hold on;
        mean_ = plot(fre_prof, '-bs', 'LineWidth',2);
        sd_ = plot(fre_prof-fre_std', 'r', 'LineWidth',2);
        plot(fre_prof+fre_std', 'r', 'LineWidth',2);
        sub_mag_ = plot(tx_sub_mag, '-gs', 'LineWidth', 2);
        title('Scenario: Brick Plane');
        xlabel('Subcarrier No.');
        ylabel('Spectral Magnitude');
        legend([mean_ sd_ sub_mag_], 'RX FP Mean', 'RX FP SD', 'TX FP');
        grid on;
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\brick_plane_fp_w_sd')
        
        figure; hold on;
        plot(mean(i_data(:,100:500)), 'c', 'LineWidth', 2);
        plot(mean(q_data(:,100:500)), 'g', 'LineWidth', 2);
        plot(mean(data(:,100:500)), 'b', 'LineWidth', 2);
        grid on;
        xlabel('Sample No.');
        ylabel('Magnitude');
        title('Scenario: Brick Plane');
        legend('I','Q','Combined');
        set(gcf,'PaperPositionMode','manual')
        set(gcf,'paperposition',[1 1 7 3])
        print(gcf,'-dpng','7_23\brick_plane_time')
    end

end

%params for this run, a lot of these are no longer used with "scene"
%detection
h0data = {brick_plane_profile; drywall_plane_profile; corner_profile;...
    drywall_cyl_profile; drywall_human_profile; drywall_corner_profile; brick_corner_profile};
h1data = brick_cyl_profile;
%h0data_titles = {'Drywall Cylinder'; 'Drywall Plane'; 'Brick Plane';...
    %'Drywall Corner'; 'Brick Corner'; 'Drywall Human'; 'Brick Cylinder'};
testFP = h1data;%%%%%
%h0_title = 'Drywall Angle-Frequency Profile';
h1_title = 'Brick Cyl Angle-Frequency Profile';
%live_data_title = 'Corner Example Data';%%%%%
pd_or_pfa = 'pd';%%%%%detection or false alarm, scene detection only uses pd
%fig1 = ['6_11\drywall_fp'];
fig2 = ['7_23\brick_cyl_fp' ];
%fig3 = ['7_23\live_corner_fp' ];%%%%%%%
fig4 = ['7_23\h0everything_h1brick_cyl_pd'];%%%%%%%
%fig5 = ['7_23\corner_td'];%%%%%%%%


tic;
%multiple runs for statistics
fpm_temp = zeros(numruns,length(Mvec));
angles_temp = zeros(numruns, length(Mvec));
wrong = zeros(7,1);
for run = 1:numruns
%for run = 100
    clc; toc; tic;
    disp([num2str(run) ' of ' num2str(numruns)]);
    
    %for different sizes of data
    for mind = 1:length(Mvec)
        disp(mind)
        M = Mvec(mind);
        
        for angle = 1:length(angles)
            %first create set of data for test FP
            %test_at_this_angle = sqeeze(testdata(angle,:,:));
            %testFP(angle,:) = mean(test_at_this_angle);
            
            %create h0 and h1 FPs
            for ii = 1:length(h0data)
                h0_at_this_angle{ii} = squeeze(h0data{ii}(angle,:,:));
            end
            h1_at_this_angle = squeeze(h1data(angle,:,:));
            
            %only use M pulses to create FPs
            [datalen numpulse] = size(h1_at_this_angle);
            inds = randperm(numpulse);
            inds = sort(inds(1:M));
            
            %grab M pulses
            for ii = 1:length(h0data)
                h0_at_this_angle{ii} = h0_at_this_angle{ii}(:,inds);
            end
            h1_at_this_angle = h1_at_this_angle(:,inds);
            
            %construct FPMs for h0 and h1
            for ii = 1:length(h0data)
                h0FP{ii}(angle,:) = mean(h0_at_this_angle{ii},2);
            end
            h1FP(angle,:) = mean(h1_at_this_angle,2);
        end
        
        %for each piece of live data, choose an array of angles
        target_is_there = 0;
        target_isnt_there = 0;
        correct_angles = 0;
        for live_data_no = 1:size(testFP,3)
            start_index = randi(size(testFP,1)-3);
            
            singleFP = testFP(start_index:start_index+3,:,live_data_no);
            
            %compute least squares for both
            for ii = 1:length(h0data)
                h0_least_squares(ii,:) = two_d_corr(singleFP, h0FP{ii});
            end
            h1_least_squares = two_d_corr(singleFP, h1FP);
            
            %first check whether h0 or h1 had the lower least squares
            if(min(min((h0_least_squares))) < min(h1_least_squares))
                %h0 is true
                target_isnt_there = target_isnt_there + 1;
                %debugging human drywall stuff
%                 if(M == 100 && run == 100)
%                     [~, ix_which_profile] = min(min(h0_least_squares'));
%                     [~, ix_which_angle] = min(min(h0_least_squares));
%                     wrong(ix_which_profile) = wrong(ix_which_profile) + 1;

%                     %%%%%%%%%%
%                     figure; hold on;
%                     subplot(2,2,1);
%                     surf(linspace(0,500,16), (angles(start_index:start_index+3))-90, singleFP);
%                     title('Drywall Human Live Data');
%                     xlabel('Frequency (MHz)');
%                     ylabel('Incidence angle (degrees)');
%                     zlabel('Spectral Magnitude');
%                     %%%%%%%%%%
% %                     figure;
% %                     subplot(3,2,2);
% %                     surf(linspace(0,500,16), (angles-90), h1FP);
% %                     title('Drywall Human Full Profile');
% %                     xlabel('Frequency (MHz)');
% %                     ylabel('Incidence angle (degrees)');
% %                     zlabel('Spectral Magnitude');
% %                     %%%%%%%%%
%                     subplot(2,2,2);
%                     surf(linspace(0,500,16), (angles(start_index:start_index+3))-90,...
%                         h1FP(start_index:start_index+3,:));
%                     title('Drywall Human FP data (should have matched)');
%                     xlabel('Frequency (MHz)');
%                     ylabel('Incidence angle (degrees)');
%                     zlabel('Spectral Magnitude');
%                     %%%%%%%%%%
%                     subplot(2,2,3);
%                     surf(linspace(0,500,16), (angles(ix_which_angle:ix_which_angle+3))-90,...
%                         h0FP{ix_which_profile}(ix_which_angle:ix_which_angle+3,:));
%                     title(['False Detection Data (' h0data_titles{ix_which_profile}...
%                         ', matched this instead)']);
%                     xlabel('Frequency (MHz)');
%                     ylabel('Incidence angle (degrees)');
%                     zlabel('Spectral Magnitude');
%                     %%%%%%%%%%%%
%                     subplot(2,2,4); hold on;
%                     h0_ =  plot(angles(1:11)-90,h0_least_squares(ix_which_profile,:),...
%                         '-bs','linewidth',2,'markersize',8);
%                     h1_ = plot(angles(1:11)-90,h1_least_squares,'-rs','linewidth',2,'markersize',8);
%                     correct_ = plot(angles(start_index)-90, h1_least_squares(start_index), '-g*','linewidth',2,'markersize',16);
%                     false_ = plot(angles(ix_which_angle)-90, h0_least_squares(ix_which_profile,ix_which_angle), '-c*','linewidth',2,'markersize',16);
%                     xlabel('Incidence angle (degrees)');
%                     ylabel('Least Squares Magnitude');
%                     title(['Least Squares']);
%                      grid on;
% 
%                     %set(gcf,'units','normalized','outerposition',[0 0 1 1])
% 
%                     set(gcf,'PaperPositionMode','manual')
%                     set(gcf,'paperposition',[1 1 13 19])%7 3
%                     legger = legend([h0_, h1_, correct_, false_], ['h0(' h0data_titles{ix_which_profile}...
%                         ')'], 'h1 (Human Drywall)', 'Correct Match', 'Incorrect Match', 'location','NorthEast'); 
%                     fig_name = ['6_17\h0everything_h1drywall_human_wrong_guess_' num2str(target_isnt_there)];
%                     print(gcf,'-dpng',fig_name);
%                     
%                 end
%                 %done debugging human drywall stuff
            else
                %h1 is true
                target_is_there = target_is_there + 1;
                %check if correct angle found
                [val, ix] = min(h1_least_squares);
                if(ix == start_index)
                    correct_angles = correct_angles + 1;
                end
            end
        end
        fpm_temp(run, mind) = target_is_there / size(testFP,3);
        angles_temp(run, mind) = correct_angles / size(testFP,3);
    end
end
fpm_out = mean(fpm_temp);
angles_out = mean(angles_temp);

%Angle Frequency profiles for h0 and h1
% figure;
% surf(linspace(0,500,16),((angles)-90),h0FP{ii});
% title(h0_title);
% xlabel('Frequency (MHz)');
% ylabel('Incidence angle (degrees)');
% zlabel('Spectral Magnitude');
%
% set(gcf,'PaperPositionMode','manual')
% set(gcf,'paperposition',[1 1 7 3])
% print(gcf,'-dpng',fig1)
%
figure;
surf(linspace(0,500,16),((angles)-90),h1FP);
title(h1_title);
xlabel('Frequency (MHz)');
ylabel('Incidence angle (degrees)');
zlabel('Spectral Magnitude');

set(gcf,'PaperPositionMode','manual')
set(gcf,'paperposition',[1 1 7 3])
print(gcf,'-dpng',fig2)

%Example "live" data
% figure;
% surf(linspace(0,500,16), (angles(start_index:start_index+3))-90, singleFP);
% title(live_data_title);
% xlabel('Frequency (MHz)');
% ylabel('Incidence angle (degrees)');
% zlabel('Spectral Magnitude');
% 
% set(gcf,'PaperPositionMode','manual')
% set(gcf,'paperposition',[1 1 7 3])
% print(gcf,'-dpng',fig3)

%styles = {'-r*','-b*','-k*','-m*','-g*','-rs','-bs','-ks','-ms','-gs'};
%legger{length(subsToUse)+1} = 'FPM';
%
figure('Name','Detection Results')
plot(Mvec,fpm_out,'-bs','linewidth',2,'markersize',8)
xlabel('Model Size');
if(strcmp(pd_or_pfa,'pd'))
    ylabel('P(correct detection)')
else
    ylabel('P(false alarm)')
end

set(gcf,'PaperPositionMode','manual')
set(gcf,'paperposition',[1 1 7 3])
print(gcf,'-dpng',fig4)


%
% ylabel('P(false alarm)')
%legend(legger,'location','best')

% figure('Name','Detection Angle Results')
% plot(Mvec,angles_out,'-bs','linewidth',2,'markersize',8)
% xlabel('Model Size');
% ylabel('P(correct detection angle')
%legend(legger,'location','best')



%
% glrt_temp = zeros(length(Mvec),length(subsToUse),numruns);
% fpm_temp = zeros(numruns,length(Mvec));
% tic;
%
% for run = 1:numruns
%     clc; toc; tic;
%     disp([num2str(run) ' of ' num2str(numruns)]);
%
%
%     for mind = 1:length(Mvec)
%         disp(mind)
%         M = Mvec(mind);             % How many test pulses to use in models
%         [h0spec, ~] = ...
%             downsamp_to_subcarriers(modelRead(baseID,R,M), fs_scope, fs_tx, Nsub);
%         [h1spec resampfn] = ...
%             downsamp_to_subcarriers(modelRead(testID,R,M), fs_scope, fs_tx, Nsub);
%         h0FP = mean(h0spec,2);
%         h1FP = mean(h1spec,2);
%
%
%
% %         figure;
% %         subplot(2,1,1);
% %         plot(h0FP);
% %         subplot(2,1,2);
% %         plot(h1FP);
% %         pause(.5);
%         % FPM DETECTION
%         h0error = zeros(1,numtest);
%         h1error = zeros(1,numtest);
%         for i = 1:numtest
%             h0error(i) = sum((testspec(:,i)-h0FP).^2);
%             h1error(i) = sum((testspec(:,i)-h1FP).^2);
%         end; clear i;
%         %h0error
%         %h1error
%         FPM = h1error < h0error;
%         %FPM
%
%         fpm_temp(run,mind) = sum(FPM)/numtest;
%
%
%
%
%         % GLRT USING ACTUAL RETURN PDFs OBTAINED FROM KSDENSITY
%         numer = zeros(Nsub,numtest);
%         denom = zeros(Nsub,numtest);
%         gamma = zeros(1,Nsub);
%         for subnum = 1:Nsub
%             % Use KSDENSITY() to get PDFs and CDFs
%             [h0pdf h0x] = ksdensity(h0spec(subnum,:),'npoints',1000);
%             h0dx = h0x(2)-h0x(1);
%             h0cdf = cumsum(h0pdf)*h0dx;
%
%             [h1pdf h1x] = ksdensity(h1spec(subnum,:),'npoints',1000);
%             h1dx = h1x(2)-h1x(1);
%             h1cdf = cumsum(h1pdf)*h1dx;
%
%             % Find where the H0 CDF is closest to (1-CFAR)
%             threshInd = find(min(abs((1-h0cdf)-CFAR)) == abs((1-h0cdf)-CFAR));
%             thresh = h0x(threshInd);
%
%             % Evaluate H0 and H1 PDFs at the threshold
%             h0thresh = pdfVal(h0pdf,h0x,thresh);
%             h1thresh = pdfVal(h1pdf,h1x,thresh);
%
%             % Per-subcarrier threshold
%             gamma(subnum) = h0thresh/h1thresh;
%
%
%             % Evaluate H0 and H1 for each of the test pulses.
%             for i = 1:numtest
%                 testval = testspec(subnum,i);
%                 numer(subnum,i) = pdfVal(h0pdf,h0x,testval);
%                 denom(subnum,i) = pdfVal(h1pdf,h1x,testval);
%             end; clear i;
%
%
%         end; clear subnum;
%
%         LRvals = numer./denom;      % LR for each of the sub-carriers and test pulses
%
%
%
%         for subind = 1:length(subsToUse)
%             % Select the best < Nsub sub-carriers to use
%             goodsubs = findBestSubcarriers_binary(h0spec,h1spec,subsToUse(subind),Nsub);
%
%             % Composite threshold and LRs for only the selected
%             % sub-carriers
%             GAMMA = prod(gamma(goodsubs == 1));
%             LR = zeros(1,numtest);
%             for ind = 1:numtest
%                 LRtemp = LRvals(:,ind)';
%                 LR(ind) = prod(LRtemp(goodsubs == 1));
%             end
%
%             GLRT = LR < GAMMA;              % Actual decisions (1 = correct)
%
%             glrt_temp(mind,subind,run) = sum(GLRT)/numtest;
%         end
%
%     end
%
% end
%
% % Average the results across the realizations
% glrt_out = mean(glrt_temp,3);
% fpm_out = mean(fpm_temp);
%
%
%
% % PLOTTING
% styles = {'-r*','-b*','-k*','-m*','-g*','-rs','-bs','-ks','-ms','-gs'};
% legger{length(subsToUse)+1} = 'FPM';
%
% figure('Name','Detection Results')
% for i = 1:length(subsToUse)
%     plot(Mvec,glrt_out(:,i),styles{i},'linewidth',2,'markersize',8);
%     legger{i} = ['GLRT: ' num2str(subsToUse(i)) ' sub-carriers'];
%     hold on
% end
% plot(Mvec,fpm_out,'-bs','linewidth',2,'markersize',8)
% xlabel('Model Size');
%
% ylabel('P(false alarm)')
% %ylabel('P(correct detection)')
% legend(legger,'location','best')
%
%





rmpath scripts










