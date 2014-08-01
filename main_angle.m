% A rework/fork of the previously-implemented detection methods. This was
% completely rewritten (again) to make use of the "final" data format being
% output from the Tek scope using Bob's code. This lets everything be
% significantly easier/faster to do because the data is already in
% CSV/plaintext.
%
% Additionally, the Weibull MLE has been replaced with KSDENSITY(), an
% inbuilt Matlab function that generates a discrete PDF of the input data.
% These PDFs have significantly more freedom, and tend to follow the data
% in a more realistic fashion. A linear interpolation function was also
% written so that PDF values could be estimated from these new
% non-closed-form discrete PDFs.
%
% Semantics:
%       "base" data is the null hypothesis. In this case, it is the wall's
%       response alone
%
%       "test" data is the H1 hypothesis, the "wall + reflector". It is
%       also the set of data that will be fed into each model



clear all; close all; clc;

addpath scripts

% Set file/scenario IDs and the corresponding range
%baseID = 'brick';                   % H0
%testID = 'brick_corner';            % H1
%baseID = 'stand_gage';
%testID = 'corner_gage';
%baseID = 'drywall';
%testID = 'drywall_corner';
%baseID = 'drywall_corner';
%testID = 'brick_corner';
%baseID = 'corner_perm';
%testID = 'drywall_corner_perm';
%R = 100;                            % IN CENTIMETERS
R = 50;

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
scenario = '2ft_angles_v3';
angles = 20:10:150;
for angle = 1:length(angles)
    %brick
    data = model_read_angle(scenario,'brick',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub);
    brick_profile(angle,:,:) = testspec;
    %fre_prof = mean(testspec, 2);
    %brick_profile = [brick_profile ; fre_prof'];
    
    %brick w/ reflector
    data = model_read_angle(scenario,'brick_corner',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub);
    brick_corner_profile(angle,:,:) = testspec;
    %fre_prof = mean(testspec, 2);
    %brick_corner_profile = [brick_corner_profile ; fre_prof'];
    
    %drywall
    data = model_read_angle(scenario,'drywall',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub);
    drywall_profile(angle,:,:) = testspec;
    %fre_prof = mean(testspec, 2);
    %drywall_profile = [drywall_profile ; fre_prof'];
    
    %drywall w/ reflector
    data = model_read_angle(scenario,'drywall_corner',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub);
    drywall_corner_profile(angle,:,:) = testspec;
    %fre_prof = mean(testspec, 2);
    %drywall_corner_profile = [drywall_corner_profile ; fre_prof'];
    
    %stand
    data = model_read_angle(scenario,'stand',0);
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub);
    stand_profile(angle,:,:) = testspec;
    %fre_prof = mean(testspec, 2);
    %stand_profile = [stand_profile ; fre_prof'];
    
    %stand w/ reflector
    data = model_read_angle(scenario,'corner',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub);
    corner_profile(angle,:,:) = testspec;
    %fre_prof = mean(testspec, 2);
    %corner_profile = [corner_profile ; fre_prof'];
    
    %drywall w/ plane
    data = model_read_angle(scenario,'drywall_plane',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub);
    drywall_plane_profile(angle,:,:) = testspec;
    
    %drywall w/ bob
    data = model_read_angle(scenario,'drywall_human',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub);
    drywall_human_profile(angle,:,:) = testspec;
    
    %drywall w/ cyl
    data = model_read_angle(scenario,'drywall_cyl',angles(angle));
    [testspec, ~] = downsamp_to_subcarriers(data, fs_scope, fs_tx, Nsub);
    drywall_cyl_profile(angle,:,:) = testspec;
end

h0data = drywall_plane_profile;
h1data = drywall_cyl_profile;
testFP = h0data;%%%%%
h0_title = 'Drywall Angle-Frequency Profile';
h1_title = 'Drywall Human Corner Angle-Frequency Profile';
live_data_title = 'Drywall Human Example Data';%%%%%
pd_or_pfa = 'pfa';%%%%%
fig1 = ['6_11\drywall_fp'];
fig2 = ['6_11\drywall_human_fp' ];
fig3 = ['6_11\live_drywall_human_fp' ];%%%%%%%
fig4 = ['6_11\h0drywall_plane_h1drywall_cyl_pfa'];%%%%%%%


tic;
%multiple runs for statistics
fpm_temp = zeros(numruns,length(Mvec));
angles_temp = zeros(numruns, length(Mvec));
for run = 1:numruns
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
           h0_at_this_angle = squeeze(h0data(angle,:,:));
           h1_at_this_angle = squeeze(h1data(angle,:,:));
           
           %only use M pulses to create FPs
           [datalen numpulse] = size(h0_at_this_angle);
           inds = randperm(numpulse);
           inds = sort(inds(1:M));
           
           %grab M pulses
           h0_at_this_angle = h0_at_this_angle(:,inds);
           h1_at_this_angle = h1_at_this_angle(:,inds);
           
           %construct FPMs for h0 and h1
           h0FP(angle,:) = mean(h0_at_this_angle,2);
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
           h0_least_squares = two_d_corr(singleFP, h0FP);
           h1_least_squares = two_d_corr(singleFP, h1FP);
           
           %first check whether h0 or h1 had the lower least squares
           if(min(h0_least_squares) < min(h1_least_squares))
              %h0 is true
              target_isnt_there = target_isnt_there + 1;
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
figure;
surf(linspace(0,500,16),((angles)-90),h0FP);
title(h0_title);
xlabel('Frequency (MHz)');
ylabel('Incidence angle (degrees)');
zlabel('Spectral Magnitude');

set(gcf,'PaperPositionMode','manual')
set(gcf,'paperposition',[1 1 7 3])
print(gcf,'-dpng',fig1)

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
figure;
surf(linspace(0,500,16), (angles(start_index:start_index+3))-90, singleFP);
title(live_data_title);
xlabel('Frequency (MHz)');
ylabel('Incidence angle (degrees)');
zlabel('Spectral Magnitude');

set(gcf,'PaperPositionMode','manual')
set(gcf,'paperposition',[1 1 7 3])
print(gcf,'-dpng',fig3)

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










