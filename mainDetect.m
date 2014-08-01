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
baseID = 'brick';                   % H0
%testID = 'brick_corner';            % H1
%baseID = 'stand_gage';
%testID = 'corner_gage';
%baseID = 'drywall';
testID = 'drywall_corner';
%baseID = 'drywall_corner';
%testID = 'brick_corner';
%baseID = 'corner_perm';
%testID = 'drywall_corner_perm';
%R = 100;                            % IN CENTIMETERS
R = 100;

% System parameters
fs_scope = 5e9;                     % Tek sampling frequency
fs_tx = 1e9;                        % GageTx sampling frequency
Nsub = 16;                          % Number of sub-carriers in tx
CFAR = 1e-3;                        % Constant false alarm rate (for GLRT)


% Load all pulses we wish to form detection decisions for
%testdat = modelRead(testID,R);%test probability of correct detection
testdat = modelRead(baseID,R);%test probability of false alarm
%testdat = modelRead('brick',R);
%testdat = modelRead('nothing',R);
%testdat = modelRead('drywall',R);
%testdat = modelRead('drywall_corner',R);
%testdat = modelRead('brick_gage_card',R);
[numtest, ~] = size(testdat);
[testspec, ~] = downsamp_to_subcarriers(testdat, fs_scope, fs_tx, Nsub);


% Brute force parameters
numruns = 100;
%Mvec = 25:25:1000;
%Mvec = 20:10:numtest;
Mvec = 100;
subsToUse = 8:2:Nsub;


glrt_temp = zeros(length(Mvec),length(subsToUse),numruns);
fpm_temp = zeros(numruns,length(Mvec));
tic;

for run = 1:numruns
    clc; toc; tic;
    disp([num2str(run) ' of ' num2str(numruns)]);
   

    for mind = 1:length(Mvec)
        disp(mind)
        M = Mvec(mind);             % How many test pulses to use in models
        [h0spec, ~] = ...
            downsamp_to_subcarriers(modelRead(baseID,R,M), fs_scope, fs_tx, Nsub);
        [h1spec resampfn] = ...
            downsamp_to_subcarriers(modelRead(testID,R,M), fs_scope, fs_tx, Nsub);
        h0FP = mean(h0spec,2);
        h1FP = mean(h1spec,2);
        
        fuck = 1;
        
        
%         figure;
%         subplot(2,1,1);
%         plot(h0FP);
%         subplot(2,1,2);
%         plot(h1FP);
%         pause(.5);
        % FPM DETECTION
        h0error = zeros(1,numtest);
        h1error = zeros(1,numtest);
        for i = 1:numtest
            h0error(i) = sum((testspec(:,i)-h0FP).^2);
            h1error(i) = sum((testspec(:,i)-h1FP).^2);
        end; clear i;
        %h0error
        %h1error
        FPM = h1error < h0error;
        %FPM
        
        fpm_temp(run,mind) = sum(FPM)/numtest;
        
        
        
        
        % GLRT USING ACTUAL RETURN PDFs OBTAINED FROM KSDENSITY
        numer = zeros(Nsub,numtest);
        denom = zeros(Nsub,numtest);
        gamma = zeros(1,Nsub);
        for subnum = 1:Nsub
            % Use KSDENSITY() to get PDFs and CDFs
            [h0pdf h0x] = ksdensity(h0spec(subnum,:),'npoints',1000);
            h0dx = h0x(2)-h0x(1);
            h0cdf = cumsum(h0pdf)*h0dx;
            
            [h1pdf h1x] = ksdensity(h1spec(subnum,:),'npoints',1000);
            h1dx = h1x(2)-h1x(1);
            h1cdf = cumsum(h1pdf)*h1dx;
            
            % Find where the H0 CDF is closest to (1-CFAR)
            threshInd = find(min(abs((1-h0cdf)-CFAR)) == abs((1-h0cdf)-CFAR));
            thresh = h0x(threshInd);
            
            % Evaluate H0 and H1 PDFs at the threshold
            h0thresh = pdfVal(h0pdf,h0x,thresh);
            h1thresh = pdfVal(h1pdf,h1x,thresh);
            
            % Per-subcarrier threshold
            gamma(subnum) = h0thresh/h1thresh;
            
            
            % Evaluate H0 and H1 for each of the test pulses.
            for i = 1:numtest
                testval = testspec(subnum,i);
                numer(subnum,i) = pdfVal(h0pdf,h0x,testval);
                denom(subnum,i) = pdfVal(h1pdf,h1x,testval);
            end; clear i;
            
            
        end; clear subnum;
        
        LRvals = numer./denom;      % LR for each of the sub-carriers and test pulses
        
        
        
        for subind = 1:length(subsToUse)
            % Select the best < Nsub sub-carriers to use
            goodsubs = findBestSubcarriers_binary(h0spec,h1spec,subsToUse(subind),Nsub);
            
            % Composite threshold and LRs for only the selected
            % sub-carriers
            GAMMA = prod(gamma(goodsubs == 1));
            LR = zeros(1,numtest);
            for ind = 1:numtest
                LRtemp = LRvals(:,ind)';
                LR(ind) = prod(LRtemp(goodsubs == 1));
            end
            
            GLRT = LR < GAMMA;              % Actual decisions (1 = correct)
            
            glrt_temp(mind,subind,run) = sum(GLRT)/numtest;
        end
        
    end
    
end

% Average the results across the realizations
glrt_out = mean(glrt_temp,3);
fpm_out = mean(fpm_temp);



% PLOTTING
styles = {'-r*','-b*','-k*','-m*','-g*','-rs','-bs','-ks','-ms','-gs'};
legger{length(subsToUse)+1} = 'FPM';

figure('Name','Detection Results')
for i = 1:length(subsToUse)
    plot(Mvec,glrt_out(:,i),styles{i},'linewidth',2,'markersize',8);
    legger{i} = ['GLRT: ' num2str(subsToUse(i)) ' sub-carriers'];
    hold on
end
plot(Mvec,fpm_out,'-bs','linewidth',2,'markersize',8)
xlabel('Model Size');

ylabel('P(false alarm)')
%ylabel('P(correct detection)')
legend(legger,'location','best')







rmpath scripts










