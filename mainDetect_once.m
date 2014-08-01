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
testID = 'brick_corner';            % H1
R = 100;                            % IN CENTIMETERS

% System parameters
fs_scope = 5e9;                     % Tek sampling frequency
fs_tx = 1e9;                        % GageTx sampling frequency
Nsub = 16;                          % Number of sub-carriers in tx
CFAR = 1e-3;                        % Constant false alarm rate (for GLRT)



% Load all pulses we wish to form detection decisions for
testdat = modelRead(testID,R);
[numtest, ~] = size(testdat);
[testspec, ~] = downsamp_to_subcarriers(testdat, fs_scope, fs_tx, Nsub);


% Generate average frequency response (FREQUENCY PROFILES) for each case,
% using M training pulses in the actual FP
M = 100;                            % How many test pulses to use in models
[h0spec, ~] = ...
    downsamp_to_subcarriers(modelRead(baseID,R,M), fs_scope, fs_tx, Nsub);
[h1spec resampfn] = ...
    downsamp_to_subcarriers(modelRead(testID,R,M), fs_scope, fs_tx, Nsub);
h0FP = mean(h0spec,2);
h1FP = mean(h1spec,2);



% figure('Name','Frequency Profiles')
% plot(resampfn, h0FP, '-b*')
% hold on
% plot(resampfn, h1FP, '-rx');
% legend(['H0: ' baseID], ['H1: ' testID])



% FPM DETECTION
h0error = zeros(1,numtest);
h1error = zeros(1,numtest);
for i = 1:numtest
    h0error(i) = sum((testspec(:,i)-h0FP).^2);
    h1error(i) = sum((testspec(:,i)-h1FP).^2);
end; clear i;
FPM = h1error < h0error;


% allpdf = figure('Name','All PDFs with CFAR Thresholds');
% pdfcdf = figure('Name','PDF and CDF');
somepdf = figure('Name','Some PDFs with CFAR Thresholds');

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
    
    
    % Plot PDF and CDFs, including threshold and evaluated points.
%     figure(pdfcdf)
%     subplot(211)
%     plot(h0x,h0pdf,'-b'); hold on; plot(h1x,h1pdf,'-r')
%     line([thresh thresh], [0 1.2*max([h0pdf h1pdf])],'color','k')
%     plot(thresh,h0thresh,'bs')
%     plot(thresh,h1thresh,'rs'); hold off;
%     ylabel('PDF(x)')
%     legend('H_0','H_1','CFAR Threshold')
%     title(['Subcarrier ' num2str(subnum) ': ' num2str(round(resampfn(subnum)/1e6)) ' MHz'])
%     
%     subplot(212)
%     plot(h0x,h0cdf,'-b'); hold on; plot(h1x,h1cdf,'-r'); hold off;
%     xlabel('x'); ylabel('CDF(x)')
    
    
%     figure(allpdf)
%     subplot(4,4,subnum);
%     plot(h0x,h0pdf,'-b');
%     hold on;
%     plot(h1x,h1pdf,'-r');
%     line([thresh thresh], [0 1.2*max([h0pdf h1pdf])],'color','k')
%     plot(thresh,h0thresh,'-bs')
%     plot(thresh,h1thresh,'-rs')
%     title([num2str(round(resampfn(subnum)/1e6)) ' MHz'])
%     if subnum == 13
%         ylabel('PDF(x)')
%         xlabel('x')
%     elseif subnum == 16
%         legend('H_0','H_1','CFAR Threshold','H_0(\gamma)','H_1(\gamma)')
%     end
    
    if subnum > 4 && subnum < 9
        
        figure(somepdf)
        subplot(2,2,subnum-4)
        plot(h0x,h0pdf,'-b','linewidth',2)
        hold on
        plot(h1x,h1pdf,'-r','linewidth',2)
        line([thresh thresh], [0 1.2*max([h0pdf h1pdf])],'color','k','linewidth',2)
        plot(thresh,h0thresh,'bs','linewidth',2)
        plot(thresh,h1thresh,'rs','linewidth',2)
        title(['Sub-carrier ' num2str(subnum)])
        grid
        if subnum == 7
            ylabel('PDF(x)')
            xlabel('x')
        elseif subnum == 8
            legend({'$$H_0$$','$$H_1$$','$$\gamma$$','$$H_0 ( \gamma )$$','$$H_1 ( \gamma )$$'},'interpreter','latex')
        end
    end
    
    
%     pause;
    
    
end; clear subnum;

LRvals = numer./denom;      % LR for each of the sub-carriers and test pulses

% Composite threshold and LRs
GAMMA = prod(gamma);
LR = prod(LRvals);

GLRT = LR < GAMMA;


fprintf('GLRT Correct: %d/%d\nFPM Correct: %d/%d\n', ...
    sum(GLRT), numtest, sum(FPM), numtest);







rmpath scripts










