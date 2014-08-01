close all; clear all; clc;
addpath scripts;

Nsub = 16;
fs_tx = 1e9;
fs_scope = 5e9;
which_signals = 37:78;
which_signals = [42, 46, 47];
subtract_ambient = 1;

%read in ambient signal
ambient_signal = load('data\different_subcarriers\subcarrier_testing\i_data_all_zeros_2.txt');
ambient_signal = mean(ambient_signal);



%read in 20 different TX signals
for ii = which_signals
    tx_signal(ii,:) = load(['data\different_subcarriers\tx_signal_' num2str(ii) ...
        '.txt'], 'r');
    %normalize
    tx_signal(ii,:) = tx_signal(ii,:) ./ max(abs(tx_signal(ii,:)));
    
    %make zero mean
    tx_signal(ii,:) = tx_signal(ii,:) - mean(tx_signal(ii,:));
    
    %find FP
    tx_profiles(ii,:) = fft(tx_signal(ii,1:end-1));
    
    %normalize
    %tx_profiles(ii,1) = 0;
    tx_profiles(ii,:) = tx_profiles(ii,:) ./ max(abs(tx_profiles(ii,:)));
    tx_profiles(ii,:) = abs(tx_profiles(ii,:));
    
end

%resample all 20 TX signals to 5GS/s
timescale = fs_scope/fs_tx;         % Resample multiplier
txlen = 2*Nsub + 1;                 % Number of actual samples in transmit waveform

for ii = which_signals
    % Resample the transmitted signal to match the 5Gsample/sec of the scope
    refdat = tx_signal(ii,1:txlen);
    refdat = resample(refdat,timescale,1);
    
    % Remove DC offset
    refdat = refdat - mean(refdat);
    
    %make it 5001 samples long
    filler = zeros(1,5001-length(refdat));
    tx_signal_resampled(ii,:) = [refdat filler];
    tx_signal_resampled_copy(ii,:) = refdat;
    
end

tx_resampled_profiles = downsamp_to_subcarriers(tx_signal_resampled, fs_scope, fs_tx, Nsub, 0);
tx_resampled_profiles = tx_resampled_profiles';

for ii = which_signals
    tx_resampled_profiles(ii,:) = tx_resampled_profiles(ii,:) ./ ...
        max(abs(tx_resampled_profiles(ii,:)));
end


%read in all the different signals RX from the scope
for ii = which_signals
    rx_signal = load(['data\different_subcarriers\subcarrier_testing\i_data_signal_' ...
        num2str(ii) '.txt']);
    
    [num_data, ~] = size(rx_signal);
    
    %optional plot of time domain data
    %     time_domain_average = mean(rx_signal);
    %     time_domain_average = time_domain_average(100:end);
    %     figure; hold on;
    %     plot(time_domain_average, 'LineWidth', 2);
    %     grid on; pause;
    
    for jj = 1:num_data
        if(subtract_ambient == 1)
            temp_sig = rx_signal(jj,:) - ambient_signal;
        else
            temp_sig = rx_signal(jj,100:end);
        end
        data_length = length(temp_sig);
        
        temp_sig = temp_sig - mean(temp_sig);
        
        if(data_length < 5001)
            filler = zeros(1,5001-data_length);
            temp_sig = [temp_sig filler];
        end
        
        rx_filled(jj,:) = temp_sig;
    end
    
    [fp, ~] = downsamp_to_subcarriers(rx_filled, fs_scope, fs_tx, Nsub, 0);
    
    rx_profiles(ii,:) = mean(fp,2);
    %rx_profiles(ii,1) = 0;
    rx_profiles(ii,:) = rx_profiles(ii,:) ./ max(abs(rx_profiles(ii,:)));
end

%plot FPs (only positive frequencies)
for ii = which_signals
   figure(ii);hold on;
   plot(tx_profiles(ii,2:17),'-bs', 'LineWidth', 2);
   plot(rx_profiles(ii,:), '-rs', 'LineWidth', 2);
   plot(tx_resampled_profiles(ii,:), '-cs', 'LineWidth', 2);
   grid on;
   legend('Specified', 'Real Data', 'Resampled');
   pause;

   %save figure
%    set(gcf,'PaperPositionMode','manual')
%    set(gcf,'paperposition',[1 1 7 3])
%    print(gcf,'-dpng',['6_30\without_first_sub_' num2str(ii)]);
end


%save certain signals
% for ii = which_signals
%     tx_signal = load(['data\different_subcarriers\tx_signal_' num2str(ii) ...
%         '.txt'], 'r');
%     
%     
%     fid = fopen(['good_signals\tx_signal_' num2str(ii) ...
%         '.txt'], 'w');
%     for jj = 1:400
%        fprintf(fid, '%i, ', 0); 
%     end
%     for jj = 1:length(tx_signal)
%         fprintf(fid, '%i, ', tx_signal(jj));
%     end
%     fprintf(fid, '%i', 0);
%     
% end

%find average error on each subcarrier
% for ii = which_signals
% %     diff(ii,:) = (abs(tx_profiles(ii,2:17) - rx_profiles(ii,:)));
%     diff(ii,:) = (abs(tx_resampled_profiles(ii,:) - rx_profiles(ii,:)));
% end
% error = (mean(diff(37:end,:)));
% figure; hold on;
% plot(error, '-bs', 'LineWidth', 2);
% grid on;
% xlabel('Subcarrier No.');
% ylabel('Average error');
% set(gcf,'PaperPositionMode','manual')
% set(gcf,'paperposition',[1 1 7 3])
% print(gcf,'-dpng',['7_3\subcarrier_error']);

%Peak to average power calculations
for ii = which_signals
    %papr for 1GS/s signal
    max_squared = (max(abs(tx_signal(ii,:)))) ^ 2;
    average_squared = mean(tx_signal(ii,:) .^ 2);
    papr(ii) = max_squared / average_squared;
    
    %papr for 5GS/s signal
    max_squared = (max(abs(tx_signal_resampled_copy(ii,:)))) ^ 2;
    average_squared = mean(tx_signal_resampled_copy(ii,:) .^ 2);
    papr_resampled(ii) = max_squared / average_squared;
end

%average error
% error = mean(diff,2);
% figure; hold on; grid on;
% scatter(papr(37:end), error(37:end), 'b');
% scatter(papr_resampled(37:end), error(37:end), 'r');
% xlabel('PAPR'); ylabel('Average Error');
% legend('1 GS/s', '5 GS/s');
% set(gcf,'PaperPositionMode','manual')
% set(gcf,'paperposition',[1 1 7 3])
% print(gcf,'-dpng',['7_3\PAPR_vs_error']);

%find ideal signals
% num_signals_to_find = 3;
% found = 0;
% [sorted, ix] = sort(error);
% for ii = 1:length(ix)
%     if(ix(ii) > 36)%valid signal
%        papr(ix(ii))
%        found = found + 1;
%        figure; hold on; title(['PAPR = ' num2str(papr(ix(ii)))]);
%        plot(tx_profiles(ix(ii),2:17),'-bs', 'LineWidth', 2);
%        plot(rx_profiles(ix(ii),:), '-rs', 'LineWidth', 2);
%        plot(tx_resampled_profiles(ix(ii),:), '-cs', 'LineWidth', 2);
%        xlabel('Subcarrier No.');
%        ylabel('Magnitude');
%        grid on;
%        legend('Specified', 'Real Data', 'Resampled');
%        set(gcf,'PaperPositionMode','manual')
%        set(gcf,'paperposition',[1 1 7 3])
%        print(gcf,'-dpng',['7_3\random_good_signal_' num2str(found)]);
%        if(num_signals_to_find == found)
%           break; 
%        end
%     end
% end



