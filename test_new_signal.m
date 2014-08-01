close all; clear all; clc;
addpath scripts;

Nsub = 16;
fs_tx = 1e9;
fs_scope = 5e9;
which_signals = [42 46 47];

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
    rx_signal = load(['good_signals\good_subcarriers\i_data_signal_' ...
        num2str(ii) '.txt']);
    
    [num_data, ~] = size(rx_signal);
    
    %optional plot of time domain data
        time_domain_average = mean(rx_signal);
        time_domain_average = time_domain_average(100:end);
        figure; hold on;
        real_plot = time_domain_average(95:265);
        res_plot = tx_signal_resampled_copy(ii,:);
        real_plot = real_plot ./ max(abs(real_plot));
        res_plot = res_plot ./ max(abs(res_plot));
        plot(real_plot, 'b', 'LineWidth', 2);
        plot(res_plot, 'r' , 'LineWidth', 2);
        grid on; pause;
    
    for jj = 1:num_data
        temp_sig = rx_signal(jj,:);

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
