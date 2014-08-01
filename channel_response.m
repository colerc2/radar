close all; clear all; clc;
addpath scripts;

Nsub = 16;
fs_tx = 1e9;
fs_scope = 5e9;
which_signals = 42;

%steps i need to do
%1.  read in time domain channel response
% .  read in daters
%2.  calculate FPs for specified, resampled, and real data with time-domain
%channel response subtracted and not subtracted
%3.  calculate fft of time-domain channel response Ho(f)
%4.  calculate fft of each received time-domain signal H(f)
%5.  H(f) / Ho(f)
%6.  convert back to time-domain
%7.  calculate FP of new signal
%8.  compare to the time-domain subtraction method

%Step 1
channel_response_in_time = load('data\different_subcarriers\subcarrier_testing\i_data_all_zeros_2.txt');
%ambient_signal = mean(ambient_signal);

%Step 2
%find FPs of specified signals
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
    tx_profiles(ii,1) = 0;%TODO WHY AM I HERE
    tx_profiles(ii,:) = tx_profiles(ii,:) ./ max(abs(tx_profiles(ii,:)));
    tx_profiles(ii,:) = abs(tx_profiles(ii,:)); 
end
%find FPs for resampled signals
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

%find transfer function of channel response
for ii = 1:size(channel_response_in_time,1)
   temp_sig = channel_response_in_time(ii,:) - mean(channel_response_in_time(ii,:));
   %try zero padding
   filler = zeros(1,5001-length(temp_sig));
   temp_sig = [temp_sig filler];
   channel_fft(ii,:) = (fft(temp_sig));
end
channel_fft = mean(channel_fft);

%find FPs for RX signals
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
        %if(subtract_ambient == 1)
         %   temp_sig = rx_signal(jj,:) - ambient_signal;
        %else
            temp_sig = rx_signal(jj,100:end);
        %end
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
    rx_profiles(ii,1) = 0;
    rx_profiles(ii,:) = rx_profiles(ii,:) ./ max(abs(rx_profiles(ii,:)));
    
    %and again with time-domain channel response subtracted
    for jj = 1:num_data
        temp_sig = rx_signal(jj,:) - mean(channel_response_in_time);
        data_length = length(temp_sig);
        
        temp_sig = temp_sig - mean(temp_sig);
        
        if(data_length < 5001)
            filler = zeros(1,5001-data_length);
            temp_sig = [temp_sig filler];
        end
        
        rx_filled(jj,:) = temp_sig;
    end
    
    [fp, ~] = downsamp_to_subcarriers(rx_filled, fs_scope, fs_tx, Nsub, 0);
    
    rx_profiles_channel_time(ii,:) = mean(fp,2);
    %rx_profiles_channel_time(ii,1) = 0;
    rx_profiles_channel_time(ii,:) = rx_profiles_channel_time(ii,:) ./ ...
        max(abs(rx_profiles_channel_time(ii,:)));
    
    %and again with frequency-domain channel response "subtracted"
    %(division)
    for jj = 1:num_data
        %first, find fft of signal
        temp_zero_mean_rx = rx_signal(jj,:) - mean(rx_signal(jj,:));
        %try zero padding
        filler = zeros(1,5001-length(temp_zero_mean_rx));
        temp_zero_mean_rx = [temp_zero_mean_rx filler];
        fft_rx_data = (fft(temp_zero_mean_rx));
        
        %divide rx by fft of channel response
        fft_after_division = (fft_rx_data - channel_fft);
        
%         figure; hold on;
%         plot(20*log10(fftshift(abs(fft_rx_data)/max(abs(fft_rx_data)))),'r');
%         plot(20*log10(fftshift(abs(channel_fft)/max(abs(channel_fft)))),'b');
%         plot(20*log10(fftshift(abs(fft_after_division)/max(abs(fft_after_division)))),'g');
%         pause;
        
        %ifft her to get back to time domain
        temp_sig = abs(ifft((fft_after_division)));
        
%         figure; hold on;
%         plot(abs(temp_sig), 'b','LineWidth', 2);
%         plot(rx_signal(jj,:), 'r', 'LineWidth', 2);
%         grid on;
%         legend('Modified', 'OG');
%         pause;
        
        data_length = length(temp_sig);
        
        temp_sig = temp_sig - mean(temp_sig);
        
        if(data_length < 5001)
            filler = zeros(1,5001-data_length);
            temp_sig = [temp_sig filler];
        end
        
        rx_filled(jj,:) = temp_sig;
    end
    
    [fp, ~] = downsamp_to_subcarriers(rx_filled, fs_scope, fs_tx, Nsub, 0);
    
    rx_profiles_channel_fft(ii,:) = mean(fp,2);
    rx_profiles_channel_fft(ii,1) = 0;
    rx_profiles_channel_fft(ii,:) = rx_profiles_channel_fft(ii,:) ./ ...
        max(abs(rx_profiles_channel_fft(ii,:)));
end





%plot FPs (only positive frequencies)
for ii = which_signals
   figure(ii);hold on;
   plot(tx_profiles(ii,2:17),'-bs', 'LineWidth', 2);
   plot(tx_resampled_profiles(ii,:), '-cs', 'LineWidth', 2);
   plot(rx_profiles(ii,:), '-rs', 'LineWidth', 2);
   plot(rx_profiles_channel_time(ii,:), '-ks', 'LineWidth', 2);
   plot(rx_profiles_channel_fft(ii,:), '-gs', 'LineWidth', 2);
   grid on;
   legend('Specified', 'Resampled', 'Real Data', 'Time Subtraction',...
       'Frequency Sub.');
   pause;

   %save figure
%    set(gcf,'PaperPositionMode','manual')
%    set(gcf,'paperposition',[1 1 7 3])
%    print(gcf,'-dpng',['6_30\without_first_sub_' num2str(ii)]);
end

