% Isolate the subcarrier magnitudes at each of the Nsub sub-carriers. 

function [resamp resampfn] = downsamp_to_subcarriers(dat, fs_scope, fs_tx, Nsub, tx_magnitudes)

fs_scale = fs_scope/fs_tx;


[num_dats datlen] = size(dat);


resamp = zeros(Nsub,num_dats);
for i = 1:num_dats
    
    
    vec = dat(i,:);
    RX = fft(vec);
    
    posRX = abs(RX(1:1+(length(RX)-1)/2));
    
    
    N = min(floor(length(posRX)/fs_scale)+1,length(posRX));%changed to min in case fs_scale = 1
    posfn = linspace(0,fs_tx*(N-1)/N/2,N);
    
    posRX = posRX(1:min(N,length(posRX)));%changed from posRX(1:N)
    
    indices = 1:Nsub;
    indices = N/(Nsub+.5) + (indices-1)*(N/(Nsub+.5));
    %indices_ = indices;
    indices = ceil(indices);
    
    resampfn = posfn(indices);
    
    %linear interpolate
%     fract = indices_ - floor(indices_);
%     lower = floor(indices_);
%     upper = ceil(indices_);
%     diff = posRX(upper)-posRX(lower);
%     resamp(:,i) = posRX(lower) + diff.*fract;
    
    resamp(:,i) = posRX(indices);
    
    %normalize by TX subcarrier magnitudes
    
    %resamp(:,i) = resamp(:,i) ./ max(abs(resamp(:,i)));
    %tx_magnitudes = tx_magnitudes ./max(abs(tx_magnitudes));
    
     %resamp(:,i) = resamp(:,i) - mean(resamp(:,i));
    %tx_magnitudes = tx_magnitudes - mean(tx_magnitudes);
%     figure(1); hold on;
%     plot(resamp(:,i));
%     plot(tx_magnitudes);
%     hold off;
    
    
    
    %resamp(:,i) = abs(resamp(:,i)) ./ abs((tx_magnitudes'));
end


end