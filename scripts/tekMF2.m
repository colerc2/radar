% Read in two vectors:
%       1. A data matrix that was read from a tek scope O/P file. Each ROW
%       is a unique pulse collection
%       2. A data vector containing the raw transmitted waveform
%
% Output the correlation of these two, after resampling the transmit
% waveform to match the recevied data



function corrdat = tekMF2(rxdat, refdat, fs_scope, fs_tx, Nsub)

[numfiles, rxlen] = size(rxdat);

timescale = fs_scope/fs_tx;         % Resample multiplier
txlen = 2*Nsub + 1;                 % Number of actual samples in transmit waveform

% Resample the transmitted signal to match the 5Gsample/sec of the scope
refdat = refdat(1:txlen);
refdat = resample(refdat,timescale,1);

% Remove DC offset
refdat = refdat - mean(refdat);

% THIS MF/CORRELATION IS USING TIME DOMAIN CORRELATION (XCORR), THUS IT IS
% STUPIDLY SLOW. IF AT ANY POINT YOU CARE ABOUT SPEED, REMOVE THIS CRAP
corrdat = zeros(size(rxdat));
for i = 1:numfiles
    temp = abs(xcorr(rxdat(i,:),refdat));
    corrdat(i,:) = temp(length(rxdat):length(temp));
end

