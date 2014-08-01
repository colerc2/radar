
function [modeldat, idat, qdat] = model_read_angle(scenario, sub_scenario, angle, M)

if nargin ~= 2 && nargin ~= 3
    error('Invalid input to modelRead.m. Revise plz.'); 
end

ifile = ['data\' scenario '\i_data_' sub_scenario '_' num2str(angle) 'deg.txt'];
qfile = ['data\' scenario '\q_data_' sub_scenario '_' num2str(angle) 'deg.txt'];
%ifile = ['data\models\i_data_' scenario '_' num2str(range) '.txt'];
%qfile = ['data\models\q_data_' scenario '_' num2str(range) '.txt'];
%ifile = ['data\fpm_gage_test\i_data_' scenario '_' num2str(range) '.txt'];
%qfile = ['data\fpm_gage_test\q_data_' scenario '_' num2str(range) '.txt'];
%ifile
%qfile

% Load the data from .txt if it exists. Error if it doesn't.
if exist(ifile,'file') == 2 && exist(qfile,'file') == 2
    idat = load(ifile);
    qdat = load(qfile);
else
    error('Improper scenario and/or range given. Check existing models.')
end

% Make sure the I and Q data are the same in shape.
if size(idat) ~= size(qdat)
    error('I and Q data sizes do not match.');
end

[numpulse datlen] = size(idat);
%display(numpulse);

% No M was input, use all files
if nargin ~= 4
    M = numpulse;
elseif M > numpulse
    error('M was specified as too large.');
end

% take a randperm() selection of pulses, use that as the model
if M < numpulse
    inds = randperm(numpulse);
    inds = sort(inds(1:M));
    
    idat = idat(inds,:);
    qdat = qdat(inds,:);
end

% for i = 1:M
%     idat(i,:) = idat(i,:) - mean(idat(i,:));
%     qdat(i,:) = qdat(i,:) - mean(qdat(i,:));
% end


% Combine I and Q
%modeldat = sqrt(idat.^2 + qdat.^2);
%modeldat = idat;

% Zero-mean each signal... NOT SURE IF NEED
for i = 1:M
    modeldat(i,:) = sqrt(idat(i,:).^2 + qdat(i,:).^2);
    modeldat(i,:) = modeldat(i,:) - mean(modeldat(i,:));
    idat(i,:) = idat(i,:) - mean(idat(i,:));
    qdat(i,:) = qdat(i,:) - mean(qdat(i,:));

end

% MAKE IT 5001 samples long. NOT SURE IF NEED OR WANT
if datlen < 5001
    filler = zeros(M,5001-datlen);
    modeldat = [modeldat filler];
    idat = [idat filler];
    qdat = [qdat filler];
elseif datlen > 5001
    modeldat = modeldat(:,1:5001);
end

end
