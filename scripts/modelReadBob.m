% Take in a scenario and range (determines filename) and the number of training pulses to use (M). 

function modeldat = modelReadBob(scenario, range, M)

if nargin ~= 2 && nargin ~= 3
    error('Invalid input to modelRead.m. Revise plz.'); 
end

ifile = [scenario '.txt'];
qfile = [range '.txt'];
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

% No M was input, use all files
if nargin ~= 3
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



% Combine I and Q
modeldat = sqrt(idat.^2 + qdat.^2);

% Zero-mean each signal... NOT SURE IF NEED
for i = 1:M
    modeldat(i,:) = modeldat(i,:) - mean(modeldat(i,:));
end

% MAKE IT 5001 samples long. NOT SURE IF NEED OR WANT
if datlen < 5001
    filler = zeros(M,5001-datlen);
    modeldat = [modeldat filler];
elseif datlen > 5001
    modeldat = modeldat(:,1:5001);
end


end
