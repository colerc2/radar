% fork of modelRead(). Nothing really different, just different directory
% structure.

function [dat truth] = lateralRead(scenario, posnum)

if nargin ~= 2
    error('Invalid input to modelRead.m. Revise plz.'); 
end

ifile = ['data\lateral\' scenario '\i_data_pos' num2str(posnum) '.txt'];
qfile = ['data\lateral\' scenario '\q_data_pos' num2str(posnum) '.txt'];
truthfile = ['data\lateral\' scenario '\truepos.txt'];

% Load the data from .txt if it exists. Error if it doesn't.
if exist(ifile,'file') == 2 && exist(qfile,'file') == 2
    idat = load(ifile);
    qdat = load(qfile);
    truth = load(truthfile);
else
    error('Improper scenario and/or range given. Check existing models.')
end

% Make sure the I and Q data are the same in shape.
if size(idat) ~= size(qdat)
    error('I and Q data sizes do not match.');
end

[numpulse datlen] = size(idat);


% Combine I and Q
dat = sqrt(idat.^2 + qdat.^2);

% Zero-mean each signal... NOT SURE IF NEED
for i = 1:numpulse
    dat(i,:) = dat(i,:) - mean(dat(i,:));
end

% MAKE IT 5001 samples long. NOT SURE IF NEED OR WANT
if datlen < 5001
    filler = zeros(numpulse,5001-datlen);
    dat = [dat filler];
elseif datlen > 5001
    dat = dat(:,1:5001);
end


end
