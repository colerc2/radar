% Read in ambient data (tx off). Split into a subfunction so that the user 
% (me) can be lazy when calling loads. Relatively robust error handling.

function ambdat = ambientRead(filename)

if nargin == 0          % Default ambient noise. Most likely case?
    filename = 'default.txt';
elseif nargin > 1
    error('Too many inputs to ambientRead.m');
elseif isempty(strfind(filename,'.txt'))    % If noise is determined to be environment-driven, feed it other data here.
    filename = [filename '.txt'];
end


if exist(['data\ambient\' filename],'file') == 2
    % Load the data if it exists.
    ambdat = load(['data\ambient\' filename]);
else
    error(['Data file (' filename ') does not exist. Fix plz']);
end

ambdat = ambdat(:,1:5001);      % MAKE IT 5001 LENGTH TO MATCH RX DATA

end