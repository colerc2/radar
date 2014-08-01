function numpos = calcNumPos(scenario)

file = ['data\lateral\' scenario '\i_data_pos1.txt'];

if exist(file,'file') ~= 2
    error('Invalid scenario/data. Revise.');
else
    numpos = 1;
end

while(exist(file,'file') == 2)
    numpos = numpos + 1;
    file =  ['data\lateral\' scenario '\i_data_pos' num2str(numpos) '.txt'];
end

numpos = numpos - 1;        % lazy coding, but it works.

end

