ind = randi(length(gclID),10,1);
name = [];
for a = 1:length(ind)
    id = ind(a);
    name = [name,gclID{id}];
end


%spkSubs2 = cellfun(@(x) round(x.*fs), sortedData(goods,2),...
%'UniformOutput', false);
fileID = fopen([dataDir, name, '.txt'], 'w');
formatSpec = '%f \n';
for i = 1:length(ind)
    row = ind(i);
    fprintf(fileID, formatSpec, spontSpks{row,:});
    fprintf(fileID, ';\n');
end
fclose(fileID);