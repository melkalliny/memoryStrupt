function cell2csv(cellArray,saveDir)
%- Save a cell array (with both strings and numbers) into csv format
%- by AJ, 2017ish?

assert(~isempty(cellArray), 'Trying to write an empty cell array to CSV!\n');

fid = fopen(saveDir, 'w+') ;

assert(fid > 0, 'Do you have the file open? Cant open %s', saveDir);

for row = 1:size(cellArray,1)
    toWrite = '';
    for col = 1:size(cellArray,2)
        thisVal = cellArray{row,col};
        if isnumeric(thisVal), thisVal = num2str(thisVal); end
        if isempty(thisVal),   thisVal = '-';              end
        toWrite = strcat(toWrite,thisVal);
        if col < size(cellArray,2), toWrite = strcat(toWrite,','); end
    end
    fprintf(fid, '%s\n', toWrite);
end

fclose(fid);




