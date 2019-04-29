function cellArray = csv2cell(csfFileDir)
%Take a csv file with combined numbers and strings, and output a cell array

cellArray = {};
moreToRead = true;
fid    = fopen(csfFileDir,'r');
if fid == -1
    error('Cannot open file: %s', csfFileDir);
end

while moreToRead
    thisLine            = fgetl(fid);
    if ~ischar(thisLine); moreToRead = false; end
    
    if moreToRead
        thisLineCell = {}; thisLineOrig = thisLine;
        % Read this line separated by commas
        isDone = false;
        while ~isDone
            commas_i = strfind(thisLine,',');
            if isempty(commas_i), isDone = true; thisLineCell = cat(2,thisLineCell,thisLine);
            else
                toAdd = thisLine(1:commas_i(1)-1);
                % If the value is empty, add a '-'
                if isempty(toAdd)
                    toAdd = '-';
                end
                thisLineCell = cat(2,thisLineCell,toAdd);
                thisLine = thisLine(commas_i(1)+1:end);
            end
        end
        % For this row, if the last few columns were empty, need to fill it
        % in. Or else we wouldn't be able to concatenate
        if length(thisLineCell) < size(cellArray,2)
            addHowMany = size(cellArray,2)-length(thisLineCell);
            for j = 1:addHowMany, thisLineCell = cat(2,thisLineCell,'-'); end
            
        elseif ~isempty(cellArray) && (length(thisLineCell) > size(cellArray,2))
            % If there's a comma in any of the cells, this code won't work
            % (because it's comma-separated).
            fprintf('ERROR: Check if there''s a comma anywhere in the cells:\n');
            fprintf('%s\n',thisLineOrig);
            keyboard;
        end
        cellArray = cat(1,cellArray,thisLineCell);
    end
    
end

fclose(fid);


end

