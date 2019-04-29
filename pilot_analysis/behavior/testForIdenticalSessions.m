

numSubjects = 49;
identicalCount = cell(length(numSubjects),3);

for subjNum=1:length(numSubjects)
    try
        load(sprintf('/Volumes/Shares/FRNU/data/eeg/NIH0%.2d/behavioral/palRam/events.mat',numSubjects(subjNum)))
    catch
        fprintf('\nno PAL for subj %d\n',numSubjects(subjNum))
        continue
    end
    
    [yesOrNo, numSessions]  = testSubject(events);
    identicalCount{subjNum,1} = yesOrNo;
    identicalCount{subjNum,2} = numSessions;
    identicalCount{subjNum,3} = numSubjects(subjNum);
    fprintf('\ndone, %d\n',numSubjects(subjNum))
end





keyboard








function [yesOrNo, numSessions] = testSubject(events)
% takes an events.mat from PAL and tests to see if
% subsequent sessions were identica
differentSessions = unique([events.session]);

numSessions = NaN; 

% if one session, no repeats
if length(differentSessions) <= 1
    yesOrNo = 0;
    return
end

% populate encoding events for different sessions
encEvents = cell(1,length(differentSessions));
numListsCompleted = NaN(1,length(differentSessions));
for i=1:length(differentSessions)
    eventsSess = events(find([events.session] == differentSessions(i)));
    numListsCompleted(1,i) = length(unique([eventsSess.list]));
    encEvents{1,i} = eventsSess(strcmp({eventsSess.type},'STUDY_PAIR'));
end

% see if lists are identical
wordsList = cell(min(numListsCompleted),length(differentSessions),6,2);
for i=1
    % see if all words are conserved
    for j=1:length(differentSessions)
        % first, get words from this list, then assign
        temp = {encEvents{1,j}(find([encEvents{1,j}.list]==i)).study_1};
        wordsList(i,j,1:length(temp),1) = temp;
        temp = {encEvents{1,j}(find([encEvents{1,j}.list]==i)).study_2};
        wordsList(i,j,1:length(temp),2) = temp;
    end
    
    % if num of unique words in first list, across all sessions, is more than
    % 6, then these were not identical sessions
    temp = squeeze(wordsList(i,:,:,:));
    allWordsList = reshape(temp,1,size(temp,1)*size(temp,2)*size(temp,3));
    allWordsList = allWordsList(find(~cellfun(@isempty,allWordsList)));
    uniqueWordsThisList = length(unique(allWordsList));
    if uniqueWordsThisList > 6
        yesOrNo = 0;
    else
        yesOrNo = 1;
    end
    
    
end

numSessions = length(differentSessions);
end
