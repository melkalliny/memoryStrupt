function [allTables] = consolidation_extractSimilarity_LowAndHigh(subjectNames, subjectsToRun, params)

%% tweaks
% - change points from one group to the other depending on their numbe
%   repeat, as opposed to their session #


% identical: 50, 51 (half-half), 52, 53, 55 (5 LTC, 4 PHC, 3 HC, 3 F) -
%   53 has faulty performance records in events.mat
% invalid: 49, 58 (1,7)

%% fixed
evTypeStrings = {'PROBE_START'};
locationStrings = {'MTL', 'ATL', 'PTL', 'Ventral Frontal', 'Dorsal Frontal', 'PPC'};
freqs = logspace(log10(3),log10(200),70);
studyTypeToRun = 1;
typeOfData = 1; %1=each event zscored to baseline, 2=1-3 s, z-score across sessions

%% changeable parameters
stringEvent = params(1,1);
locationToRun = params(1,2);
freqToDo = params(1,3);
windowToDo = params(1,4);
if stringEvent==2; windowToDo = 90-(20*windowToDo):90; %90 = 3s (response). go x seconds before
elseif stringEvent==1; windowToDo = 30:30+(20*windowToDo);
end

if stringEvent==2; stringEvent = 'RESPONSE';
elseif stringEvent==1; stringEvent = 'PROBE_START';
end

numShuffles = 1;
pairwiseOrEnds = 'pairwise'; % 'pairwise' or 'ends'
numWindows = length(windowToDo);

allTables = cell(length(subjectsToRun),numWindows);



%for freqs=1:70
for n = subjectsToRun  % 2,4,6 = identical-session subjects, with performance marked in events.mat
    
    load(sprintf('/Volumes/FRNU_SVN/consolidation_SVN/consolidationProject/intermediates/processedPower_LowAndHigh_monopolar/%s/allPowData.mat',subjectNames{n}))
    
    
    % lets build up something that describes which unique word pairs are on
    % which events, across sessions, and the evolution of performance
    % across them
    sessions = unique([evsOfInterest.session]);
    practiceEvents = find([evsOfInterest.list] == 0);
    
    evs_grab = find(strcmp({evsOfInterest.type},stringEvent));
    
    
    for location=locationToRun
        
        
        evs_index = cell(1,length(sessions));
        
        % initialize some variables
        power = cell(4,length(sessions));
        behavior = cell(1,length(sessions));
        power_raw = cell(3,length(sessions));
        
        for i=1:length(sessions)
            evs_index{1,i} = find([evsOfInterest.session] == sessions(i));
            
            %                pull_evs = evsOfInterest(intersect(evs_grab,find([evsOfInterest.session] == sessions(i))));
            %                 [a,b] = unique({pull_evs.probe_word});
            %                 pull_evs = pull_evs(b);
            pull_evs_index = intersect(evs_grab,find([evsOfInterest.session] == sessions(i)));
            temp_data = processed_power{typeOfData,location}(pull_evs_index,:,:,:); % first dimension of processed_power has different measures
            behavior{1,i} = [evsOfInterest(pull_evs_index).correct];
            %temp_data = squeeze(nanmean(temp_data,2));
            
            power{4,i} = temp_data;
            power_raw{3,i} = [evsOfInterest(pull_evs_index).correct];
            
        end
        
        % matrix of 156 words x sessions with 1s and 0s
        useThesePresentationsForRT = zeros(156,size(power,2));
        findThesePresentations = NaN(156,size(power,2));
        
        reactionTimes = NaN(156,size(power,2));
        accuracy = NaN(156,size(power,2));
        evsResponse = evsOfInterest(find(strcmp({evsOfInterest.type},'RESPONSE')));
        
        
        wordStrings = {};
        for sess=1:length(sessions)
            evs_session = evsResponse(find([evsResponse.session] == sessions(sess)));
            wordStrings = cat(2,wordStrings,unique({evs_session.probe_word}));
        end
        % wordStrings(find(isempty(wordStrings))) = [];
        wordStrings = unique(wordStrings);
        
        % ok, we cant assume that the word pairs are shown in the same
        % order.. lets match up words and assign data that way
        for word=1:156
            word_session_performance = NaN(1,length(sessions));
            for sess=1:length(sessions)
                evs_session = evsResponse(find([evsResponse.session] == sessions(sess)));
                wordStrings_sess = {evs_session.probe_word};
                match = find(strcmp(wordStrings_sess,wordStrings(word)));
                if size(match,2) ~= 0
                    word_session_performance(1,sess) = evs_session(match(1)).correct;
                    reactionTimes(word,sess) = evs_session(match(1)).RT;
                    accuracy(word,sess) = evs_session(match(1)).correct;
                    findThesePresentations(word,sess) = match(1);
                else
                    % lets make sure that cue direction wasnt flipped
                    wordStrings_sess = {evs_session.expecting_word};
                    match = find(strcmp(wordStrings_sess,wordStrings(word)));
                    if size(match,2) ~= 0
                        word_session_performance(1,sess) = evs_session(match(1)).correct;
                        reactionTimes(word,sess) = evs_session(match(1)).RT;
                        accuracy(word,sess) = evs_session(match(1)).correct;
                        findThesePresentations(word,sess) = match(1);
                    else
                        continue
                    end
                end
            end
            if length(find(word_session_performance==1)) == length(find(~isnan(word_session_performance)))
                useThesePresentationsForRT(word,find(~isnan(word_session_performance))) = 1;
            elseif length(find(word_session_performance==1)) <=1
                continue
            else
                % find longest sequence of consecutive corrects
                f = find(diff([false,word_session_performance==1,false])~=0);
                [m,ix] = max(f(2:2:end)-f(1:2:end-1));
                startIndex = f(2*ix-1);
                if m>1
                    useThesePresentationsForRT(word,startIndex:startIndex+m-1) = 1;
                end
            end
        end
        
        % lets do the similarity calculation in a data-driven way
        
        
        if strcmpi(pairwiseOrEnds,'ends')
            
            fprintf('not set up')
            keyboard
            
        elseif strcmpi(pairwiseOrEnds,'pairwise')
            
            
            if strcmpi(subjectNames{n},'NIH049')
                doTheseSess = 2:size(power,2);
            elseif strcmpi(subjectNames{n},'NIH051')
                doTheseSess = 1:size(power,2)-1;
            elseif strcmpi(subjectNames{n},'NIH043')
                doTheseSess = 4:size(power,2);
            else
                doTheseSess = 1:size(power,2);
            end
            
            % 156 max word pairs, 3 window sizes, # sessions
            consistencyByTime = NaN(156,91,2,size(power,2)-1); % 3 window sizes. 100, 250, and 500
            consistencyByTime_patterns = cell(156,91,2,size(doTheseSess,2));
            
            % lets do pairwise sessions here
            for sess=doTheseSess
                windowsToGrab = 15:105;
                for i=1:91 
                    freqBands = [[1:10];[11:20]];
                    for freqs=freqToDo
                        if size(power{4,sess},2)==0; continue; end
                        power_initial = squeeze(power{4,sess}(:,:,windowsToGrab(i):windowsToGrab(i)+10,freqBands(freqs,:)));
                        power_initial = squeeze(nanmean(power_initial,3)); %avg across time
                        power_initial = reshape(power_initial,size(power_initial,1),size(power_initial,2)*size(power_initial,3));
                        %                             power_second = squeeze(power{4,sess+1}(:,:,windowsToGrab(i):windowsToGrab(i)+10,freqBands(freqs,:)));
                        %                             power_second = squeeze(nanmean(power_second,3)); %avg across time
                        %                             power_second = reshape(power_second,size(power_second,1),size(power_second,2)*size(power_second,3));
                        %
                        for word=7:156 % start at 7 in order to exclude data from practice lists
                            if (~isnan(findThesePresentations(word,sess)))
                                initial_pattern = power_initial(findThesePresentations(word,sess),:);
                                tempRemove = find(isnan(initial_pattern)); initial_pattern(tempRemove) = [];
                                consistencyByTime_patterns{word,i,freqs,sess} = initial_pattern;
                            end
                        end

                    end
                end
            end
            % consistencyByTime = squeeze(nanmean(consistencyByTime,3));
            
        end
        
        
        temp = find(~cellfun(@isempty,consistencyByTime_patterns));
        if length(temp) == 0 %#ok<ISMT>
            continue
        end
        
        for freqs=freqToDo
            
            for window=1:numWindows
                
                dataSubject = {};
                
                consistencyRaw = squeeze(consistencyByTime_patterns(:,windowToDo(window),freqs,:));
                dataRT = reactionTimes;
                dataAcc = accuracy;
                zidx = useThesePresentationsForRT==0;
                dataRT(zidx) = NaN;
                
                sessionLabels = NaN(size(dataRT,1),size(dataRT,2));
                for i=1:size(sessionLabels,2)
                    sessionLabels(:,i) = i;
                end
                
                dataRT = reshape(dataRT,1,size(dataRT,1)*size(dataRT,2));
                dataAcc = reshape(dataAcc,1,size(dataAcc,1)*size(dataAcc,2));
                sessionLabels = reshape(sessionLabels,1,size(sessionLabels,1)*size(sessionLabels,2));
                
                consistency_scores = NaN(1,size(dataRT,2));
                for word=1:size(consistencyRaw,1)
                    temp = find(~cellfun(@isempty,consistencyRaw(word,:)));
                    if length(temp)>1
                        similarity_all = NaN(1,length(temp)-1);
                        for sess=1:length(doTheseSess)-1
                            pattern_first = consistencyRaw{word,doTheseSess(sess)};
                            pattern_second = consistencyRaw{word,doTheseSess(sess)+1};
                            if (length(pattern_first)==length(pattern_second)) && (~isempty(pattern_first))
                                similarity_all(1,sess) = (1-pdist2(pattern_first,pattern_second,'cosine'))+1;
                            end
                        end
                        
                        
                        nComparisons = size(similarity_all,1);
                        similarity_total = nanmean(similarity_all);
                        %consistency_scores(1,word:156:end) = (nComparisons*similarity_total) / (1+(nComparisons-1)*similarity_total);
                        consistency_scores(1,word:156:end) = similarity_total;
                    end
                end
                
                
                % lets go through each word, and calculate a consistency score
                
                %consistencyRaw = reshape(consistencyRaw,1,size(consistencyRaw,1)*size(consistencyRaw,2));
                %                     dataSubject.wordString = repmat(wordStrings_sess,1,length(find(unique(sessionLabels))));
                %                     dataSubject.subject = repmat({subjectNames{n}},1,size(dataRT,2));
                %                     dataSubject.reactionTime = num2cell(dataRT);
                %                     dataSubject.accuracy = num2cell(dataAcc);
                %                     dataSubject.session = num2cell(sessionLabels);
                %                     dataSubject.window = num2cell(repmat(window,1,size(dataRT,2)));
                
                allWords = cell(1,156);
                allWords(1:length(wordStrings_sess)) = wordStrings_sess;
                
                dataSubject(1,:) = repmat(allWords,1,length(find(unique(sessionLabels))));
                dataSubject(2,:) = repmat({subjectNames{n}},1,size(dataRT,2));
                dataSubject(3,:) = num2cell(dataRT);
                dataSubject(4,:) = num2cell(dataAcc);
                dataSubject(5,:) = num2cell(sessionLabels);
                dataSubject(6,:) = num2cell(repmat(window,1,size(dataRT,2)));
                dataSubject(7,:) = num2cell(consistency_scores);
                dataSubject(8,:) = num2cell(repmat(1:156,1,length(find(unique(sessionLabels)))));
                dataSubject(9,:) = num2cell(repmat({doTheseSess},1,size(dataRT,2)));
                
                dataSubject = cell2table(dataSubject','VariableNames', {'WordPairs', 'Subject', 'ReactionTimes', 'Accuracy', 'Session', 'Window', 'Consistency', 'TrialNum', 'SessionsUsed'})   ;
                
                allTables{n,window} = dataSubject;
                
            end
            

            
        end
        
        
        try
            fprintf('\n done with subject %s, location %s\n',subjectNames{1,n},locationStrings{1,location})
        catch
            continue
        end
        
    end
    

    
end



end



% close all
% freqs = logspace(log10(3),log10(200),70);
% figure
% clf
% freqCorrelation(find(isnan(freqCorrelation))) = [];
% plot(1:length(freqCorrelation),freqCorrelation)
% xticks([1 10 20 30 40 50 60 70])
% xticklabels(round(freqs([1 10 20 30 40 50 60 70])))
% ylabel('corr, mean pairwise consistency vs mean pairwise reaction time change')
% xlabel('frequencies')
% st = suptitle(sprintf('subj %s, location %s',subjectNames{n},locationStrings{location}));
% set(st,'FontSize',18)
% set(gca, 'FontName', 'Helvetica');
% set(gca,'FontSize',18)
% keyboard



