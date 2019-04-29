function [corr_ConsistencyVsBehavior] = consolidation_calculateSimilarity_LowAndHigh_biowulf(subjectNum, subjectsToRun, stringEvent, windowType)

%% tweaks
% - change points from one group to the other depending on their numbe
%   repeat, as opposed to their session #

subjectNames = {'NIH043', 'NIH050', 'NIH051', 'NIH052', 'NIH053', 'NIH055'};
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
locationToRun = 1:5;
numShuffles = 100;
pairwiseOrEnds = 'pairwise'; % 'pairwise' or 'ends'

corr_ConsistencyVsBehavior = cell(6,6,1); % 6 regions, 6 subjects, raw data

%for freqs=1:70
subjectsToRun = str2num(subjectsToRun); %#ok<*ST2NM>
windowType = str2num(windowType);

i_full=[111;51;21];


% baseline = {11:20,6:10,4:5};
% adjust = {21:111,11:51,6:21};


for n = subjectsToRun  % 2,4,6 = identical-session subjects, with performance marked in events.mat
    
    fprintf('/data/elkallinymm/consolidationProject/inputs/processedPower_LowAndHigh_monopolar_windows/processedPower_LowAndHigh_monopolar_window%d/%s/allPowData.mat',windowType,subjectNames{n})
    load(sprintf('/data/elkallinymm/consolidationProject/inputs/processedPower_LowAndHigh_monopolar_windows/processedPower_LowAndHigh_monopolar_window%d/%s/allPowData.mat',windowType,subjectNames{n}))
    
    
    % lets build up something that describes which unique word pairs are on
    % which events, across sessions, and the evolution of performance
    % across them
    sessions = unique([evsOfInterest.session]);
    practiceEvents = find([evsOfInterest.list] == 0);
    
    for studyType=studyTypeToRun
        if studyType==1
            evs_grab = find(strcmp({evsOfInterest.type},stringEvent));
        else
            fprintf('\nneeds to be 1\n')
            keyboard
        end

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
                consistencyByTime = NaN(156,i_full(windowType),2,size(power,2)-1); % 3 window sizes. 100, 250, and 500
                consistencyByTime_patterns = cell(156,i_full(windowType),2,size(doTheseSess,2));
                
                % lets do pairwise sessions here
                for sess=doTheseSess
                    windowsToGrab = 1:i_full(windowType);
                    for i=1:i_full(windowType)
                        freqBands = [[1:10];[11:20]];
                        for freqs=1:2
                            if size(power{4,sess},2)==0; continue; end
                            power_initial = squeeze(power{4,sess}(:,:,windowsToGrab(i),freqBands(freqs,:)));
                            %power_initial = squeeze(nanmean(power_initial,3)); %avg across time
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
                            
                            
                            %                             minNumWords = min([size(power_initial,1) size(power_second,1)]);
                            %                             for word=1:156
                            %                                 % only clculate pairwise similarity if the sessions actually
                            %                                 % contained both words
                            %                                 if (~isnan(findThesePresentations(word,sess))) && (~isnan(findThesePresentations(word,sess+1)))
                            %                                     initial_pattern = power_initial(findThesePresentations(word,sess),:);
                            %                                     final_pattern = power_second(findThesePresentations(word,sess+1),:);
                            %                                     consistencyByTime_patterns{word,i,freqs,sess} = (1-pdist2(initial_pattern,final_pattern,'cosine'));
                            %                                 end
                            %                             end
                        end
                    end
                end
                % consistencyByTime = squeeze(nanmean(consistencyByTime,3));
                
            end
            
            
            dataSubject = cell(1,6);
            
            corr_RT= NaN(2,i_full(windowType),3);
            corr_Content = NaN(2,i_full(windowType),3);
            corr_Context = NaN(2,i_full(windowType),3);
            corr_ContentCorrect = NaN(2,i_full(windowType),3);
            corr_ContentIncorrect = NaN(2,i_full(windowType),3);
            
            for freqs=1:2
                
                
                for window=1:i_full(windowType)
                    consistencyRaw = squeeze(consistencyByTime_patterns(:,window,freqs,:));
                    idealizedItemLevel = zeros(size(consistencyRaw,1),size(consistencyRaw,2));
                    for i=1:156
                        idealizedItemLevel(i,:) = i;
                    end
                    dataRT = reactionTimes;
                    dataAcc = accuracy;
                    zidx = useThesePresentationsForRT==0;
                    dataRT(zidx) = NaN;
                    zidx = find(accuracy==0);
                    dataAcc(zidx) = NaN;

                    sessionLabels = NaN(size(dataRT,1),size(dataRT,2));
                    for i=1:size(sessionLabels,2)
                        sessionLabels(:,i) = i;
                    end
                    
                    dataAcc = reshape(dataAcc,1,size(dataAcc,1)*size(dataAcc,2));
                    dataRT = reshape(dataRT,1,size(dataRT,1)*size(dataRT,2));
                    consistencyRaw = reshape(consistencyRaw,1,size(consistencyRaw,1)*size(consistencyRaw,2));
                    idealizedItemLevel = reshape(idealizedItemLevel,1,size(idealizedItemLevel,1)*size(idealizedItemLevel,2));
                    sessionLabels = reshape(sessionLabels,1,size(sessionLabels,1)*size(sessionLabels,2));
                    
                    similarity_matrix = NaN(length(consistencyRaw),length(consistencyRaw));
                    similarity_matrix_corrects = NaN(length(consistencyRaw),length(consistencyRaw));
                    similarity_matrix_incorrects = NaN(length(consistencyRaw),length(consistencyRaw));
                    % ^ only populated if diff session, two corrects being compared
                    idealized_content_matrix = NaN(length(consistencyRaw),length(consistencyRaw));
                    % ^ simply has 1 if same item, diff session, or 0 if
                    % diff item, diff session. NaN otherwise
                     RT_matrix = NaN(length(consistencyRaw),length(consistencyRaw));
                     % only fill the RT matrix with same words, diff sessions
                    idealized_context_matrix = NaN(length(consistencyRaw),length(consistencyRaw));
                    acc_matrix = NaN(length(consistencyRaw),length(consistencyRaw));
                    
                    
                
                    for i=1:length(consistencyRaw)
                        for j=1:length(consistencyRaw)
                            if i==j
                                continue % never compare identical patterns
                            end
                            
                            
                            pattern_first = consistencyRaw{1,i};
                            pattern_second = consistencyRaw{1,j};
                            if ~isempty(pattern_first) && ~isempty(pattern_second)
                                similarity_matrix(i,j) = (1-pdist2(pattern_first,pattern_second,'cosine'));
                            end
                            
                            
                            if sessionLabels(1,i)==sessionLabels(1,j)
                                idealized_context_matrix(i,j) = 1;
                            elseif abs(sessionLabels(1,i)-sessionLabels(1,j))==1 % 1 session apart
                                idealized_context_matrix(i,j) = 0;
                            elseif abs(sessionLabels(1,i)-sessionLabels(1,j))>=2 % 2 sessions apart
                                idealized_context_matrix(i,j) = 0;
                            end
                            
                            % if it was the same word
                            if idealizedItemLevel(1,i) == idealizedItemLevel(1,j)
                                
                                % if same word, expect correlation
                                idealized_content_matrix(i,j) = 1;
                                % only take RT if we're going to subtract
                                % RT(earlier)-RT(later) -- (so that higher
                                % value indicates learning)
                                if (sessionLabels(1,i)>sessionLabels(1,j))==-1
                                    RT_matrix(i,j) = dataRT(1,i) - dataRT(1,j);
                                end
                                
                            else
                                
                               idealized_content_matrix(i,j) = 0;
                               
                            end
                            
                            acc_matrix(i,j) = nansum([dataAcc(1,i),dataAcc(1,j)]);
                            if ~isempty(pattern_first) && ~isempty(pattern_second)
                                if (acc_matrix(i,j)==2) && (sessionLabels(1,i)~=sessionLabels(1,j))
                                    similarity_matrix_corrects(i,j) = (1-pdist2(pattern_first,pattern_second,'cosine'));
                                elseif (acc_matrix(i,j)==0) && (sessionLabels(1,i)~=sessionLabels(1,j))
                                    similarity_matrix_incorrects(i,j) = (1-pdist2(pattern_first,pattern_second,'cosine'));
                                end
                            end
                            
                            
                        end
                    end
                
                    
                    
                    ii=ones(size(similarity_matrix));
                    idx=(tril(ii,-1));
                    similarity_matrix(~idx)=nan;
                    
%                     % lets z-score the item level contrast weights
%                     tempmean = nanmean(idealized_content_matrix(:));
%                     tempstd = nanstd(idealized_content_matrix(:));
%                     idealized_content_matrix = (idealized_content_matrix-tempmean) / tempstd;
%                     % and the session level contrast weights
%                     tempmean = nanmean(idealized_context_matrix(:));
%                     tempstd = nanstd(idealized_context_matrix(:));
%                     idealized_context_matrix = (idealized_context_matrix-tempmean) / tempstd;
                    
                    content_for_corrects_matrix = idealized_content_matrix; %to be used later
                    content_for_incorrects_matrix = idealized_content_matrix; %to be used later
                    
                    % lets remove the words for which we have no consistency metric
                    removeThese = find(isnan(similarity_matrix));
                    similarity_matrix(removeThese) = [];
                    idealized_content_matrix(removeThese) = [];
                    RT_matrix(removeThese) = [];
                    idealized_context_matrix(removeThese) = [];
                    
                    try
                    % if it fails, then there's no data in this brain region
                    [rho,~] = corr(similarity_matrix',idealized_content_matrix','rows','complete','Type','Spearman');
                    corr_Content(freqs,window) = rho;
                    % produce shuffled data
                    temp_matrix1 = similarity_matrix; temp_matrix2 = idealized_content_matrix;
                    temp_remove = find(isnan(temp_matrix1));
                    temp_matrix1(temp_remove) = []; temp_matrix2(temp_remove) = [];
                    temp_remove = find(isnan(temp_matrix2));
                    temp_matrix1(temp_remove) = []; temp_matrix2(temp_remove) = [];  
                    rhoShuffles = NaN(1,numShuffles);
                    for shuffle=1:numShuffles
                       temp_matrix1 = temp_matrix1(randperm(length(temp_matrix1)));
                       [rhoShuffles(1,shuffle),~] = corr(temp_matrix1',temp_matrix2','rows','complete','Type','Spearman');
                    end
                    errorStd = nanstd(rhoShuffles);
                    errorMean = nanmean(rhoShuffles);
                    corr_Content(freqs,window,2) = errorMean;
                    corr_Content(freqs,window,3) = errorStd;
                    
                    
                    catch
                        continue
                    end
                    
                    % similarity vs reaction time changes
                    try
                    [rho,~] = corr(similarity_matrix',RT_matrix','rows','complete','Type','Spearman');
                    corr_RT(freqs,window,1) = rho;
                    
                    % produce shuffled data
                    temp_matrix1 = similarity_matrix; temp_matrix2 = RT_matrix;
                    temp_remove = find(isnan(temp_matrix1));
                    temp_matrix1(temp_remove) = []; temp_matrix2(temp_remove) = [];
                    temp_remove = find(isnan(temp_matrix2));
                    temp_matrix1(temp_remove) = []; temp_matrix2(temp_remove) = [];  
                    rhoShuffles = NaN(1,numShuffles);
                    for shuffle=1:numShuffles
                       temp_matrix1 = temp_matrix1(randperm(length(temp_matrix1)));
                       [rhoShuffles(1,shuffle),~] = corr(temp_matrix1',temp_matrix2','rows','complete','Type','Spearman');
                    end
                    errorStd = nanstd(rhoShuffles);
                    errorMean = nanmean(rhoShuffles);
                    corr_RT(freqs,window,2) = errorMean;
                    corr_RT(freqs,window,3) = errorStd;
                    catch
                        continue
                        end
                    
                    
                    % similarity vs session null model
                    [rho,~] = corr(similarity_matrix',idealized_context_matrix','rows','complete','Type','Spearman');
                    corr_Context(freqs,window,1) = rho;
                    
                    % produce shuffled data
                    temp_matrix1 = similarity_matrix; temp_matrix2 = idealized_context_matrix;
                    temp_remove = find(isnan(temp_matrix1));
                    temp_matrix1(temp_remove) = []; temp_matrix2(temp_remove) = [];
                    temp_remove = find(isnan(temp_matrix2));
                    temp_matrix1(temp_remove) = []; temp_matrix2(temp_remove) = [];  
                    rhoShuffles = NaN(1,numShuffles);
                    for shuffle=1:numShuffles
                       temp_matrix1 = temp_matrix1(randperm(length(temp_matrix1)));
                       [rhoShuffles(1,shuffle),~] = corr(temp_matrix1',temp_matrix2','rows','complete','Type','Spearman');
                    end
                    errorStd = nanstd(rhoShuffles);
                    errorMean = nanmean(rhoShuffles);
                    corr_Context(freqs,window,2) = errorMean;
                    corr_Context(freqs,window,3) = errorStd;
                    
                    
                    
                    % similarity vs content-correct 
                    % need to limit the similarity_matrix to corrects
                    
                    removeThese = find(isnan(similarity_matrix_corrects));
                    similarity_matrix_corrects(removeThese) = [];
                    content_for_corrects_matrix(removeThese) = [];
                    
                    
                    [rho,~] = corr(similarity_matrix_corrects',content_for_corrects_matrix','rows','complete','Type','Spearman');
                    corr_ContentCorrect(freqs,window,1) = rho;
                    
                    % produce shuffled data
                    temp_matrix1 = similarity_matrix_corrects; temp_matrix2 = content_for_corrects_matrix;
                    temp_remove = find(isnan(temp_matrix1));
                    temp_matrix1(temp_remove) = []; temp_matrix2(temp_remove) = [];
                    temp_remove = find(isnan(temp_matrix2));
                    temp_matrix1(temp_remove) = []; temp_matrix2(temp_remove) = [];  
                    rhoShuffles = NaN(1,numShuffles);
                    for shuffle=1:numShuffles
                       temp_matrix1 = temp_matrix1(randperm(length(temp_matrix1)));
                       [rhoShuffles(1,shuffle),~] = corr(temp_matrix1',temp_matrix2','rows','complete','Type','Spearman');
                    end
                    errorStd = nanstd(rhoShuffles);
                    errorMean = nanmean(rhoShuffles);
                    corr_ContentCorrect(freqs,window,2) = errorMean;
                    corr_ContentCorrect(freqs,window,3) = errorStd;  
                    
                    [rho,~] = corr(similarity_matrix_corrects',content_for_corrects_matrix','rows','complete','Type','Spearman');
                    corr_ContentCorrect(freqs,window,1) = rho;
                    
                    % produce shuffled data
                    temp_matrix1 = similarity_matrix_corrects; temp_matrix2 = content_for_corrects_matrix;
                    temp_remove = find(isnan(temp_matrix1));
                    temp_matrix1(temp_remove) = []; temp_matrix2(temp_remove) = [];
                    temp_remove = find(isnan(temp_matrix2));
                    temp_matrix1(temp_remove) = []; temp_matrix2(temp_remove) = [];
                    rhoShuffles = NaN(1,numShuffles);
                    for shuffle=1:numShuffles
                        temp_matrix1 = temp_matrix1(randperm(length(temp_matrix1)));
                        [rhoShuffles(1,shuffle),~] = corr(temp_matrix1',temp_matrix2','rows','complete','Type','Spearman');
                    end
                    errorStd = nanstd(rhoShuffles);
                    errorMean = nanmean(rhoShuffles);
                    corr_ContentCorrect(freqs,window,2) = errorMean;
                    corr_ContentCorrect(freqs,window,3) = errorStd;
                    
                    
                    removeThese = find(isnan(similarity_matrix_incorrects));
                    similarity_matrix_incorrects(removeThese) = [];
                    content_for_incorrects_matrix(removeThese) = [];
                    
                    
                    [rho,~] = corr(similarity_matrix_incorrects',content_for_incorrects_matrix','rows','complete','Type','Spearman');
                    corr_ContentIncorrect(freqs,window,1) = rho;
                    
                    % produce shuffled data
                    temp_matrix1 = similarity_matrix_incorrects; temp_matrix2 = content_for_incorrects_matrix;
                    temp_remove = find(isnan(temp_matrix1));
                    temp_matrix1(temp_remove) = []; temp_matrix2(temp_remove) = [];
                    temp_remove = find(isnan(temp_matrix2));
                    temp_matrix1(temp_remove) = []; temp_matrix2(temp_remove) = [];
                    rhoShuffles = NaN(1,numShuffles);
                    for shuffle=1:numShuffles
                        temp_matrix1 = temp_matrix1(randperm(length(temp_matrix1)));
                        [rhoShuffles(1,shuffle),~] = corr(temp_matrix1',temp_matrix2','rows','complete','Type','Spearman');
                    end
                    errorStd = nanstd(rhoShuffles);
                    errorMean = nanmean(rhoShuffles);
                    corr_ContentIncorrect(freqs,window,2) = errorMean;
                    corr_ContentIncorrect(freqs,window,3) = errorStd;
                                       
                    %                 end
                   % fprintf('\n done with window %0.2f\n',window)
                end
            
            
            
            
            
            
            end
            
            dataSubject{1,1} = corr_Context;
            dataSubject{1,2} = corr_Content;
            dataSubject{1,3} = corr_RT;
            dataSubject{1,4} = corr_ContentCorrect;
            dataSubject{1,5} = corr_ContentIncorrect; %this is uninterpretable if stringEvent==1, bc in that case
            % 'response' locked is actually arbitrary
            
            %% remember to filter out incorrects
            reactionTimesOutput = reactionTimes;
            reactionTimesOutput = reactionTimesOutput(find(accuracy==1));
            reactionTimesOutput(find(reactionTimesOutput==-999)) = [];
            dataSubject{1,6} = nanmedian(reactionTimesOutput(:));
            
            corr_ConsistencyVsBehavior{location,n,1} = dataSubject;
            
            
            fprintf('\n done with subject %s, location %s\n',subjectNames{1,n},locationStrings{1,location})
          
            
        end
        
        fprintf('/data/elkallinymm/consolidationProject/outputs/windows/%s/dataSubject_%s_%s.mat',subjectNames{n},stringEvent,num2str(windowType))
        
        save(sprintf('/data/elkallinymm/consolidationProject/outputs/windows/%s/dataSubject_%s_%s.mat',subjectNames{n},stringEvent,num2str(windowType)),'dataSubject','corr_ConsistencyVsBehavior')
        
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



