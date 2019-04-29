function [corr_ConsistencyVsReaction] = consolidation_calculateSimilarity_retrieval(subjectNames, subjectsToRun,freqsOfInterest)

%% tweaks
% - change points from one group to the other depending on their numbe
%   repeat, as opposed to their session #


% identical: 50, 51 (half-half), 52, 53, 55 (5 LTC, 4 PHC, 3 HC, 3 F) -
%   53 has faulty performance records in events.mat
% invalid: 49, 58 (1,7)

%% fixed
evTypeStrings = {'RESPONSE'};
locationStrings = {'Hippocampus', 'Entorhinal', 'Lateral temporal', 'Frontal', 'Parahippocampal'};
freqs = logspace(log10(3),log10(200),70);
studyTypeToRun = 1;
typeOfData = 1; %1=each event zscored to baseline, 2=1-3 s, z-score across sessions

%% changeable parameters
locationToRun = 3:5;
pairwiseOrEnds = 'pairwise'; % 'pairwise' or 'ends'

corr_ConsistencyVsReaction = cell(5,6,5); % 5 regions, 6 subjects, corr and # samples

%for freqs=1:70
for n = subjectsToRun  % 2,4,6 = identical-session subjects, with performance marked in events.mat
    
    load(sprintf('/Volumes/Shares/FRNU/dataWorking/consolidationProject/intermediates/processedPower_allFrequencies_responses/%s/allPowData.mat',subjectNames{n}))
    
    
    % lets build up something that describes which unique word pairs are on
    % which events, across sessions, and the evolution of performance
    % across them
    sessions = unique([evsOfInterest.session]);
    practiceEvents = find([evsOfInterest.list] == 0);
    
    for studyType=studyTypeToRun
        if studyType==1
            evs_grab = find(strcmp({evsOfInterest.type},'RESPONSE'));
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
            useThesePresentations = zeros(156,size(power,2));
            findThesePresentations = NaN(156,size(power,2));
            
            reactionTimes = NaN(156,size(power,2));
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
                        findThesePresentations(word,sess) = match(1);
                    else
                        % lets make sure that cue direction wasnt flipped
                        wordStrings_sess = {evs_session.expecting_word};
                        match = find(strcmp(wordStrings_sess,wordStrings(word)));
                         if size(match,2) ~= 0
                        word_session_performance(1,sess) = evs_session(match(1)).correct;
                        reactionTimes(word,sess) = evs_session(match(1)).RT;
                        findThesePresentations(word,sess) = match(1);
                         else
                             continue
                         end
                    end
                end
                if length(find(word_session_performance==1)) == length(find(~isnan(word_session_performance)))
                    useThesePresentations(word,find(~isnan(word_session_performance))) = 1;
                elseif length(find(word_session_performance==1)) <=1
                    continue
                else
                    % find longest sequence of consecutive corrects
                    f = find(diff([false,word_session_performance==1,false])~=0);
                    [m,ix] = max(f(2:2:end)-f(1:2:end-1));
                    startIndex = f(2*ix-1);
                    if m>1
                        useThesePresentations(word,startIndex:startIndex+m-1) = 1;
                    end
                end
            end
            
            
                        
            % lets do the similarity calculation in a data-driven way
            
            
            if strcmpi(pairwiseOrEnds,'ends')
               
                fprintf('not set up')
                keyboard
                
            elseif strcmpi(pairwiseOrEnds,'pairwise')
                
                % 156 max word pairs, 3 window sizes, # sessions
                consistencyByTime = NaN(156,3,size(power,2)); % 3 window sizes. 100, 250, and 500
                
                if strcmpi(subjectNames{n},'NIH049')
                    doTheseSess = 2:size(power,2)-1;
                elseif strcmpi(subjectNames{n},'NIH051')
                    doTheseSess = 1:size(power,2)-2;
                else
                    doTheseSess = 1:size(power,2)-1;
                end
                    
                % lets do pairwise sessions here
                for sess=1:doTheseSess
                    for i=1:3
                        windowSizes = [5 10 20];
                        power_initial = squeeze(power{4,sess}(:,:,end-windowSizes(i):end,freqsOfInterest));
                        power_initial = squeeze(nanmean(power_initial,4)); %avg across time
                        power_initial = squeeze(nanmean(power_initial,3)); %avg across freqs
                        power_second = squeeze(power{4,sess+1}(:,:,end-windowSizes(i):end,freqsOfInterest));
                        power_second = squeeze(nanmean(power_second,3)); %avg across time
                        power_second = squeeze(nanmean(power_second,3)); %avg across freqs
                        minNumWords = min([size(power_initial,1) size(power_second,1)]);
                        for word=1:minNumWords
                            % only clculate pairwise similarity if the sessions actually
                            % contained both words
                            if (~isnan(findThesePresentations(word,sess))) && (~isnan(findThesePresentations(word,sess+1)))
                            initial_pattern = power_initial(findThesePresentations(word,sess),:);
                            final_pattern = power_second(findThesePresentations(word,sess+1),:);
                            consistencyByTime(word,i,sess) = (1-pdist2(initial_pattern,final_pattern,'cosine'));
                            end
                        end
                    end
                end
                % consistencyByTime = squeeze(nanmean(consistencyByTime,3));
                
            end
            
           
            
            gcf = figure('Units', 'Normalized', 'OuterPosition', [0 0 0.8 0.8]);
            clf
            
            for window=3
                dataKept = squeeze(consistencyByTime(:,window,:));
                % [xfind, yfind] = find(useThesePresentations==1);
                zidx = useThesePresentations==0;
                dataKept(zidx) = NaN;
                reactionTimes(zidx) = NaN;
                %[xRemove, yRemove] = find(isnan(dataKept));
                %dataKept(xRemove, yRemove) = [];
                
                dataKept = nanmean(dataKept,2);
                
                reactionTimeChange = NaN(1,length(dataKept));
                for word=1:156
                    temp_reactionTimes = reactionTimes(word,:);
                    temp_reactionTimes = temp_reactionTimes(find(~isnan(temp_reactionTimes)));
                    meanReactionDiff = nanmean(diff(temp_reactionTimes));
                    reactionTimeChange(1,word) = meanReactionDiff;
                end
                
                subplot(1,1,1)
                
                scatter(dataKept',reactionTimeChange,'o');
                lsline
                
                
                [a,pval] = corr(dataKept,reactionTimeChange','rows','complete','type','Spearman');
                
                xlabel('mean pairwise consistency')
                ylabel('mean reaction time change')
                
                windows = {'250 ms', '500 s', '1s'};
                title(sprintf('window size %s, corr %0.2f',windows{window},a));
                
                set(gca, 'FontName', 'Helvetica');
                set(gca,'FontSize',18)
                hold on
                
                
                % save option
                %             writeFigsHere = '/Volumes/FRNU_SVN/consolidationProject/intermediates/figures/visualizePower/';
                %             if ~exist(sprintf('%s%s',writeFigsHere,subjectNames{n})); mkdir(sprintf('%s%s',writeFigsHere,subjectNames{n})); end
                %             fig2pngSimple(gcf,sprintf('%s%s/%s',writeFigsHere,subjectNames{n},sprintf('%s_studyType%d_pca.png',locationStrings{location},studyType)));
                %
            end
            st = suptitle(sprintf('subj %s, location %s',subjectNames{n},locationStrings{location}));
            set(st,'FontSize',18)
            
            close all
            
            % 3 electrode minimum, at least 5 potentially consolidated words
            if (length(find(~isnan(initial_pattern))) > 3) && (length(find(~isnan(dataKept))) > 5)
            corr_ConsistencyVsReaction{location,n,1} = a;
            corr_ConsistencyVsReaction{location,n,2} = length(find(~isnan(dataKept)));
            corr_ConsistencyVsReaction{location,n,3} = pval;
            corr_ConsistencyVsReaction{location,n,4} = dataKept;
            corr_ConsistencyVsReaction{location,n,5} = reactionTimeChange';
            end
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



