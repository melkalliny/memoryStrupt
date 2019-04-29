function [pointsOfInterest] = consolidation_visualizePower_allFrequencies()

%% tweaks
% - change points from one group to the other depending on their numbe
%   repeat, as opposed to their session #


subjects = [49 50 51 52 53 55 58];
subjectNames = {'NIH049', 'NIH050', 'NIH051', 'NIH052', 'NIH053', 'NIH055', 'NIH058'};
npats = [1 2 3 4 5 6 7];

% identical: 50, 51 (half-half), 52, 53, 55 (5 LTC, 4 PHC, 3 HC, 3 F) -
%   53 has faulty performance records in events.mat
% invalid: 49, 58 (1,7)

%% fixed
evTypeStrings = {'STUDY', 'TEST'};
locationStrings = {'Hippocampus', 'Entorhinal', 'Lateral temporal', 'Frontal', 'Parahippocampal'};
%% changeable parameters
locationToRun = [1,3,5];
studyTypeToRun = [1,2];
typeOfData = 1; %1=each event zscored to baseline, 2=1-3 s, z-score across sessions
subjectsToRun = [4 6];

pointsOfInterest = cell(2,5,7);

for n = subjectsToRun  % 2,4,6 = identical-session subjects, with performance marked in events.mat
    
    load(sprintf('/Volumes/Shares/FRNU/dataWorking/consolidationProject/intermediates/processedPower_allFrequencies/%s/allPowData.mat',subjectNames{n}))
    
    
    % lets build up something that describes which unique word pairs are on
    % which events, across sessions, and the evolution of performance
    % across them
    sessions = unique([evsOfInterest.session]);
    practiceEvents = find([evsOfInterest.list] == 0);
    
    evs_ret = evsOfInterest(find(strcmp({evsOfInterest.type},'PROBE_START')));
    clear all_word_pairs word_pairs_category
    for session=1:length(sessions)
        evs_sess = evs_ret(find([evs_ret.session] == sessions(session)));
        unique_pairs = unique({evs_sess.probe_word});
        
        if session==1 % if first session
            all_word_pairs = unique_pairs;
            word_pairs_category = NaN(length(unique_pairs),length(sessions));
            
        elseif session>1
            if length(unique_pairs) <= length(all_word_pairs)
                continue
            else
                keyboard
            end
            %check to see if we have any words in this
            % session that arent already included in the all_word_pairs
            % variable. if not, continue
        end
        
        for i=1:length(all_word_pairs)
            temp_event_matches = [];
            % find all events in which this word was used
            matches = find(strcmp({evs_ret.probe_word},unique_pairs(i)));
            temp_event_matches = cat(2,temp_event_matches,matches);
            matches = find(strcmp({evs_ret.expecting_word},unique_pairs(i)));
            temp_event_matches = cat(2,temp_event_matches,matches);
            
            % check that there are 2 or greater, and then calculate which of
            % the categories it is in
            if size(temp_event_matches,2) > 2
                temp_events = evs_ret(temp_event_matches);
                temp_performance = [temp_events.correct];
                diffs = diff(temp_performance);
                
                if length(find(diff([temp_events.session]) < 0) > 0)
                    fprintf('\nwords are not in sequential order, resorting\n')
                    [a,resort] = sort([temp_events.session]);
                    temp_performance = temp_performance(resort);
                    diffs = diff(temp_performance);
                end
                
                word_pairs_category(i,1:size(temp_performance,2)) = temp_performance;
                %                if length(find(diffs<0)) > 0 %#ok<*ISMT> %negative change
                %                    word_pairs_category(i,:) = -1;
                %                continue
                %                end
                %
                %                if length(find(temp_performance==1)) == length(temp_performance)
                %                    word_pairs_category(i) = 0; %perfect performance
                %                    continue
                %                else
                %                    temp_performance
                %                    word_pairs_category(i) = 1; %strictly increasing
                %                    continue
                %                end
                
            else
                continue
                % word_pairs_category(i,1) = NaN; % word pair had no repeats
            end
        end
        
    end
    
    % lets find how much time separated each session
    events = cell(1,length(sessions));
    differenceInDays = NaN(1,length(sessions)-1);
    for i=1:length(sessions)-1
        temp_evs = evsOfInterest(find([evsOfInterest.session] == sessions(i)));
        temp_evs_nextSession = evsOfInterest(find([evsOfInterest.session] == sessions(i+1)));
        msTimeDiff = temp_evs_nextSession(1).mstime - temp_evs(1).mstime;
        differenceInDays(i) = msTimeDiff / 1000 / 60 / 60 / 24;
    end
    allDifferences = '';
    for i=1:length(differenceInDays)
        if i==1
            allDifferences = cat(2,allDifferences,sprintf('%s',num2str(differenceInDays(i))));
        else
            allDifferences = cat(2,allDifferences,sprintf(', %s',num2str(differenceInDays(i))));
        end
    end
    
    
    for studyType=studyTypeToRun
        if studyType==1
            evs_grab = find(strcmp({evsOfInterest.type},'STUDY_PAIR_START'));
            % lets grab retrieval events even during encoding, so we that can visualize performance
            % (which is only accurately marked during retrieval)
            evs_grab_retrieval = find(strcmp({evsOfInterest.type},'PROBE_START'));
        elseif studyType==2
            evs_grab = find(strcmp({evsOfInterest.type},'PROBE_START'));
            evs_grab_retrieval = find(strcmp({evsOfInterest.type},'PROBE_START'));
        end
        
        
        
        
        for location=locationToRun
            
            
            evs_index = cell(1,length(sessions));
            
            % initialize some variables
            power = cell(4,length(sessions));
            behavior = cell(1,length(sessions));
            power_raw = cell(3,length(sessions));
            percentCorrect = NaN(3,length(sessions));
            power_by_category = cell(3,length(sessions));
            
            for i=1:length(sessions)
                evs_index{1,i} = find([evsOfInterest.session] == sessions(i));
                
                pull_evs = evsOfInterest(intersect(evs_grab,find([evsOfInterest.session] == sessions(i))));
                category_word = NaN(length(pull_evs),length(sessions));
                session_words = NaN(1,length(pull_evs));
                for word=1:length(pull_evs)
                    temp_word = pull_evs(word).probe_word;
                    match = find(strcmp(all_word_pairs,temp_word));
                    if size(match,2)==0
                        temp_word = pull_evs(word).expecting_word;
                        match = find(strcmp(all_word_pairs,temp_word));
                    end
                    category_word(word,:) = word_pairs_category(word,:);
                    
                    % now lets use the performance in consecutive word
                    % presentations to determine whether to assign this
                    % word - in this session - into consolidated or
                    % non-consolidated
                    
                    if i==1
                        if (category_word(word,i) == 1) && (category_word(word,i+1) == 1)
                            session_words(word) = 1; % consolidated category
                        elseif (category_word(word,i) == 1) && (category_word(word,i+1) == 0)
                            session_words(word) = -1; % unsure category
                        elseif (category_word(word,i) == 0) && (category_word(word,i+1) == 1)
                            session_words(word) = 1; % consolidated category
                        elseif (category_word(word,i) == 0) && (category_word(word,i+1) == 0)
                            session_words(word) = -1; % not consolidated category
                        end
                        continue
                    else
                        if (category_word(word,i) == 1) && (category_word(word,i-1) == 0)
                            session_words(word) = 1; % consolidated category
                        elseif (category_word(word,i) == 0) && (category_word(word,i-1) == 1)
                            session_words(word) = -1; % not consolidated category
                        elseif (category_word(word,i) == 1) && (category_word(word,i-1) == 1)
                            session_words(word) = 1; % consolidated category
                        elseif (category_word(word,i) == 0) && (category_word(word,i-1) == 0)
                            session_words(word) = -1; % unsure category
                        end
                        
                        continue
                    end
                    
                    
                    
                end
                
                
                
                pull_evs_index = intersect(evs_grab,find([evsOfInterest.session] == sessions(i)));
                temp_data = processed_power{typeOfData,location}(pull_evs_index,:,:,:); % first dimension of processed_power has different measures
                behavior{1,i} = [evsOfInterest(pull_evs_index).correct];
                %temp_data = squeeze(nanmean(temp_data,2));
                power{1,i} = (temp_data(find(session_words==-1),:,:,:));
                power{2,i} = (temp_data(find(session_words==0),:,:,:));
                power{3,i} = (temp_data(find(session_words==1),:,:,:));
                power{4,i} = temp_data;
                power_raw{3,i} = [evsOfInterest(pull_evs_index).correct];
                
                % and let write out the performance for each event, not
                % just the averages, so that we can do some stuff with it
                % later. if the corrects arent properly labeled during
                % encoding then this will not be accurate
                temp = [evsOfInterest(pull_evs_index(find(session_words==-1))).correct]; temp(find(temp==-999)) = [];
                percentCorrect(2,i) = sum(temp) / length(temp);% non-consolidated
                temp = [evsOfInterest(pull_evs_index(find(session_words==1))).correct];  temp(find(temp==-999)) = [];
                percentCorrect(3,i) = sum(temp) / length(temp); % consolidated
                temp = [evsOfInterest(pull_evs_index).correct]; temp(find(temp==-999)) = [];
                percentCorrect(1,i) = sum(temp) / length(temp); % all
                
                
            end
            
            %             if length(find(isnan(power))) >= (size(power,1)*size(power,2))
            %                 fprintf('\nno data, for %s, in %s\n',subjectNames{n},locationStrings{location})
            %                 continue
            %             end
            
            gcf = figure('Units', 'Normalized', 'OuterPosition', [0 0 0.8 0.8]);
            set(gcf,'Color','w')
            
            markersize = 0.3;
            
            
            gcf = figure('Units', 'Normalized', 'OuterPosition', [0 0 0.8 0.8]);
            set(gcf,'Color','w')
            
            tStatMat = cell(1,length(sessions));
            
            all_correct_observs = [];
            all_incorrect_observs = [];
            
            for i=1:length(sessions)
                %                subplot(length(sessions),3,1*i)
                power_plot = power{4,i}(find(behavior{1,i}==1),:,:,:);
                
                all_correct_observs = cat(1,all_correct_observs,squeeze(nanmean(power_plot,2)));
                
                power_plot = squeeze(nanmean(power_plot,1));
                power_plot = squeeze(nanmean(power_plot,1));
                %imagesc(flipud(power_plot'))
                correct_plot = power_plot;
                %                 yticks([10:10:70]); xticks([30 50 70 90 110])
                %                 freqs = logspace(log10(3),log10(200),70);
                %                 yticklabels(fliplr(round(freqs(1:10:70))))
                %                 xticklabels([0 1 2 3 4])
                %                 ylabel('frequencies'); xlabel('time (s)');
                %                 title('correct')
                %                 box off
                %                 colorbar
                %                 set(gca, 'FontName', 'Helvetica');
                %                 set(gca,'FontSize',18)
                
                %                subplot(1,length(sessions),i)
                power_plot = power{4,i}(find(behavior{1,i}==0),:,:,:);
                
                all_incorrect_observs = cat(1,all_incorrect_observs,squeeze(nanmean(power_plot,2)));
                
                power_plot = squeeze(nanmean(power_plot,1));
                power_plot = squeeze(nanmean(power_plot,1));
                %                 imagesc(flipud(power_plot'))
                %                 yticks([10:10:70]); xticks([30 50 70 90 110])
                %                 freqs = logspace(log10(3),log10(200),70);
                %                 yticklabels(fliplr(round(freqs(1:10:70))))
                %                 xticklabels([0 1 2 3 4])
                %                 ylabel('frequencies'); xlabel('time (s)');
                %                 title('incorrect')
                %                 box off
                %                 set(gca, 'FontName', 'Helvetica');
                %                 set(gca,'FontSize',18)
                
                subplot(1,length(sessions),i)
                imagesc(flipud(correct_plot' - power_plot'))
                yticks([10:10:70]); xticks([30 50 70 90 110])
                freqs = logspace(log10(3),log10(200),70);
                yticklabels(fliplr(round(freqs(1:10:70))))
                xticklabels([0 1 2 3 4])
                ylabel('frequencies'); xlabel('time (s)');
                title('corr-incorr')
                box off
                caxis([-0.5 0.5])
                set(gca, 'FontName', 'Helvetica');
                set(gca,'FontSize',18)
                
                
               % pointsOfInterest{
                
                
            end
            % colorbar
            st = suptitle(sprintf('%s, %s, %s',subjectNames{n},locationStrings{1,location},evTypeStrings{1,studyType}));
            set(st,'FontSize',18);
            
            [~,~,~,d] = ttest2(all_correct_observs,all_incorrect_observs);
            tstats = flipud(squeeze(d.tstat));
            pointsOfInterest{studyType,location,n} = tstats;
            % save option
            %             writeFigsHere = '/Volumes/FRNU_SVN/consolidationProject/intermediates/figures/visualizePower/';
            %             if ~exist(sprintf('%s%s',writeFigsHere,subjectNames{n})); mkdir(sprintf('%s%s',writeFigsHere,subjectNames{n})); end
            %             fig2pngSimple(gcf,sprintf('%s%s/%s',writeFigsHere,subjectNames{n},sprintf('%s_studyType%d_pca.png',locationStrings{location},studyType)));
            %
        end
        
        
        
        % now, do the same but for hippocampus/cortex ratio of power
    end
    
end

end