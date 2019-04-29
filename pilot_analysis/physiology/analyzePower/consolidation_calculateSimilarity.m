function [pointsOfInterest] = consolidation_calculateSimilarity(pointsOfInterest)

%% tweaks
% - change points from one group to the other depending on their numbe
%   repeat, as opposed to their session #

if ~exist('pointsOfInterest')
    pointsOfInterest = [];
end


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
locationToRun = [3,5];
studyTypeToRun = [2];
typeOfData = 1; %1=each event zscored to baseline, 2=1-3 s, z-score across sessions
subjectsToRun = [4,6];
tStatThresh = 3;
pairwiseOrEnds = 'pairwise'; %'pairwise' or 'ends'

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
                
                
            else
                continue
                % word_pairs_category(i,1) = NaN; % word pair had no repeats
            end
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
                
                
                pull_evs_index = intersect(evs_grab,find([evsOfInterest.session] == sessions(i)));
                temp_data = processed_power{typeOfData,location}(pull_evs_index,:,:,:); % first dimension of processed_power has different measures
                behavior{1,i} = [evsOfInterest(pull_evs_index).correct];
                %temp_data = squeeze(nanmean(temp_data,2));
                
                power{4,i} = temp_data;
                power_raw{3,i} = [evsOfInterest(pull_evs_index).correct];
                
            end
            
            freqs = logspace(log10(3),log10(200),70);
            freqsOfInterest = 60:70;
            
            % lets do the similarity calculation in a data-driven way
            
            % lets add a toggle here for measuring pairwise sessions or
            % first to last session
            
            if strcmpi(pairwiseOrEnds,'ends')
                % lets do first to last session here
                consistencyByTime = NaN(156,length(30:110));
                consistencyROI = NaN(1,156);
                for i=25:105
                    power_initial = squeeze(power{4,1}(:,:,i-5:i+5,freqsOfInterest));
                    power_initial = squeeze(nanmean(power_initial,3)); %avg across time
                    power_initial = squeeze(nanmean(power_initial,3)); %avg across freqs
                    power_final = squeeze(power{4,end}(:,:,i-5:i+5,freqsOfInterest));
                    power_final = squeeze(nanmean(power_final,3)); %avg across time
                    power_final = squeeze(nanmean(power_final,3)); %avg across freqs
                    minNumWords = min([size(power_initial,1) size(power_final,1)]);
                    for word=1:minNumWords
                        initial_pattern = power_initial(word,:);
                        final_pattern = power_final(word,:);
                        consistencyByTime(word,i-24) = (1-pdist2(initial_pattern,final_pattern,'cosine'));
                    end
                end
                % now, lets incorporate correct-vs-incorrect in choosing what
                % to calculate similarity for
                tStats = pointsOfInterest{studyType,location,n};
                [xFind,yFind] = find((tStats) > tStatThresh);
                % eliminating certain frequencies and time points
                %             yFind(find(xFind<30)) = []; yFind(find(xFind>110)) = [];
                %             xFind(find(xFind<30)) = []; xFind(find(xFind>110)) = [];
                %             xFind(find(yFind<60)) = []; yFind(find(yFind<60)) = [];
                power_initial = squeeze(power{4,1}(:,:,xFind,yFind));
                power_final = squeeze(power{4,end}(:,:,xFind,yFind));
                minNumWords = min([size(power_initial,1) size(power_final,1)]);
                for word=1:minNumWords
                    initial_pattern = power_initial(word,:,:,:);
                    initial_pattern = reshape(initial_pattern,1,size(initial_pattern,2)*size(initial_pattern,3)*size(initial_pattern,4));
                    final_pattern = power_final(word,:);
                    final_pattern = reshape(final_pattern,1,size(final_pattern,2)*size(final_pattern,3)*size(final_pattern,4));
                    consistencyROI(1,word) = (1-pdist2(initial_pattern,final_pattern,'cosine'));
                end
                
            else
                
                consistencyByTime = NaN(156,length(30:110),size(power,2));
                consistencyROI = NaN(size(power,2),156);
                
                % lets do pairwise sessions here
                for sess=1:size(power,2)-1
                for i=25:105
                    power_initial = squeeze(power{4,sess}(:,:,i-5:i+5,freqsOfInterest));
                    power_initial = squeeze(nanmean(power_initial,3)); %avg across time
                    power_initial = squeeze(nanmean(power_initial,3)); %avg across freqs
                    power_final = squeeze(power{4,sess+1}(:,:,i-5:i+5,freqsOfInterest));
                    power_final = squeeze(nanmean(power_final,3)); %avg across time
                    power_final = squeeze(nanmean(power_final,3)); %avg across freqs
                    minNumWords = min([size(power_initial,1) size(power_final,1)]);
                    for word=1:minNumWords
                        initial_pattern = power_initial(word,:);
                        final_pattern = power_final(word,:);
                        consistencyByTime(word,i-24,sess) = (1-pdist2(initial_pattern,final_pattern,'cosine'));
                    end
                end
                % now, lets incorporate correct-vs-incorrect in choosing what
                % to calculate similarity for
                tStats = pointsOfInterest{studyType,location,n};
                [xFind,yFind] = find((tStats) > tStatThresh);
                % eliminating certain frequencies and time points
                %             yFind(find(xFind<30)) = []; yFind(find(xFind>110)) = [];
                %             xFind(find(xFind<30)) = []; xFind(find(xFind>110)) = [];
                %             xFind(find(yFind<60)) = []; yFind(find(yFind<60)) = [];
                power_initial = squeeze(power{4,sess}(:,:,xFind,yFind));
                power_final = squeeze(power{4,sess+1}(:,:,xFind,yFind));
                minNumWords = min([size(power_initial,1) size(power_final,1)]);
                for word=1:minNumWords
                    initial_pattern = power_initial(word,:,:,:);
                    initial_pattern = reshape(initial_pattern,1,size(initial_pattern,2)*size(initial_pattern,3)*size(initial_pattern,4));
                    final_pattern = power_final(word,:);
                    final_pattern = reshape(final_pattern,1,size(final_pattern,2)*size(final_pattern,3)*size(final_pattern,4));
                    consistencyROI(sess,word) = (1-pdist2(initial_pattern,final_pattern,'cosine'));
                end       
                end
                consistencyROI = squeeze(nanmean(consistencyROI,1));
                consistencyByTime = squeeze(nanmean(consistencyByTime,3));
                
            end
            
            
            % lets find relationship between consistency and performance,
            % for all time points
            slope_high = NaN(size(power,2),size(consistencyByTime,2));
            slope_low = NaN(size(power,2),size(consistencyByTime,2));

            
            for sess=1:size(power,2)
            for time=1:size(consistencyByTime,2)
                consistencyThisTime = consistencyByTime(:,time);
                
                lowCutOff = prctile(consistencyThisTime,33);
                highCutOff = prctile(consistencyThisTime,66);
                
                lowCutOff = 0;
                highCutOff = 0;
                
                highConsistency = (find(consistencyThisTime > highCutOff));
                lowConsistency = (find(consistencyThisTime < lowCutOff));
                
                performance_initial = behavior{1,sess}; performance_initial(find(performance_initial==-999)) = 0;
                highConsistency(find(highConsistency>length(performance_initial))) = []; % remove words which werent shown in this session
                highConsistency_performance_initial = nansum(performance_initial(highConsistency)) / length(find(~isnan(performance_initial(highConsistency))));
                lowConsistency(find(lowConsistency>length(performance_initial))) = []; % remove words which werent shown in this session
                lowConsistency_performance_initial = nansum(performance_initial(lowConsistency)) / length(find(~isnan(performance_initial(lowConsistency))));
                
                performance_second = behavior{1,sess}; performance_second(find(performance_second==-999)) = 0;
                highConsistency(find(highConsistency>length(performance_second))) = []; % remove words which werent shown in this session
                highConsistency_performance_second = nansum(performance_second(highConsistency)) / length(find(~isnan(performance_second(highConsistency))));
                lowConsistency(find(lowConsistency>length(performance_second))) = []; % remove words which werent shown in this session
                lowConsistency_performance_second = nansum(performance_second(lowConsistency)) / length(find(~isnan(performance_second(lowConsistency))));
                
                slope_high(sess,time) = (highConsistency_performance_second - highConsistency_performance_initial);
                slope_low(sess,time) = (lowConsistency_performance_second - lowConsistency_performance_initial);
            end
            end
            
            
            
            
            gcf = figure('Units', 'Normalized', 'OuterPosition', [0 0 0.8 0.8]);
            clf
            subplot(2,2,1)
            slope_subtraction = slope_high - slope_low;
            plot(1:size(slope_high,1),slope_subtraction)
            xticks([0 20 40 60 80])
            xticklabels([0 1 2 3 4])
            xlabel('time (s)')
            ylabel('performance post-pre')
            ylim([-0.5 0.5])
            set(gca, 'FontName', 'Helvetica');
            set(gca,'FontSize',18)
            
            max_point = find(slope_subtraction==max(slope_subtraction));
            time = max_point;
            hold on
            plot(time,slope_subtraction(time),'r.')
            
            subplot(2,2,2)
            
            consistencyThisTime = consistencyByTime(:,time);
            lowCutOff = prctile(consistencyThisTime,33);
            highCutOff = prctile(consistencyThisTime,66);
            
            lowCutOff = 0;
            highCutOff = 0;
            
            highConsistency = (find(consistencyThisTime > highCutOff));
            lowConsistency = (find(consistencyThisTime < lowCutOff));
            
            performance_initial = behavior{1,1}; performance_initial(find(performance_initial==-999)) = 0;
            highConsistency_performance_initial = nansum(performance_initial(highConsistency)) / length(find(~isnan(performance_initial(highConsistency))));
            lowConsistency_performance_initial = nansum(performance_initial(lowConsistency)) / length(find(~isnan(performance_initial(lowConsistency))));
            
            performance_second = behavior{1,end}; performance_second(find(performance_second==-999)) = 0;
            highConsistency_performance_second = nansum(performance_second(highConsistency)) / length(find(~isnan(performance_second(highConsistency))));
            lowConsistency_performance_second = nansum(performance_second(lowConsistency)) / length(find(~isnan(performance_second(lowConsistency))));
            
            
            bar(1.2:2.2,[highConsistency_performance_initial highConsistency_performance_second],0.2,'r')
            hold on
            bar(1:2,[lowConsistency_performance_initial lowConsistency_performance_second],0.2,'k')
            legend('high consistency', 'low consistency')
            title(sprintf('at max difference',max_point));
            ylim([0 1])
            ylabel('performance'); xlabel('sessions');
            set(gca, 'FontName', 'Helvetica');
            set(gca,'FontSize',18)
            
            subplot(2,2,3)
            temp = consistencyByTime(:,time);
            hist(temp,20);
            title(sprintf('thresholds: below %0.2f, above %0.2f',lowCutOff,highCutOff))
            ylabel('counts')
            xlabel('cosine similarity, -1 to 1')
            set(gca, 'FontName', 'Helvetica');
            set(gca,'FontSize',18)
            
            
            
            % now, lets find relationship between performance and
            % similarity, for time points that showed correct vs incorrect
            % differences
            subplot(2,2,4)
            
            
            highConsistency = (find(consistencyROI > 0));
            lowConsistency = (find(consistencyROI  < 0));
            
             performance_initial = behavior{1,1}; performance_initial(find(performance_initial==-999)) = 0;
            highConsistency_performance_initial = nansum(performance_initial(highConsistency)) / length(find(~isnan(performance_initial(highConsistency))));
            lowConsistency_performance_initial = nansum(performance_initial(lowConsistency)) / length(find(~isnan(performance_initial(lowConsistency))));
            
            performance_second = behavior{1,end}; performance_second(find(performance_second==-999)) = 0;
            highConsistency_performance_second = nansum(performance_second(highConsistency)) / length(find(~isnan(performance_second(highConsistency))));
            lowConsistency_performance_second = nansum(performance_second(lowConsistency)) / length(find(~isnan(performance_second(lowConsistency))));
            
            slope_high = (highConsistency_performance_second - highConsistency_performance_initial);
            slope_low = (lowConsistency_performance_second - lowConsistency_performance_initial);
            
            bar(1.2:2.2,[highConsistency_performance_initial highConsistency_performance_second],0.2,'r')
            hold on
            bar(1:2,[lowConsistency_performance_initial lowConsistency_performance_second],0.2,'k')
            legend('high consistency', 'low consistency')
            title(sprintf('all correct-vs-incorr tStat>%0.2f. performance diff %0.2f',tStatThresh,slope_high-slope_low))
            set(gca, 'FontName', 'Helvetica');
            set(gca,'FontSize',18)
            
            st = suptitle(sprintf('subj %s, location %s, study type %s',subjectNames{n},locationStrings{location},evTypeStrings{studyType}));
            set(st,'FontSize',18)
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