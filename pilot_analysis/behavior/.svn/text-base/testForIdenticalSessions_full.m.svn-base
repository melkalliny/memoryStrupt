clear all;
collectTrialInfo = 0;

% create the subj_beh_info.csv
% contains info on which subjects have repeated PAL sessions

beh_dir = '/Volumes/Shares/FRNU/data/eeg';

if ~exist(beh_dir, 'dir')
    error(['You are not connected to ' beh_dir]);
end

subj_info_output_dir = '/Volumes/Shares/FRNU/dataWorking/consolidationProject/data/behavior';
trial_info_output_dir = '/Volumes/Shares/FRNU/dataWorking/consolidationProject/data/behavior/trial_info';
list_info_output_dir = '/Volumes/Shares/FRNU/dataWorking/consolidationProject/data/behavior/all_list_info';

if ~exist(subj_info_output_dir, 'dir')
    error(['directory doesnt exist: ' subj_info_output_dir]);
end

if ~exist(trial_info_output_dir, 'dir')
    error(['directory doesnt exist: ' trial_info_output_dir]);
end

if ~exist(list_info_output_dir, 'dir')
    error(['directory doesnt exist: ' list_info_output_dir]);
end

skipSess_runs = [0 1 2 3];

for iSkipSess = 1:length(skipSess_runs)
    
    skipSess = skipSess_runs(iSkipSess);

    subj_info_fpath = [subj_info_output_dir sprintf('/subj_info_s%d.csv', skipSess)];

    subj_info_lines = {};
    subj_info_lines{1} = strjoin({'subj' ...
                                  'numSess' ...
                                  'listsPerSess' ...
                                  'trialsPerSess' ...
                                  'trialsInListPerSess' ...
                                  'varyingListLength' ...
                                  'minCommonListsAcrossSess' ...
                                  }, ',');

     for l = 0:25
        subj_info_lines{1} = [ subj_info_lines{1} ',' sprintf('list%d', l) ]; 
     end

    numSubjects = 10:60;

    for subjNum=1:length(numSubjects)

        current_subj = numSubjects(subjNum);

        fprintf('NIH0%.2d\n', current_subj);
        trial_info_fpath = [trial_info_output_dir sprintf('/NIH0%.2d_trial_info.csv', current_subj)];

        subj_events_fpath_PAL = sprintf('%s/NIH0%.2d/behavioral/palRam/events.mat',beh_dir, current_subj);
        subj_events_fpath_pa3 = sprintf('%s/NIH0%.2d/behavioral/pa3/events.mat',beh_dir, current_subj);

        % info to be output on each line of the csv 
        numSess = 0;
        listsPerSess = []; % array containing number of lists started in each session
        trialsPerSess = []; % array containing number of word pairs tested in each session
    %     trialsInListPerSess = []; % array containing the number trials that made up each list in each session
        trialsInListPerSess_str = '';
        varyingListLength = 0;

        commonListsAcrossSess = []; % the number of common lists across all sessions

        subj_list_array_matches = '';
        subj_list_matches = [];

        if exist(subj_events_fpath_PAL, 'file') && exist(subj_events_fpath_pa3, 'file')

            fprintf('subject did both PAL and pa3\n');
            keyboard;

        elseif exist(subj_events_fpath_PAL, 'file')

            subj_events_fpath = subj_events_fpath_PAL;

        elseif exist(subj_events_fpath_pa3, 'file')

            subj_events_fpath = subj_events_fpath_pa3;

        else

            subj_events_fpath = '';
        end

        % subject only gets a line in the csv if they have a PAL session
        if exist(subj_events_fpath, 'file')

            fprintf('%s\n', subj_events_fpath);
            load(subj_events_fpath);
            events_table = struct2table(events);
            events = events_table;

            %create table with all study_pair_start events AND exclude events
            %that have -999 in serialpos or probepos, these trials were not completed
            study_pair_filt = cellfun(@(x) isequal('STUDY_PAIR_START', x) || isequal('STUDY_PAIR', x)  , events.type);
            serialpos_filt = events.serialpos ~= -999;
            probepos_filt = events.probepos ~= -999;

            %apply all filters, make table
            all_studyPairs = events(study_pair_filt & serialpos_filt & probepos_filt , : );


            %get number of sessios in subj
            subjUniqSess = unique([events.session]);
            numSess = length(subjUniqSess);
            sessCellSize = max(subjUniqSess) + 1;

            sessCell = cell(sessCellSize, 1);

            studyPairs_bySess = cell(sessCellSize, 1);

            for iSess = 0:(sessCellSize-1)

                %fprintf('sorting study pairs for session %d\n', iSess);
                currentSess_studyPairs = all_studyPairs([all_studyPairs.session] == iSess, :);
                studyPairs_bySess{iSess+1} = currentSess_studyPairs;

                %session might have been skipped
                if ~isempty(currentSess_studyPairs)

                    currentSess_lists = unique(currentSess_studyPairs{:, 'list'});
                    currentSess_listCount = length(currentSess_lists);

                    listsPerSess(length(listsPerSess) + 1) = currentSess_listCount;
                    trialsPerSess(length(trialsPerSess) + 1) = length(currentSess_studyPairs{:, 'list'});

                    currentSess_maxList = max(currentSess_lists)+1;

                    sessCell{iSess+1} = cell(currentSess_maxList, 1);

                    for iList = 0:(currentSess_maxList-1)

                        %fprintf('\tsorting study pairs for session %d - list %d \n', iSess, iList);

                        currentSess_currentList_studyPairs = currentSess_studyPairs([currentSess_studyPairs.list] == iList, :);
                        sessCell{iSess+1}{iList+1} = currentSess_currentList_studyPairs;
                    end

                end

            end

            %are all the list the same length across sessions

            for iSess = 0:(length(sessCell)-1)

                if ~isempty(sessCell{iSess+1})
                    %fprintf('checking list lengths for session %d\n', iSess);

                    currentSess_listLengths = [];

                    for iList = 0:(length(sessCell{iSess+1})-1)

                        %fprintf('\tchecking list lengths for session %d - list %d : %d\n', iSess, iList, size(sessCell{iSess+1}{iList+1}, 1));            

                        if size(sessCell{iSess+1}{iList+1}, 1) > 0 
                            currentSess_listLengths(length(currentSess_listLengths) + 1) = size(sessCell{iSess+1}{iList+1}, 1);
                        end
                    end

                    unique_list_length = unique(currentSess_listLengths);

                    if length(unique_list_length) > 1
                        varyingListLength = 1;
                    end

                    trialsInListPerSess_str = [ trialsInListPerSess_str ' '];

                    for iLL = 1:length(unique_list_length)

                        LL = unique_list_length(iLL);

                        LL_count = sum(currentSess_listLengths == LL);

                        trialsInListPerSess_str = [ trialsInListPerSess_str sprintf('%d(%d)', LL, LL_count )];
                    end

    %                 if any(currentSess_listLengths - currentSess_listLengths(1))
    % 
    %                     fprintf('list length changed in the middle of a session\n');
    % %                     keyboard;
    %                     
    %                     trialsInListPerSess(length(trialsInListPerSess) + 1) = currentSess_listLengths(1) + 0.1;
    % 
    %                 end
    % 
    %                 trialsInListPerSess(length(trialsInListPerSess) + 1) = currentSess_listLengths(1);
    %                 
    %             else
    %                 trialsInListPerSess(length(trialsInListPerSess) + 1) = -1;
                end

            end

    %         trialsInListPerSess_realSess = trialsInListPerSess( trialsInListPerSess ~= -1 );

    %         if any(trialsInListPerSess - trialsInListPerSess(1))
    %             varyingListLength = 1;
    %         end


            %calculate minCommonListsAcrossSess

            for iList = 0:(max(listsPerSess)-1)

                listPresentInAllSess = 1;

                for iSess = 0:(length(sessCell)-1)

                    if ~isempty(sessCell{iSess+1})
                        if (iList+1) > length(sessCell{iSess+1}) || ...
                           ((iList+1) <= length(sessCell{iSess+1}) && isempty(sessCell{iSess+1}{iList+1}))
                            listPresentInAllSess = 0;
                        end
                    end
                end

                if listPresentInAllSess

                    commonListsAcrossSess(length(commonListsAcrossSess) + 1) = iList;
                end


                list_match_array = '['; % + if study1 and study2 remain the same, if study order is the same, and if cue direction is same, and if test order is differnet
                                       % - otherwise
                                       % # if no list exists for that session
                                       % 1 marks the list used as comparison
                list_match_bin = [];

                % find first session where list exists and check for subsequent
                % matches
                listToMatch = [];
                skippedSess = 0;

                for iSess = 0:(length(sessCell)-1)

                    if ~isempty(sessCell{iSess+1})
                        if  (iList+1) > length(sessCell{iSess+1}) || ...
                           ((iList+1) <= length(sessCell{iSess+1}) && isempty(sessCell{iSess+1}{iList+1}))

                            list_match_array = [ list_match_array ' # ' ];
                            list_match_bin = [ list_match_bin 0 ];
                            skippedSess = skippedSess + 1;

                        elseif isempty(listToMatch)  

                            if skippedSess >= skipSess
                                listToMatch = sessCell{iSess+1}{iList+1};
                                list_match_array = [ list_match_array ' 1 ' ];
                                list_match_bin = [ list_match_bin 1 ];                        
                            else
                               skippedSess = skippedSess + 1;
                               list_match_array = [ list_match_array ' s ' ];                            
                               list_match_bin = [ list_match_bin 0 ];                        

                            end

                        else

                                    %check those following lists to see if they match the
                                    %first one

                                    first_sess_study_1 = listToMatch{:, 'study_1'};
                                    first_sess_study_2 = listToMatch{:, 'study_2'};
                                    first_sess_serialpos = listToMatch{:, 'serialpos'};
                                    first_sess_probpos = listToMatch{:, 'probepos'};
                                    first_sess_probe_word = listToMatch{:, 'probe_word'};

                                    check_sess_study_1 = sessCell{iSess+1}{iList+1}{:, 'study_1'};
                                    check_sess_study_2 = sessCell{iSess+1}{iList+1}{:, 'study_2'};
                                    check_sess_serialpos = sessCell{iSess+1}{iList+1}{:, 'serialpos'};
                                    check_sess_probpos = sessCell{iSess+1}{iList+1}{:, 'probepos'};
                                    check_sess_probe_word = sessCell{iSess+1}{iList+1}{:, 'probe_word'};

                                    %check for exact match
                                    exact_check = isequal(first_sess_study_1, check_sess_study_1) ...
                                        && isequal(first_sess_study_2, check_sess_study_2) ...
                                        && isequal(first_sess_serialpos, check_sess_serialpos) ...
                                        && isequal(first_sess_probe_word, check_sess_probe_word);

                                    %check for non exact match
                                    %check if the pairs are in first_sess are in
                                    %check_sess, in any order and any pairing order
                                    pairs_present = [];

                                    for f1_idx = 1:length(first_sess_study_1)

                                        f1 = first_sess_study_1{f1_idx};
                                        f2 = first_sess_study_2{f1_idx};

                                        c1_match_idx = -1;
                                        c2_match_idx = -1;

                                        for c1_idx = 1:length(check_sess_study_1)
                                            if isequal(f1, check_sess_study_1{c1_idx})
                                                c1_match_idx = c1_idx;
                                            end
                                        end

                                        %if c1_match_idx ~= -1 then the word was
                                        %found and the pair should be in check_sess_study_2

                                        if c1_match_idx ~= -1

                                            for c2_idx = 1:length(check_sess_study_2)
                                                if isequal(f2, check_sess_study_2{c2_idx})
                                                    c2_match_idx = c2_idx;
                                                end
                                            end

                                        end

                                        %if c1_match_idx == -1 then the word was
                                        %not found and we should check if its in check_sess_study_2

                                        if c1_match_idx == -1

                                            for c2_idx = 1:length(check_sess_study_2)
                                                if isequal(f1, check_sess_study_2{c2_idx})
                                                    c2_match_idx = c2_idx;
                                                end
                                            end

                                            %if still no match, then c2_match_idx will
                                            %remain -1

                                            %if there is a match, the c2_match_idx will
                                            %be positive and we should look for the f2
                                            %match in c1

                                            if c2_match_idx ~= -1

                                                for c1_idx = 1:length(check_sess_study_1)
                                                    if isequal(f2, check_sess_study_1{c1_idx})
                                                        c1_match_idx = c1_idx;
                                                    end
                                                end
                                            end

                                        end

                                        if c1_match_idx ~= -1 && c2_match_idx ~= -1 && c1_match_idx == c2_match_idx
                                            pairs_present = [ pairs_present 1 ];
                                        else
                                            pairs_present = [ pairs_present 0 ];
                                        end

                                    end



                            non_exact_match = all(pairs_present);

                            if exact_check

                                list_match_array = [ list_match_array ' x ' ];
                                list_match_bin = [ list_match_bin 1 ];

                            elseif non_exact_match

                                list_match_array = [ list_match_array ' + ' ];
                                list_match_bin = [ list_match_bin 1 ];                      

                            else

                                list_match_array = [ list_match_array ' - ' ];
                                list_match_bin = [ list_match_bin 0 ];

                            end
                        end
                    end    
                end

                list_match_array = [ list_match_array ']' ];

                if isempty(subj_list_array_matches)
                    subj_list_array_matches = list_match_array;
                else  
                    subj_list_array_matches = [ subj_list_array_matches ',' list_match_array ];
                end

                subj_list_matches = [ subj_list_matches ; list_match_bin ];
             end


            if collectTrialInfo
                %measure the commonalities of trials across session

                trial_info_lines = {};

                current_trial_info_line = [ 'session' ',' 'study_1' ',' 'study_2' ',' 'list_num' ',' 'serialpos' ',' 'probepos' ',' 'probe_word' ]; 
                trial_info_lines{1} = current_trial_info_line;

                firstSess_studyPairs = studyPairs_bySess{1};

                current_listxSess = cell(length(sessCell), 1);
                current_list = -1;

                %loop over each study pair 
                for iPair = 1:size(firstSess_studyPairs, 1)

                    study1 = firstSess_studyPairs{iPair, 'study_1'};
                    study2 = firstSess_studyPairs{iPair, 'study_2'};

                    list_num = firstSess_studyPairs{iPair, 'list'};
                    serialpos = firstSess_studyPairs{iPair, 'serialpos'};
                    probepos = firstSess_studyPairs{iPair, 'probepos'};
                    probe_word = firstSess_studyPairs{iPair, 'probe_word'};

                    if current_list == list_num

                        firstPairInList = 0;
                        current_listxSess{1} = [ current_listxSess{1} ; firstSess_studyPairs(iPair, :)];
                    else

                        firstPairInList = 1;
                        current_list = list_num;
                        current_listxSess{1} = firstSess_studyPairs(iPair, :);
                    end

                    study_pair = [ study1{1} study2{1} ];
                    study_pair_reverse = [ study2{1} study1{1} ];

                    current_trial_info_line = [ num2str(subjUniqSess(1)) ',' study1{1} ',' study2{1} ',' num2str(list_num) ',' num2str(serialpos) ',' num2str(probepos) ',' probe_word{1} ]; 
                    trial_info_lines{length(trial_info_lines) + 1} = current_trial_info_line;

                    fprintf('\n');
                    fprintf('%s : trial %d of %d ( %s %s )', num2str(current_subj), iPair, size(firstSess_studyPairs, 1), study1{1}, study2{1} );

                    % check if trial exists in all following sessions, within
                    % existing lists

                    trialMatchesCommonSess = [];
                    trialMatchesAllSess = [];

                    for jSess = 1:(length(sessCell)-1)

                       currentSessNum = jSess;

                       fprintf(' ... checking session%d', currentSessNum);

                       trialMatchFound = 0; 

                       % assuming that lists will be in the same order
                       if ~isempty(sessCell{jSess+1}) 

                           if list_num < length(sessCell{jSess+1}) && ~isempty(sessCell{jSess+1}{list_num+1})

                               currentSess_studyPairs = studyPairs_bySess{jSess+1};

                               for jPair = 1:size(currentSess_studyPairs, 1)

                                   current_study1 = currentSess_studyPairs{jPair, 'study_1'};
                                   current_study2 = currentSess_studyPairs{jPair, 'study_2'};

                                   current_list_num = currentSess_studyPairs{jPair, 'list'};
                                   current_serialpos = currentSess_studyPairs{jPair, 'serialpos'};
                                   current_probepos = currentSess_studyPairs{jPair, 'probepos'};
                                   current_probe_word = currentSess_studyPairs{jPair, 'probe_word'};

                                   if current_list == current_list_num && firstPairInList
                                        %currentSess_studyPairs(jPair, :);
                                        current_listxSess{jSess+1} = [ current_listxSess{jSess+1} ; currentSess_studyPairs(jPair, :)];
                                   end

                                   current_study_pair = [ current_study1{1} current_study2{1} ];

                                   %found the matching trial in this session
                                   if isequal(study_pair, current_study_pair) || isequal(study_pair_reverse, current_study_pair) 

                                        trialMatchFound = 1;

                                        current_trial_info_line = [ num2str(currentSessNum) ',' current_study1{1} ',' current_study2{1} ',' num2str(current_list_num) ',' num2str(current_serialpos) ',' num2str(current_probepos) ',' current_probe_word{1} ]; 
                                        trial_info_lines{length(trial_info_lines) + 1} = current_trial_info_line;

                                   end
                               end

                               % no trial match was found 
                               if trialMatchFound == 0

                                    current_trial_info_line = [ num2str(currentSessNum) ',' 'missing' ',' 'missing' ',' 'missing' ',' 'missing' ',' 'missing' ',' 'missing' ]; 
                                    trial_info_lines{length(trial_info_lines) + 1} = current_trial_info_line;
                               end

                           else

                                noList_string = sprintf('no_list_%d', list_num);
                                current_trial_info_line = [ num2str(currentSessNum) ',' noList_string ',' noList_string ',' noList_string ',' noList_string ',' noList_string ',' noList_string ]; 
                                trial_info_lines{length(trial_info_lines) + 1} = current_trial_info_line;

                           end 

                       else

                            noList_string = sprintf('no_session_%d', jSess);
                            current_trial_info_line = [ num2str(currentSessNum) ',' noList_string ',' noList_string ',' noList_string ',' noList_string ',' noList_string ',' noList_string ]; 
                            trial_info_lines{length(trial_info_lines) + 1} = current_trial_info_line;

                       end


                    end

                    % add a blank line inbetween trials
    %                 trial_info_lines{length(trial_info_lines) + 1} = ' ';
                end

            end

            if collectTrialInfo

                %write out the trial info

                trial_info_fid = fopen(trial_info_fpath, 'w');

                for iLine = 1:length(trial_info_lines)
                    fprintf(trial_info_fid, '%s\n', trial_info_lines{iLine});
                end

                fclose(trial_info_fid);

            end

            subj_list_mat_fpath = [list_info_output_dir sprintf('/NIH0%.2d_s%d.mat', current_subj, skipSess)];
            save(subj_list_mat_fpath, 'subj_list_matches');

            next_subj_line = sprintf('NIH0%.2d', current_subj);
            next_subj_line = [next_subj_line ',' num2str(numSess)];
            next_subj_line = [next_subj_line ',' '[' num2str(listsPerSess) ']'];
            next_subj_line = [next_subj_line ',' '[' num2str(trialsPerSess) ']'];
            next_subj_line = [next_subj_line ',' '[' trialsInListPerSess_str ']'];
            next_subj_line = [next_subj_line ',' num2str(varyingListLength)];
            next_subj_line = [next_subj_line ',' shrinkNumArray(commonListsAcrossSess)];
            next_subj_line = [next_subj_line ',' subj_list_array_matches];

            subj_info_lines{length(subj_info_lines) + 1} = next_subj_line;


        end


    end

    subj_info_fid = fopen(subj_info_fpath, 'w');

    for iLine = 1:length(subj_info_lines)
        fprintf(subj_info_fid, '%s\n', subj_info_lines{iLine});
    end

    fclose(subj_info_fid);
end

function s = shrinkNumArray(arr)

    if isempty(arr)
       s = '[]';
       return;
    end

    s_num = arr(1);
    s = [ '[' num2str(arr(1)) ];

    last_consecutive_num = -1;
    
    for i = 2:length(arr)
       
        if arr(i) - s_num > 1
            
            if last_consecutive_num == -1
                s = [ s ',' num2str(arr(i))];
                s_num = arr(i);
            else
                s = [ s ':' num2str(last_consecutive_num) ',' arr(i)];
                last_consecutive_num = -1;
            end
            
        else
            last_consecutive_num = arr(i);
            s_num = arr(i);
        end
    end
    
    if last_consecutive_num ~= -1
        s = [ s ':' num2str(last_consecutive_num)];
    end
    
    s = [ s ']' ];
    
end

