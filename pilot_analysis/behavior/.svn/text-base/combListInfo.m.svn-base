function combListInfo()

global user_people_path;
global project_data_path;

%load all the common smaller helper funcs, too small to deserve a file
run([user_people_path '/code/helper/small_helpers.m']); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%combListInfos%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%testForIdenticalSession.m looks at all the PAL sessions for 
% subjects in FRNU/eeg/data and creates outputs indicating which 
% sessions are indentical. This includes a subject_listInfo.mat file with 
% a (list x session) matrix with each row indicating lists that matched
% across sessions. However there are some subjects for which the
% repeated sessions do not include the first session. To detect those
% subjects, I ran testForIdenticalSession.m with varying skipSess values.
% Doing that, created multiple listInfos for each subject. This file looks
% through all the listInfos for any subject and finds the one representing
% the actual pattern of repeated sessions

listInfo_dir = [project_data_path '/data/behavior/all_list_info'];
finalListInfo_dir = [project_data_path '/data/behavior/used_list_info'];
checkDirValid(listInfo_dir);
checkDirValid(finalListInfo_dir);

bookkeeping_fname = 'used_listInfos.txt';
bookkeeping_fid = fopen([ finalListInfo_dir '/' bookkeeping_fname], 'a');

fprintf(bookkeeping_fid, '-----------------------------------------\n');
fprintf(bookkeeping_fid, '%s\n', datetime('now'));

%get all the listInfos
listInfo_dir_ls = dir(listInfo_dir);

%get the names of subjects with files in listInfo_dir

listInfo_dir_ls_names = { listInfo_dir_ls.name };
listInfo_dir_ls_unq_splits = unique(...
                                    cellfun(@(x) getCellEl(strsplit(x, '_'), 1) , listInfo_dir_ls_names, 'UniformOutput', 0));

%exclude random non NIH-prefixed files
listInfo_subjs = listInfo_dir_ls_unq_splits(...
                                            cellfun(@(x) length(x) >= 3 && isequal({'NIH'}, {x(1:3)}) , listInfo_dir_ls_unq_splits));


%check how many matches each list has across listInfos

for iSubj = 1:length(listInfo_subjs)

    current_subj = listInfo_subjs{iSubj};
    
    %which listInfos belong to this subject
    subj_listInfos_idx = cellfun(@(x) contains(x, current_subj), listInfo_dir_ls_names);
    subj_listInfos_fnames = listInfo_dir_ls_names(subj_listInfos_idx);
    
    numListInfos = sum(subj_listInfos_idx);
    
    %vector to hold the matchCount for each listInfo
    subj_listInfos = zeros(numListInfos, 1);
    
    for iInfo = 1:numListInfos
        
        load([ listInfo_dir '/' subj_listInfos_fnames{iInfo}]);
        
        crossSessMatch_perList = sum(subj_list_matches, 2);
        
        %sums across columns of only 1 indicate that no sessions contained
        %a matching list
        if any((crossSessMatch_perList - 1) > 0)
            subj_listInfos(iInfo) = sum(crossSessMatch_perList);
        else
            subj_listInfos(iInfo) = 0;         
        end
        
    end
    
    %find the listInfo index of the listInfo containing the most repeated sessions
    maxListInfo_idx = 0;
    
    if any(subj_listInfos)
        [~, maxListInfo_idx] = max(subj_listInfos);
    end
    
    %copy those listInfos
    if maxListInfo_idx ~= 0
        fprintf('subj: %s , max list info: %s\n', current_subj, subj_listInfos_fnames{maxListInfo_idx});
        
        fprintf(bookkeeping_fid, '%s --> %s\n', subj_listInfos_fnames{maxListInfo_idx}, [current_subj '.mat']);
        copyfile([ listInfo_dir '/' subj_listInfos_fnames{maxListInfo_idx}], [finalListInfo_dir '/' current_subj '.mat']);
    else
        fprintf('subj: %s , no matches of any kind\n', current_subj);        
    end
    
    
                                        
end

fclose(bookkeeping_fid);

end
