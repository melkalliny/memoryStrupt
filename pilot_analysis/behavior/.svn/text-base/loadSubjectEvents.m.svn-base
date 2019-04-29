function [subjs, subjData] = loadSubjectEvents(subjCell, typeCell)

%subjCell - cell array with subj names referring to those with a list info
%       in [project_data_path '/data/behavior/used_list_info]

%typeCell - cell array to filter list events by type. only list events
%       with a matching type to a element in typeCell will be kept

%subjs - retunrs list of subject names corresponding to subjData indices
%subjData - cell array, one index for each subj. Each cell contains a list x session cell array 
%       with list events for repeated lists, contains an empty matrix if
%       not a repeated list

% ex. [subjs, subjData] = loadSubjectEvents({'NIH050'}, {'STUDY_PAIR_START'})

global debug;
debug = 1;

global user_people_path;
global project_data_path;

%load all the common smaller helper funcs, too small to deserve a file
run([user_people_path '/code/helper/small_helpers.m']); 

list_info_path = [project_data_path '/data/behavior/used_list_info'];
events_path_string = [project_data_path '/data/%s/behavioral/palRam/events.mat'];



%check inputs
if nargin < 1
    subjCell = {'all'};
end

if nargin < 2
    typeCell = {'all'};
end

dprintf('*loadSubjectEvents()*\n');

dprintf('subjCell:');
dprintf(' %s ', subjCell{:});
dprintf('\n');

dprintf('typeCell:');
dprintf(' %s ', typeCell{:});
dprintf('\n');




%load list_infos 
dprintf('loading listInfos\n');
[subjs, list_infos] = loadAllListInfo(list_info_path);

%subselect subjects
if isequal(subjCell{1}, 'all')
    selected_list_infos = list_infos;
else
    selected_subjs_filter= cellFilterXbyY(subjs, subjCell);
    subjs = subjs(selected_subjs_filter);
    selected_list_infos = list_infos(selected_subjs_filter);
end




%populate subjData
num_listInfo = length(selected_list_infos);
subjData = cell(num_listInfo, 1);

for iSubj = 1:num_listInfo
    
    dprintf('loading events for subj: %s \n', subjs{iSubj});
   
    current_list_info = selected_list_infos{iSubj};
    
    subj_num_list = size(current_list_info, 1);
    subj_num_sess = size(current_list_info, 2);
    
    subjData{iSubj} = cell(subj_num_list, subj_num_sess);
    
    events_fpath = sprintf(events_path_string, subjs{iSubj});
    load(events_fpath);
    
    %fill list x session spots in the subjData
    for iList = 1:subj_num_list
        for iSess = 1:subj_num_sess
           
            current_sess = iSess - 1;
            current_list = iList -1;
            
            if current_list_info(iList, iSess) ~= 0
               
                session_filter = ([events.session] == current_sess);
                list_filter = ([events.list] == current_list);
                
                if ~isequal(typeCell{1}, 'all')
                   type_filter = cellFilterXbyY({events.type}, typeCell);
                else
                   type_filter = ones(length(events), 1);
                end
                
                [session_filter, list_filter, type_filter] = reorient(session_filter, list_filter, type_filter);
                
                subjData{iSubj}{iList, iSess} = events(session_filter & list_filter & type_filter);

            else
                subjData{iSubj}{iList, iSess} = [];
            end
            
        end
    end
    
end

if debug
    debug = 0;
end

end