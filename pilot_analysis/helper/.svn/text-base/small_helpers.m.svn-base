
checkDirValid = @checkDirValidFunc;
getCellEl = @getCellElFunc;
loadAllListInfo = @loadAllListInfoFunc;
cellFilterXbyY = @cellFilterXbyYFunc;
dprintf = @dprintfFunc;
reorient = @reorientFunc;

%takes in a cell array of 1d vectors and returns them all with dimensions
%aligned, all column vectors
function varargout = reorientFunc(varargin)

    varargout = cell(length(varargin), 1);
    
    for iCell = 1:length(varargin)
       
        [s1, s2] = size(varargin{iCell});
        if s2 > s1
           
            varargout{iCell} = varargin{iCell}';
        else
            varargout{iCell} = varargin{iCell};
        end
        
    end
    
end

%debugging print function, only prints when global debug is true
function dprintfFunc(varargin)
    
    global debug;

    if debug
        fprintf(varargin{:});
    end
end

%returns a binary vector indicating which elements of X are also in Y
function binary_res = cellFilterXbyYFunc(X, Y)
    binary_res = cellfun(@(x) any(cellfun(@(y) isequal(x, y) , Y)) , X);
end

%checks to make sure your directory exists, throws an error otherwise
function checkDirValidFunc(d, makeflag)
    
    if nargin < 2
        makeflag = 0;
    end
       
    if ~makeflag
        if ~exist(d, 'dir')
           error('%s is not an available directory', d); 
        end
    else
        mkdir(d);
    end
end

%returns the index of a cell array, useful with strsplit
function e = getCellElFunc(c, idx)
    e = c{idx};
end


function [subj_names, list_infos] = loadAllListInfoFunc(list_info_path)
    
    checkDirValidFunc(list_info_path);
    
    ls_dir = dir(list_info_path);
    
    all_file_names = {ls_dir.name};
    all_subj_names = cellfun( @(x) getCellElFunc(strsplit(x, '.'), 1) ,all_file_names, 'UniformOutput', 0);
    mat_file_filter = cellfun(@(x) contains(x, '.mat') , all_file_names);
    
    subj_names = all_subj_names(mat_file_filter);
    list_info_fnames = all_file_names(mat_file_filter);
    
    list_infos = cell(length(list_info_fnames), 1);
    
    for iFile = 1:length(list_info_fnames)
        loaded = load([list_info_path '/' list_info_fnames{iFile}]);
        list_infos{iFile} = loaded.subj_list_matches;
    end

end