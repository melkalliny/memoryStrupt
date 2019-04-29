function t_simple = atlas_simplify(subj, rootEEGdir, varargin)
% ATLAS_SIMPLIFY post-processes and simplifies the atlas_[bi|mono]polar.csv
%
% USAGE:
%   t_simple = atlas_simplify(subj, rootEEGdir, ...)
%
% DESCRIPTION:
%   Outputs and a simplified version of the atlas.csv files
%
%   Simplifications:
%       - One row per electrode
%       - Desikan-Killiany atlas (aparc+aseg_rank)
%       - Important columns
%
%
%
% INPUT:
%   subj
%   rootEEGdir
%
% OPTIONAL KEY-VALUE INPUT:
%   'allow_ties' - If true, if there are k "nearest" regions for an electrode, output k rows
%               for that electrode. If false (default), choose a random row
%   'fname_in'   - filepath of atlas to simplify (default subject's monopolar)
%   'fname_out'  - filepath of save file (default is filename_simple.csv in subject's tal/atlas directory)
%   'max_dist'   - maximum distance in millimeters to include in table (defualt 5)
%
% OUTPUT:
%   t_simple
%
% FILES_OUTPUT:
%   prompts the user whether they would like to save t_simple to the subject's tal/atlas directory
%
% REVISION HISTORY:
%   MST 03/18   - Created
%   MST 07/18   - Added "more human" formatting and columns
%   JHW 11/18   - softened max dist contstraint... max dist cuts out all rows for a channel then the max_dist constraint is lifted for that channel
%               - added some protection against "no aparc or aparc.2009s" entries found. this happened for 57,58,61. now also checks for +aseg in those cases
%   JHW 01/19   - removed "forced_save" option... always save without prompt if function is called
%
% See also: atlas_table, atlas_table_fsaverage, atlas_simplify

DBG = 0;

% =============== Input argument parsing =================
ip = inputParser;
ip.addParameter('allow_ties', false);
ip.addParameter('fname_in', []);
ip.addParameter('fname_out', []);
ip.addParameter('max_dist', 5);
results = ip.parse(varargin{:});

allow_ties  = ip.Results.allow_ties;
fname_in    = ip.Results.fname_in;
fname_out   = ip.Results.fname_out;
max_dist    = ip.Results.max_dist;
clear ip

if isempty(fname_in) || isempty(fname_out)
    dir_atlas = fullfile(rootEEGdir, subj, 'tal/atlas');
    assert(exist(dir_atlas,'dir') > 0, 'Atlas dir not found: %s\n', dir_atlas);
    
    if isempty(fname_in)
        fname_in = fullfile(dir_atlas, 'atlas_monopolar.csv');
    end
    if isempty(fname_out)
        [~,fname] = fileparts(fname_in);
        fname_out = fullfile(dir_atlas, [fname, '_simple.csv']);
    end
end

% Define hemisphere expansion
hems_abbr = {'lh','rh'};
hems_full = {'Left', 'Right'};
hem_expand = @(a) hems_full{contains(hems_abbr, a)};

% Get the original table
assert(exist(fname_in, 'file') > 0, 'File not found: %s\n', fname_in);
torig = readtable(fname_in);
VARS = {'chanName' 'hemisphere' 'lobe' 'dist_mm' 'atlas' 'label_trim', 'label', 'hardwareType'};
torig = torig(:, VARS);

torig = torig(sortBySubstring(torig.chanName, getTagNames(subj, rootEEGdir)),:);

simple = [];
cache = [];
ATLASES = {'aparc', 'aparc.a2009s', 'aseg'};

% For each channel, build a single row
chanNames = unique(torig.chanName, 'stable');
if contains(chanNames{1}, '-')
    % bipolar
    pairs = getLeads(subj, rootEEGdir, 'isBP',1);
    chanNames = intersect(chanNames, strcat(pairs(:,1),'-',pairs(:,2)), 'stable');
else
    % monopolar
    chanNames = intersect(chanNames, getLeads(subj, rootEEGdir, 'isBP',0), 'stable');
end
assert(~isempty(chanNames), 'Empty chanNames');
for c = 1:numel(chanNames)
    chan = chanNames{c};
    tchan = torig(ismember(torig.chanName, {chan}), :);
    if sum(tchan.dist_mm <= max_dist)>0, %- JHW 10/2018.  this was cutting out some channels that didn't meet minimum criterion
        tchan = tchan(tchan.dist_mm <= max_dist, :);
    else
        fprintf('\n warning... distance threshold lifted for channel %s (min dist=%.1f mm)',chan,min(tchan.dist_mm));
    end
    
    isDepth = height(tchan) > 0 && strcmpi(tchan.hardwareType{1}, 'depth');
    
    % do a quick pass to find how many equi-distant labels there are (ie how many rows we will need)
    num_dups = 1;
    if allow_ties
        for a = 1:numel(ATLASES)
            atlas = ATLASES{a};
            tchan_atlas = tchan(ismember(tchan.atlas, {atlas}), :);
            num_dups = max(num_dups, height(tchan_atlas));
        end
    end
    
    
    for d = 1:num_dups
        desikan_label = []; % used for lobe
        isMissing = 0;
        
        if d == 1
            % struct initializaiton (this determines the order of columns!)
            s = struct();
            s.chanName = '';
            s.hemisphere = '';
            s.lobe_desikan = '';
            s.label_desikan = '';
            s.dist_desikan = nan;
            
            s.label_destrieux = '';
            s.dist_destrieux = nan;
            
            s.label_aseg = '';
            s.dist_aseg = nan;
            
            s.label_wittig = '';
        else
            % copy first
            s(d) = s(1);
            
            if DBG
                fprintf('duplicate: %s\n', chan);
            end
        end
        
        for a = 1:numel(ATLASES)
            atlas = ATLASES{a};
            tchan_atlas = tchan(ismember(tchan.atlas, {atlas, [atlas, '_rank']}), :);
            % _rank suffix added. (note suma adds _rank, which is equivalent except for a the color code)
            
            if isempty(tchan_atlas) & ~isDepth
                warning('Chan %s missing atlas %s', chan, atlas);
                
                %- JHW 10/2018. probably just a hack... but added this to fix a bug with 57,59,61 where there were no rows with "pure" aparc
                tchan_atlas = tchan(ismember(tchan.atlas, {atlas, [atlas, '_rank'], [atlas, '+aseg']}), :);
                if ~isempty(tchan_atlas),
                    fprintf(' --> but found aseg version (can we used this?)'); %- JW added for NIH057... no entries without aseg
                else
                    isMissing = 1;
                end
            end
            if isempty(tchan_atlas)
                continue;
            end
            
            ndx = min(height(tchan_atlas), d);
            
            switch atlas
                case 'aparc.a2009s'
                    label = tchan_atlas{ndx,'label'};
                otherwise
                    label = tchan_atlas{ndx,'label_trim'};
            end
            
            
            % further simplifications
            %  - enforce lowercase labels
            %  - split label into an Aseg-column and a DK-column
            %  - strip 'left'/'right' from label
            
            label = lower(char(label));
            to_remove = {'left','right'};
            
            for j = 1:numel(to_remove)
                label = strrep(label, ['-' to_remove{j} '-'], '');
                label = strrep(label, [to_remove{j} '-'], '');
                label = strrep(label, ['-' to_remove{j}], '');
                label = strrep(label, to_remove{j}, '');
            end
            
            switch atlas
                case 'aparc' % Desikan
                    s(d).label_desikan = atlasFSTranslate('desikan', label, 'translation','full');
                    desikan_label = label;
                    s(d).dist_desikan = tchan_atlas{ndx,'dist_mm'};
                case 'aparc.a2009s' % Destrieux
                    s(d).label_destrieux = atlasFSTranslate('destrieux', label, 'translation','pretty');
                    s(d).dist_destrieux = tchan_atlas{ndx,'dist_mm'};
                case 'aseg'
                    s(d).label_aseg = atlasFSTranslate('aseg', tchan_atlas{ndx,'label'}, 'translation','pretty');
                    s(d).dist_aseg = tchan_atlas{ndx,'dist_mm'};
                otherwise
                    error('Add case: %s', atlas)
            end
            
            % populate
            s(d).chanName = chan;
            s(d).hemisphere = hem_expand(tchan_atlas{ndx,'hemisphere'});
            lobe = pretty_lobe(char(tchan_atlas{ndx,'lobe'}));
            if ~isempty(lobe)
                s(d).lobe_desikan = lobe;
            end
            assert(~isempty(s(d).chanName));
        end % atlas loop
        
        % Get lobe from Desikan atlas
        if ~isempty(desikan_label) && ~contains(lower(desikan_label), 'unknown')
            lobe = aparc2lobe(desikan_label);
            if ~isempty(lobe)
                s(d).lobe_desikan = pretty_lobe(lobe);
            end
        end
        
        % Wittig
        [s(d).label_wittig, cache] = getLabelJW(subj, rootEEGdir, chan, cache);
        
        
        
        if ~isMissing
            simple = cat(1, simple, s(d));
        end
    end % end dup loop
end % end chan_loop


%- this happens occasionally
if isempty(simple),
    fprintf('\n ERROR: no simple atlas created for subj %s',subj);
    keyboard;
    return;
end
simple = rmfield(simple, {'dist_desikan','dist_destrieux','dist_aseg'});
t_simple = struct2table(simple);


% write to file
writetable(t_simple, fname_out);
fprintf('\n Atlas simple successful, written to %s\n\n',fname_out');


end %- main function



%%% HELPER FUNCTIONS
function pretty = pretty_lobe(lobe)
lobe = strtrim(lobe);
pretty = lobe;
if isempty(lobe) || all(double(lobe))==0 || contains(lower(lobe), 'unknown')
    pretty = '';
    return
end
switch lobe
    case 'temporal'
        pretty = 'Temporal Lobe';
    case 'parietal'
        pretty = 'Parietal Lobe';
    case 'frontal'
        pretty = 'Frontal Lobe';
    case 'occipital'
        pretty = 'Occipital Lobe';
    case 'cingulate'
        pretty = 'Cingulate Cortex';
    case 'insula'
        pretty = 'Insular Cortex';
    otherwise
        warning('Make a new case for lobe %s', lobe)
        keyboard;
end
end
