function [labelJW, cache] = getLabelJW(subj, rootEEGdir, chanName, cache)
%
%  Load and convert new tal LUTmp to old format
%
cacheEmpty = isempty(cache);
if nargin < 4 || cacheEmpty
    cache.bdc = [];
    %cache.bdf = [];
    cache.LUTbp = [];
    cache.leadsXYZ = [];
    cache.LUTmp = [];
    cache.atlas = [];
    cache.bd = [];
end

bdc      = cache.bdc;
LUTbp    = cache.LUTbp;
leadsXYZ = cache.leadsXYZ;
LUTmp    = cache.LUTmp;
atlas    = cache.atlas;
bd       = cache.bd;


% in patients with ONLY depths, there is no lead to ROI LUT, so the JW label cant be identified
roiLUTbp = fullfile(rootEEGdir, subj, '/tal/roi/lead_ROIC_LUT_bipolar.mat');
roiLUTmp = fullfile(rootEEGdir, subj, '/tal/roi/lead_ROIC_LUT_monopolar.mat');
if ~exist(roiLUTbp,'file') | ~exist(roiLUTmp,'file'),
    if cacheEmpty,
        %- first call to this function for a subject cache will be empty
        fprintf('\n heads up... cant compute JW column of atlas simple if no ROIC LUT');
        %keyboard
    end
    labelJW = 'Missing ROIC LUT';
    return;
end


% load Colin standard brain... this is what the JW measurements are based on
if isempty(bdc) || length(bdc.surf)==0,
    bdc = braindata;
    bdc.loadStandard();  %- colin = standard (old school)
end


% Load Mesh indices
if isempty(LUTbp)
    temp = load(fullfile(rootEEGdir, subj, '/tal/roi/lead_ROIC_LUT_bipolar.mat'));
    LUTbp = temp.lead_roi_lut_bp;
end
assert(strcmpi(LUTbp.subjId{1}, subj));

if isempty(LUTmp) || isempty(leadsXYZ)
    leadsXYZ = readtableSafe(fullfile(rootEEGdir, subj, '/tal/leads.csv'));
    temp     = load(fullfile(rootEEGdir, subj, '/tal/roi/lead_ROIC_LUT_monopolar.mat'));
    LUTmp    = temp.lead_roi_lut;
end

if isempty(bd)
    bd = braindata();
    bd.loadSubject(subj, rootEEGdir, 'hemi','separate')
end


%monopolar/bipolar
if contains(chanName, '-')
    lut = LUTbp;
    atlasFile = fullfile(rootEEGdir, subj, '/tal/atlas/atlas_bipolar.csv');
else
    lut = LUTmp;
    atlasFile = fullfile(rootEEGdir, subj, '/tal/atlas/atlas_monopolar.csv');
end
atlas = readtableSafe(atlasFile);


%- confirm existance of this channel
temp = atlas(strcmpi(atlas.chanName, chanName), :);
assert(height(temp) >= 1, 'Did not find a row in atlas table with chanName %s', chanName);
chanHemi = char(temp{1,'hemisphere'});
if isempty(chanHemi)
    e = strsplit(chanName, '-');
    e = char(e(1));
    jack = getJackTable(subj, rootEEGdir);
    chanHemi = char(jack{strcmpi(jack.chanName, e), 'whichHemi'});
end
assert(~isempty(chanHemi));



% find "center" of electrode in standard space so simple xyz can be saved and used for anterior/posterior decision (ATL vs PTL)
%- first look to see if channel is in LUT, if not (case for depths), use leadsCSV value.
irow = find(strcmp(lut.chanName,chanName));
if length(irow)==1
    
    % load channel here
    lutrow = lut(irow, :);
    mesh_ind = lutrow.nearest_mesh_ind{1};
    
    % Mesh_ind --> xyz
    xyzC = bdc.surf.(chanHemi).vertices(mesh_ind, :);
    %xyzF = bdf.surf.(chanHemi).vertices(mesh_ind, :);
    xyzS = lutrow(1,{'x','y','z'});
    
    % attempt to boil down to a single point in xyz space
    % ...find center of biggest grouping in xyzC space...
    % decrease neighbor radius until there is a single point with
    % more neighbors in this radius than any other point
    xyzCenter=[]; xyzMax=[]; maxNeighbors=0;
    for rad=[5:-0.2:.1],
        idx = rangesearch(xyzC,xyzC,rad);
        idxL = cellfun('length',idx);
        [maxVal maxI] = max(idxL);
        if maxVal>maxNeighbors,
            maxNeighbors = maxVal;
            xyzMax = mean(xyzC(idx{maxI},:),1);
        end
        %fprintf('\n rad=%.1f: maxidxL = %d', rad, max(idxL)); disp(idxL');
        if sum(idxL==max(idxL))==1 && max(idxL)>2,
            xyzCenter = xyzC(idxL==max(idxL),:);
            %fprintf('\n %s:::chan %s (%d/%d): center next to %d of %d points (rad=%.2f)',subj,chanName,irow,numChan, max(idxL),length(idxL),rad);
            break;
        end
    end
    if isempty(xyzCenter),
        xyzCenter = xyzMax;
        %fprintf('\n %s:::chan %s (%d/%d): didnt find center, using max',subj,chanName,iChan,numChanJCK);
        %keyboard;
    end
    
else
    
    % Take a coronal slice of the surface. For depths, find the nearest pial mesh index on this slice.
    % *assume* that brain is oriented such that post/anterior axis is || to y axis
    MARGIN = 2.5; % mm
    xyz_surf = bd.surf.(chanHemi);
    ysurf = xyz_surf.vertices(:,2);
    
    if contains(chanName, '-')
        bipo = strsplit(chanName, '-');
        ileads = find(ismember(leadsXYZ.chanName, bipo));
    else
        ileads = find(strcmp(leadsXYZ.chanName,chanName));
    end
    
    xyzS = mean(leadsXYZ{ileads, {'x','y','z'}}, 1); % (euclidean mean for bp)
    slice_mask = (xyzS(2) - MARGIN <= ysurf) & (ysurf <= xyzS(2) + MARGIN);
    slice_ndx = find(slice_mask);
    [i_ndx_slice, d] = knnsearch(xyz_surf.vertices(slice_ndx,:), xyzS);
    mesh_ind = slice_ndx(i_ndx_slice);
    xyzC = bdc.surf.(chanHemi).vertices(mesh_ind,:);
    
    xyzCenter = xyzC;
end


%-find atlas info for this channel
iAtlas = find(strcmp(atlas.chanName,chanName) & contains(atlas.atlas,'aparc.a2009s+aseg'));
if isempty(iAtlas),
    iAtlas = find(strcmp(atlas.chanName,chanName));
end
if length(iAtlas)>1,
    iAtlas = iAtlas(1); %- take the first index (random pick here)
end
if ~isempty(iAtlas),
    %- confirm that lobe is a string not just "nan"
    if ~iscell(atlas.lobe) && sum(isnan(atlas.lobe))==length(atlas.lobe),
        if find(strcmp(atlas.chanName,chanName))==1,
            fprintf('\n WARNING, atlas.lobe is all NAN in %s',atlasFile);
            keyboard;
        end
        chanLobe = 'empty lobes in atlas';
    else
        chanLobe = atlas.lobe{iAtlas};
    end
else
    chanLobe = 'not found in atlas';
end


%- convert hemisphere to more readable version
if     strcmp(chanHemi,'lh'), hemiStr = 'Left ';
elseif strcmp(chanHemi,'rh'), hemiStr = 'Right';
else   fprintf('\n %s --> not lt or rt', chanHemi); hemiStr = 'Inter'; keyboard; end


%- convert lobe to more readable version
if     strcmp(chanLobe,'temporal'),  labelJW = 'Temporal Lobe';
elseif strcmp(chanLobe,'parietal'),  labelJW = 'Parietal Lobe';
elseif strcmp(chanLobe,'frontal'),   labelJW = 'Frontal Lobe';
elseif strcmp(chanLobe,'occipital'), labelJW = 'Occipital Lobe';
elseif strcmp(chanLobe,'cingulate'), labelJW = 'Cingulate Cortex';
elseif strcmpi(chanLobe, 'insula'),  labelJW = 'Insular Cortex';
elseif isempty(chanLobe) & (strcmp(atlas.hardwareType{iAtlas},'depth')|isempty(atlas.hardwareType{iAtlas})),
    labelJW =  'Sub-lobar (depth)';
else
    fprintf('chanLobe %s not found', char(chanLobe));
    labelJW = 'Unknown';
end



%-  Wittig et al 2018:  ATL vs PTL cuts to match the paper -- a 45 degree line
LOBE_SPLIT_TEMPORAL = 1;  LOBE_TRIM_FRONTAL = 1;

if LOBE_SPLIT_TEMPORAL & (strcmp(labelJW,'Temporal Lobe') || strcmp(labelJW,'Occipital Lobe')), %- push limbic into anterior/posterior
    shiftCut = 18; %- (18) positive value will move cut posterior. empirically picked by JW to match lesion
    if     xyzCenter(2)>=xyzCenter(3)-shiftCut,  %- use
        labelJW = 'Anterior Temporal Lobe';
    elseif xyzCenter(2)>=0.25*xyzCenter(3)-60,  %- cut very posterior electrodes
        labelJW = 'Posterior Temporal Lobe';
    else
        labelJW = 'Occipital Lobe'; %- new hack to better match the brodman grouping version
    end
    
end
if LOBE_TRIM_FRONTAL & (strcmp(labelJW,'Frontal Lobe'))
    if xyzCenter(2)>=5*xyzCenter(3)+8 | xyzCenter(2)>62,      %- multiplier takes line off of the 45degree angle.  <1 = more vertical
        labelJW = 'Anterior Frontal Lobe';
    elseif xyzCenter(2)<=0.6*xyzCenter(3)-12 | xyzCenter(3)>55,
        labelJW = 'Posterior Frontal Lobe';
    else
        labelJW = 'Frontal Lobe Attention';
    end
end



%- save outputs
cache.bdc = bdc;
cache.LUTbp = LUTbp;
cache.LUTmp = LUTmp;
cache.atlas = atlas;
cache.bd = bd;

return


end
