%
%   JW Brain Plotting Wrapper
%
%



root =  '/Volumes/JW24TB/data24TB/eeg_new';   % This directory should contain subjects NIH050, NIH051, and NIH052
%root =  '/Volumes/jwGDRIVE/data24TB/eeg_new';   % This directory should contain subjects NIH050, NIH051, and NIH052
subj =  'NIH050';




%% ----- Loading a subject's data -----
%   - The ez methods are the easiest way to plot a subject's data. 
%   - They don't require a brainplotter object; they make one for you and return it to you!
%       

bd = braindata2(subj, root);    %- load the subject braindata object
bpSubj = brainplotter();            %- create a brainplotter object for subject's brain
bpSubj.loadSurface(bd, 'lh pial');  %- load left  pial surface... use the bd object to provide path to subject data
bpSubj.loadSurface(bd, 'rh pial');  %- load right pial surface

bpAve = brainplotter();            %- create a brainplotter object for average brain
bpAve.loadSurface(fullfile(bpAve.dataPath, 'fsaverage/SUMA'),'lh pial');
bpAve.loadSurface(fullfile(bpAve.dataPath, 'fsaverage/SUMA'),'rh pial');


NUM_BRAINS = 1; % testing
%NUM_BRAINS = 3; % jw standard


%-----------------------------------------------------------------------------------------------------------------------------------------------------------%
%%% SUBJECT BRAIN %%%%
figure(1); clf; set(gcf,'color','w','units','pixels','position',[ -1903         895        1440         420])
if NUM_BRAINS==3, hAxList = [axes('position',[0  0.01 .33 .9]) axes('position',[.33 0.01 .33 .9]) axes('position',[.66 0.01 .33 .9])]; %- left,bottom,wdith,height
else              hAxList = axes('position',[0  0.01 .9 .9]); end

thisBP = bpSubj;     %- define the surface being plotted
bpSubjAx = hAxList;  %- keep a handle to the subject brain axes for later plotting

cla
bd.ezplot(thisBP,hAxList(1));
thisBP.view('right');
%thisBP.setOpacity(0.25);  %- transparency

if NUM_BRAINS>1,
    bd.ezplot(thisBP,hAxList(2));
    thisBP.view('ventral');
    title(subj);
    
    bd.ezplot(thisBP,hAxList(3));
    thisBP.view('left');
else
    title(subj)
end
%-----------------------------------------------------------------------------------------------------------------------------------------------------------%



%-----------------------------------------------------------------------------------------------------------------------------------------------------------%
%%% FS AVERAGE %%%%
figure(2); clf; set(gcf,'color','w','units','pixels','position',[ -1903         400        1440         420])
if NUM_BRAINS==3, hAxList = [axes('position',[0  0.01 .33 .9]) axes('position',[.33 0.01 .33 .9]) axes('position',[.66 0.01 .33 .9])]; %- left,bottom,wdith,height
else              hAxList = axes('position',[0  0.01 .9 .9]); end

thisBP = bpAve;       %- define the surface being plotted
bpAveAx = hAxList;    %- keep a unique handle to the brain axes for later plotting

cla
bd.ezplot(thisBP,hAxList(1));
thisBP.view('right');
%thisBP.setOpacity(0.25);  %- transparency

if NUM_BRAINS>1,
    bd.ezplot(thisBP,hAxList(2));
    thisBP.view('ventral');
    title('FS Average');
    
    bd.ezplot(thisBP,hAxList(3));
    thisBP.view('left');
else
    title(subj)
end
%-----------------------------------------------------------------------------------------------------------------------------------------------------------%



thisBP = bpSubj; axes(bpSubjAx(1));  thisBP.name = 'subj';
thisBP = bpAve;  axes(bpAveAx(1));   thisBP.name = 'ave';

%%% INDIVUDAL SUBJECT PLOT --- can use tal xyz to plot all electrodes, including DEPTHS
%- plot the electrodes color coded by grid/strip with lines between elements from same grid/strip. only makes sense on patient brain
PLOT_ELECTRODE_POINTS = 0;
if PLOT_ELECTRODE_POINTS,
    t = bd.tal.xyz;  %- plotting lines connecting elements of each grid only works if using an xyz table from tal... cant do that from clustered mesh
    thisBP.plotPointGroup(t, 'element_info',bd.docs.element_info,'radius',3.0);  %- use this if you want to color by hardware
    thisBP.legend('orientation','horizontal');
    thisBP.setOpacity(0.6); %- so you can see depths
 
    %thisBP.plotPoint(t, 'color', [0 1 1]); %- use this to color all the same
end

%- change the color of the entire brain --- darker gray
MAKE_BRAIN_DARKER = 0;
if MAKE_BRAIN_DARKER,
    thisBP.plotRegion([1:thisBP.stdNumVert],'color',[1 1 1]*.4,'blend',0.5,'surf','pial_lh'); %- only do this once... gets darker each time executed
    thisBP.plotRegion([1:thisBP.stdNumVert],'color',[1 1 1]*.4,'blend',0.5,'surf','pial_rh');
end

%- show atlas coloring
SHOW_ATLAS = 0;
if SHOW_ATLAS,
    if strcmp(thisBP.name,'subj'),
        atlasPath = [root '/' subj '/tal/zloc/freesurfer/' subj '/SUMA/std.141.lh.aparc.a2009s.annot.niml.dset'];
    else
        atlasPath = [thisBP.dataPath '/fsaverage/SUMA/std.141.lh.aparc.a2009s.annot.niml.dset'];
    end
    thisBP.plotAtlas(atlasPath,'surf','pial_lh')
    thisBP.plotAtlas(atlasPath,'surf','pial_rh')
end

%- now the mesh way -- cant connect the electrodes with lines.  and can't show DEPTHS
PLOT_ELECTRODES_MESH = 1;
if PLOT_ELECTRODES_MESH,
    t = bd.roi.lead_mesh;
    
    
    %- Plot a point at the mean of each cluster of mesh indicies
    LR = {'lh' 'rh'};
    for iLR=1:2,
        clust_thresh = 5.0;   %- (5) biases towards single point; 2.5 splits across gyri
        %thisBP.plotPoint(t(strcmp(t.whichHemi,LR{iLR}),:), 'color', [1 1 0], 'radius', 2.0, 'alpha', 1.0, 'mesh_clust_d', clust_thresh,'surf',['pial_' LR{iLR}]);
        thisBP.plotPointGroup(t(strcmp(t.whichHemi,LR{iLR}),:), 'radius', 2.0, 'alpha', 1.0, 'mesh_clust_d', clust_thresh,'surf',['pial_' LR{iLR}]);
    end
    %thisBP.clearPoints();
    
    
    %- plot a circle around every electrode
    for iElec=1:height(t),
        maxGeoDist = 2.5; %- 2.5 is good on subject brain; 5.0 on fsAve; mm BEYOND radius of electrode... 0 is same as "nearest_mesh_ind"
        thisBP.plotRegion(t.geodisc_mesh_ind{iElec}(find(t.geodisc_dist{iElec}<=maxGeoDist)), 'color',[1 1 1],'blend',1.0,'surf',['pial_' t.whichHemi{iElec}],'boundaryWidth',2); %- outline: blend=0; boundaryWidth>0
        %thisBP.plotRegion(t.geodisc_mesh_ind{iElec}(find(t.geodisc_dist{iElec}<=maxGeoDist)), 'color',[1 1 1],'blend',1.0,'surf','pial_lh','boundaryWidth',2); %- outline: blend=0; boundaryWidth>0
        
        maxGeoDist = 0.0; %- mm BEYOND radius of electrode... 0 is same as "nearest_mesh_ind"
        thisBP.plotRegion(t.geodisc_mesh_ind{iElec}(find(t.geodisc_dist{iElec}<=maxGeoDist)), 'color',[1 0 0],'blend',1.0,'surf',['pial_' t.whichHemi{iElec}],'boundaryWidth',0); %- fill: blend>0 (blend=opacity)
    end
    %thisBP.clearRegions();
    
end



return;

%- next step #1) plotting physiological effects
%- next step #2) combining physio effects into ROIs for an idividual subject, then combining across subjects


%% Plot data (subdural only)
bp = bd.ezplot();
xyz = bd.tal.xyz;
[~,isub] = bd.chans2subdural(xyz.chanName, bd.docs.jacksheet);
xyz = bd.tal.xyz(ismember(bd.tal.xyz.chanName, t.chanName), :);
n = height(xyz);
data = normrnd(0,2,n,1);
bp.plotPoint(xyz, 'color', data)
colorbar();

% Data to ROI
[VPerROI, valPerROI, sparseMask, rh_begin] = bd.leadData2ROI(data, xyz.chanName);

% Look at output
if 0
    whos valPerROI
    whos VPerROI
    whos sparseMask
    sum(sparseMask);
end
% Plot data at calculated ROIs
[vertexData, vertexCnt] = bp.plotRegionsData(VPerROI, valPerROI);

%% Cross-subject
bds = braindata2.xLoadSubj({'NIH050','NIH051','NIH052'}, root);
bda = braindata2();
bda.loadAverage();

% Create random data
datas = [];
xyzs = [];
chanNames = [];
for i=1:3
    xyz = bds(i).tal.xyz;
    [~,isub] = bd.chans2subdural(xyz.chanName, bds(i).docs.jacksheet);
    xyz = xyz(isub,:);
    datas{i} = normrnd(0,2,height(xyz),1);
    xyzs{i} = xyz;
    chanNames{i} = xyz.chanName;
end

% Average data
[xDataAvg, m] = braindata2.xLeadData2ROI(bds, datas, chanNames, 'minSubj',1);

% Combine hemispheres
nROI = length(m) / 2;
ml   = m(1:nROI);
mr   = m(1+nROI:end);

% Get index of rh start
rhi  = sum(ml) + 1;

% Get ROI vertices on fsaverage
VPerROI_lh = bda.roic2roi(bda.roi.ROIC_lh(ml), 'lh');
VPerROI_rh = bda.roic2roi(bda.roi.ROIC_rh(mr), 'rh');
VPerROI    = [VPerROI_lh; VPerROI_rh];

% Plot
bp = bda.ezplot();
bp.plotRegionsData(VPerROI, nanmean(xDataAvg(m,:),2), 'rh_begin',rhi);


