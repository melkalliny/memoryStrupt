
pars = consolidation_setParams();
expType = 1;

subjects = [49 50 51 52 53 55 58];
subjectNames = {'NIH049', 'NIH050', 'NIH051', 'NIH052', 'NIH053', 'NIH055', 'NIH058'};

%  subjects = [29 32 36 41 42 46 47];
%  subjectNames = {'NIH029', 'NIH032', 'NIH036', 'NIH041', 'NIH042', 'NIH046', 'NIH047', 'NIH050', 'NIH051'};


npats = [1:9];

for n = [1]
   
    fprintf('\npreparing to do %s\n',subjectNames{n})
    
    load(sprintf('/Volumes/Shares/FRNU/dataWorking/consolidationProject/intermediates/power_allFrequencies/%s/powDir/allPowData.mat',subjectNames{n}))
    %load(sprintf('/Volumes/FRNU_SVN/consolidation_SVN/consolidationProject/intermediates/power/%s/powDir/pars.mat',subjectNames{n}))
    load(sprintf('/Volumes/Shares/FRNU/dataWorking/consolidationProject/intermediates/power_allFrequencies/%s/powDir/events.mat',subjectNames{n}))
    % apply the same events selection as the one used to create the power
    evsOfInterest = evsPull;
    
    if size(powDataAll,1) ~= length(evsOfInterest)
        fprintf('\nmismatch between events and power\n')
        keyboard
    end
    
    % starting indices of 100ms time bins: -1500:50:4400
    % lets do 0-2s (30-70) and 1-3s (50-90)
    
    % currently, events x bipolars x freqs x time bins
    powDataAll = permute(powDataAll,[1 2 4 3]); %avg across freqs
    powDataAllZscored = NaN(size(powDataAll,1),size(powDataAll,2),1);
    powDataAllBaseline = NaN(size(powDataAll,1),size(powDataAll,2),41);
    
        
    powDataAllMeanPrestim = nanmean(powDataAll(:,:,10:30),3);
    powDataAllMeanPoststim = nanmean(powDataAll(:,:,50:90),3);

    
    uniqueSessions = unique([evsOfInterest.session]);
    
    
    % get data z-scored to baseline
    eventsToGrabBaselineFrom = find(strcmp({evsOfInterest.type},'PROBE_START'));
    for sess = 1:length(uniqueSessions)
        eventsInSession = find([evsOfInterest.session] == uniqueSessions(sess));
        eventsBaselineInSession = intersect(eventsToGrabBaselineFrom,eventsInSession);
        for elec=1:size(powDataAll,2)
            for freq=1:size(powDataAll,4)
                meanDistrib = nanmean(nanmean((powDataAll(eventsBaselineInSession,elec,10:30,freq))));
                stdDistrib = nanmean(nanstd(powDataAll(eventsBaselineInSession,elec,10:30,freq)));
                for ev=1:length(eventsInSession)
                    powDataAllZscored(eventsInSession(ev),elec,31:119,freq) = ((powDataAll(eventsInSession(ev),elec,31:119,freq)) - meanDistrib) / stdDistrib;
                end
            end
        end
    end

    
    
    
    % lets go ahead and load the bad_chans.mat
    dirEcog = evsOfInterest(1).eegfile;
    dirEcog = strrep(dirEcog,'noreref','processed');
    dirBadChans = fullfile(dirEcog,'bad_chans.mat');
    load(dirBadChans)
    
    % average across electrodes separately for hippocampus,
    % parahippocampal gyrus/entorhinal, and temporal cortex
    % so, first, pull all the bipolar pairs that make up the powDataAll
     elecDir = sprintf('/Volumes/Shares/FRNU/dataWorking/consolidationProject/data/%s/docs/leads_bp.txt',subjectNames{n});
    electrodes = cell(1,1); iter = 1;
    fid=fopen(elecDir);
    if fid < 0; fprintf('\nno leads_bp.txt in /docs\n'); keyboard; end
    tline = fgetl(fid);
    while ischar(tline)
        electrodes{iter,1} = tline;
        tline = fgetl(fid);
        iter = iter+1;
    end
    fclose(fid);
    
    % now, match each of these up with elementInfo
    elementInfo = csv2cell(sprintf('/Volumes/Shares/FRNU/dataWorking/consolidationProject/data/%s/docs/element_info.csv',subjectNames{n}));
    atlasInfo = csv2cell(sprintf('/Volumes/Shares/FRNU/dataWorking/consolidationProject/data/%s/tal/atlas/atlas_bipolar_simple.csv',subjectNames{n}));
    % lets load the tal atlas
    
    % crude labeling for now
    hippocamp_elements = atlasInfo(find(contains(atlasInfo(:,6),'hippocampus','IgnoreCase',true)),1);
    MTL_elements = atlasInfo(union(find(contains(atlasInfo(:,4),'parahippocampal','IgnoreCase',true)),find(contains(atlasInfo(:,4),'entorhinal','IgnoreCase',true))),1);
    entorhinal_elements = atlasInfo(find(contains(atlasInfo(:,4),'entorhinal','IgnoreCase',true)),1);
    phc_elements = atlasInfo(find(contains(atlasInfo(:,4),'parahippocampal','IgnoreCase',true)),1);
    %lateral_elements = atlasInfo(find(contains(atlasInfo(:,4),'inferior temporal gyrus','IgnoreCase',true)),1);
    lateral_elements = atlasInfo(union(union(find(contains(atlasInfo(:,4),'inferior temporal gyrus','IgnoreCase',true)),find(contains(atlasInfo(:,4),'middle temporal gyrus','IgnoreCase',true))),find(contains(atlasInfo(:,4),'superior temporal gyrus','IgnoreCase',true))),1);
    frontal_elements = atlasInfo(find(contains(atlasInfo(:,4),'frontal gyrus','IgnoreCase',true)),1);
    
    hippocamp_bipolars = [];
    entorhinal_bipolars = [];
    lateral_bipolars = [];
    frontal_bipolars = [];
    parahippocampal_bipolars = [];
     % find which group each bipolar pair belongs to (if any)
     % only looking at first electrode of each pair. assuming that they are
     % both from the same anatomical location
    for elec=1:size(electrodes,1)
        findDash = strfind(electrodes{elec},'-');
        elecOne = electrodes{elec}(1:findDash-1);
        elecTwo = electrodes{elec}(findDash+1:end);
        
        elecOne(find(isstrprop(elecOne,'digit'))) = [];
        
        % while we're here, lets load the bad_chans.mat and not allow any
        % chans marked in bad_chans.mat to pass this step
        find_bad_chan_one = length(find(contains(bad_chans,elecOne)));
        find_bad_chans_two = length(find(contains(bad_chans,elecTwo)));
        if find_bad_chan_one + find_bad_chans_two > 0
            continue
        end
        
        hippocampal_elecs = find(contains(hippocamp_elements,electrodes{elec},'IgnoreCase',true));
        if size(hippocampal_elecs,1) ~= 0
            hippocamp_bipolars = cat(2,hippocamp_bipolars,elec);
            continue
        end
        entorhinal_elecs = find(contains(entorhinal_elements,electrodes{elec},'IgnoreCase',true));
        if size(entorhinal_elecs,1) ~= 0
            entorhinal_bipolars = cat(2,entorhinal_bipolars,elec);
            continue
        end
        lateral_elecs = find(contains(lateral_elements,electrodes{elec},'IgnoreCase',true));
        if size(lateral_elecs,1) ~= 0
            lateral_bipolars = cat(2,lateral_bipolars,elec);
            continue
        end
        frontal_elecs = find(contains(frontal_elements,electrodes{elec},'IgnoreCase',true));
        if size(frontal_elecs,1) ~=0
            frontal_bipolars = cat(2,frontal_bipolars,elec);
            continue
        end
        parahippocampal_elecs = find(contains(phc_elements,electrodes{elec},'IgnoreCase',true));
        if size(parahippocampal_elecs,1) ~=0
            parahippocampal_bipolars = cat(2,parahippocampal_bipolars,elec);
            continue
        end        
        
    end
    
   % aggregate data and save it out
   processed_power = cell(5,6);
   processed_power{5,1} = 'hippocamp'; processed_power{5,2} = 'entorhinal';
   processed_power{5,3} = 'lateralCortex'; processed_power{5,4} = 'frontalCortex';
   processed_power{5,5} = 'parahippocampal';
   
   processed_power{1,1} = powDataAllZscored(:,hippocamp_bipolars,:,:);
   processed_power{1,2} = powDataAllZscored(:,entorhinal_bipolars,:,:);
   processed_power{1,3} = powDataAllZscored(:,lateral_bipolars,:,:);
   processed_power{1,4} = powDataAllZscored(:,frontal_bipolars,:,:);
   processed_power{1,5} = powDataAllZscored(:,parahippocampal_bipolars,:,:);
   processed_power{1,6} = '1-3 s, Each Event Z-Scored To Prestim';

   
   
   if ~exist(sprintf('/Volumes/Shares/FRNU/dataWorking/consolidationProject/intermediates/processedPower_allFrequencies_responses/%s/',subjectNames{n}))
       mkdir(sprintf('/Volumes/Shares/FRNU/dataWorking/consolidationProject/intermediates/processedPower_allFrequencies_responses/%s/',subjectNames{n})); 
   end
   
   save(sprintf('/Volumes/Shares/FRNU/dataWorking/consolidationProject/intermediates/processedPower_allFrequencies_responses/%s/allPowData.mat',subjectNames{n}),'processed_power','evsOfInterest','-v7.3')
    
   fprintf('\ncompleted %s\n',subjectNames{n})

   
end