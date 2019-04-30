function [out] = strupt_preprocessLowAndHigh_allTime()

out = [];
pars = strupt_setParams();
pars.freqs=cat(2,logspace(log10(4),log10(12),10),logspace(log10(80),log10(150),10));

expType = 1;

subjects = [71];
subjectNames = {'NIH071'};

for windowSize=1:3
    for n = subjects
        
        fprintf('\npreparing to do %s\n',subjectNames{n})
        
        load(sprintf('%s/%d/%s/powDir/allPowData.mat',pars.dirWriteOut,windowSize,subjectNames{n}))
        %load(sprintf('/Volumes/FRNU_SVN/consolidation_SVN/consolidationProject/intermediates/power/%s/powDir/pars.mat',subjectNames{n}))
        load(sprintf('%s/%d/%s/powDir/events.mat',pars.dirWriteOut,windowSize,subjectNames{n}))
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
        powDataAllZscored = NaN(size(powDataAll,1),size(powDataAll,2),size(powDataAll,3),size(powDataAll,4));
        
        uniqueSessions = unique([evsOfInterest.session]);
        
        baseline = {11:20,6:10,4:5};
        adjust = {21:111,11:51,6:21};
        
        % get data z-scored to baseline
        eventsToGrabBaselineFrom = find(strcmp({evsOfInterest.type},'PROBE_START'));
        for sess = 1:length(uniqueSessions)
            eventsInSession = find([evsOfInterest.session] == uniqueSessions(sess));
            eventsBaselineInSession = intersect(eventsToGrabBaselineFrom,eventsInSession);
            for elec=1:size(powDataAll,2)
                for freq=1:size(powDataAll,4)
                    meanDistrib = nanmean(nanmean((powDataAll(eventsBaselineInSession,elec,baseline{1,windowSize},freq))));
                    stdDistrib = nanmean(nanstd(powDataAll(eventsBaselineInSession,elec,baseline{1,windowSize},freq)));
                    for ev=1:length(eventsInSession)
                        powDataAllZscored(eventsInSession(ev),elec,adjust{1,windowSize},freq) = ((powDataAll(eventsInSession(ev),elec,adjust{1,windowSize},freq)) - meanDistrib) / stdDistrib;
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
            dashSplit = strfind(tline,'-');
            channelOne = tline(1:dashSplit-1);
            electrodes{iter,1} = channelOne;
            tline = fgetl(fid);
            iter = iter+1;
        end
        fclose(fid);
        
        % now, match each of these up with tal atlas, and categorize accordingly
        elementInfo = csv2cell(sprintf('/Volumes/Shares/FRNU/dataWorking/consolidationProject/data/%s/docs/element_info.csv',subjectNames{n}));
        atlasInfo = csv2cell(sprintf('/Volumes/Shares/FRNU/dataWorking/consolidationProject/data/%s/tal/atlas/atlas_monopolar_simple.csv',subjectNames{n}));
        
        MTL_elements = atlasInfo(union(union(find(contains(atlasInfo(:,6),'hippocampus','IgnoreCase',true)),find(contains(atlasInfo(:,4),'parahippocampal','IgnoreCase',true))),find(contains(atlasInfo(:,4),'entorhinal','IgnoreCase',true))),1);
        lateral_elements = atlasInfo(union(union(find(contains(atlasInfo(:,4),'inferior temporal gyrus','IgnoreCase',true)),find(contains(atlasInfo(:,4),'middle temporal gyrus','IgnoreCase',true))),find(contains(atlasInfo(:,4),'superior temporal gyrus','IgnoreCase',true))),1);
        ATL_elements = atlasInfo(find(contains(atlasInfo(:,7),'Anterior Temporal Lobe','IgnoreCase',true)),1);
        PTL_elements = setdiff(lateral_elements,ATL_elements);
        
        PPC_elements = atlasInfo(union(union(find(contains(atlasInfo(:,4),'supramarginal','IgnoreCase',true)),find(contains(atlasInfo(:,4),'parietal cortex','IgnoreCase',true))),find(contains(atlasInfo(:,4),'angular','IgnoreCase',true))),1);
        
        all_frontal_elements = atlasInfo(find(contains(atlasInfo(:,3),'frontal lobe','IgnoreCase',true)),1);
        dorsal_frontal_elements = atlasInfo(union(union(find(contains(atlasInfo(:,4),'frontal gyrus','IgnoreCase',true)),find(contains(atlasInfo(:,4),'pars opercularis','IgnoreCase',true))),find(contains(atlasInfo(:,4),'angular','IgnoreCase',true))),1);
        ventral_frontal_elements = atlasInfo(union(find(contains(atlasInfo(:,4),'rostral','IgnoreCase',true)),find(contains(atlasInfo(:,4),'orbital','IgnoreCase',true))),1);
        dorsal_frontal_elements = intersect(all_frontal_elements,dorsal_frontal_elements); %lets just make sure we got frontal lobe elecs
        ventral_frontal_elements = intersect(all_frontal_elements,ventral_frontal_elements); %lets just make sure we got frontal lobe elecs
        
        MTL_bipolars = [];
        ATL_bipolars = [];
        PTL_bipolars = [];
        ventral_frontal_bipolars = [];
        dorsal_frontal_bipolars = [];
        PPC_bipolars = [];
        % make sure this monopolar elec wasnt marked as bad
        for elec=1:size(electrodes,1)
            findDash = strfind(electrodes{elec},'-');
            elecOne = electrodes{elec};
            
            elecOne(find(isstrprop(elecOne,'digit'))) = [];
            
            % while we're here, lets load the bad_chans.mat and not allow any
            % chans marked in bad_chans.mat to pass this step
            find_bad_chan_one = length(find(contains(bad_chans,elecOne)));
            find_bad_chans_two = 0;
            if find_bad_chan_one + find_bad_chans_two > 0
                continue
            end
            
            MTL_elecs = find(contains(MTL_elements,electrodes{elec},'IgnoreCase',true));
            if size(MTL_elecs,1) ~= 0
                MTL_bipolars = cat(2,MTL_bipolars,elec);
                continue
            end
            ATL_elecs = find(contains(ATL_elements,electrodes{elec},'IgnoreCase',true));
            if size(ATL_elecs,1) ~= 0
                ATL_bipolars = cat(2,ATL_bipolars,elec);
                continue
            end
            PTL_elecs = find(contains(PTL_elements,electrodes{elec},'IgnoreCase',true));
            if size(PTL_elecs,1) ~= 0
                PTL_bipolars = cat(2,PTL_bipolars,elec);
                continue
            end
            ventral_frontal_elecs = find(contains(ventral_frontal_elements,electrodes{elec},'IgnoreCase',true));
            if size(ventral_frontal_elecs,1) ~=0
                ventral_frontal_bipolars = cat(2,ventral_frontal_bipolars,elec);
                continue
            end
            dorsal_frontal_elecs = find(contains(dorsal_frontal_elements,electrodes{elec},'IgnoreCase',true));
            if size(dorsal_frontal_elecs,1) ~=0
                dorsal_frontal_bipolars = cat(2,dorsal_frontal_bipolars,elec);
                continue
            end
            PPC_elecs = find(contains(PPC_elements,electrodes{elec},'IgnoreCase',true));
            if size(PPC_elecs,1) ~=0
                PPC_bipolars = cat(2,PPC_bipolars,elec);
                continue
            end
            
        end
        
        % aggregate data and save it out
        processed_power = cell(5,6);
        processed_power{5,1} = 'MTL'; processed_power{5,2} = 'ATL';
        processed_power{5,3} = 'PTL'; processed_power{5,4} = 'ventral_frontal';
        processed_power{5,5} = 'dorsal_frontal'; processed_power{5,6} = 'PPC';
        
        
        processed_power{1,1} = powDataAllZscored(:,MTL_bipolars,:,:);
        processed_power{1,2} = powDataAllZscored(:,ATL_bipolars,:,:);
        processed_power{1,3} = powDataAllZscored(:,PTL_bipolars,:,:);
        processed_power{1,4} = powDataAllZscored(:,ventral_frontal_bipolars,:,:);
        processed_power{1,5} = powDataAllZscored(:,dorsal_frontal_bipolars,:,:);
        processed_power{1,6} = powDataAllZscored(:,PPC_bipolars,:,:);
        processed_power{1,7} = 'Z-Scored To Prestim, -1000 to 500';
        
        
        
        if ~exist(sprintf('/Volumes/FRNU_EXT/memoryStrupt_Data/preprocessed/processedPower_LowAndHigh_monopolar_window%d/%s/',windowSize,subjectNames{n})) %#ok<EXIST>
            mkdir(sprintf('/Volumes/FRNU_EXT/memoryStrupt_Data/preprocessed/processedPower_LowAndHigh_monopolar_window%d/%s/',windowSize,subjectNames{n}));
        end
        
        save(sprintf('/Volumes/FRNU_EXT/memoryStrupt_Data/preprocessed/processedPower_LowAndHigh_monopolar_window%d/%s/allPowData.mat',windowSize,subjectNames{n}),'processed_power','evsOfInterest','-v7.3')
        
        fprintf('\ncompleted %s\n',subjectNames{n})
        
        
    end
end

