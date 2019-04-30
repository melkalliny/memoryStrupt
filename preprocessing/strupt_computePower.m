

pars = consolidation_setParams();
expType = 1;

subjects = [71];
subjectNames = {'NIH071'};

npats = [1];
% needs to be edited so that it uses the RT to pull out data centered
% around response times
for windowSize=1:3
    if windowSize==1
        pars.points_per_win = 500;
        pars.points_per_win_ret = 100;
        pars.points_per_slide = 50;   %prevent bug propagation
    elseif windowSize==2
        pars.points_per_win = 1000;
        pars.points_per_win_ret = 100;
        pars.points_per_slide = 100;   %prevent bug propagation
    elseif windowSize==3
        pars.points_per_win = 2000;
        pars.points_per_win_ret = 100;
        pars.points_per_slide = 200;   %prevent bug propagation
    end
    
    pars.dirWriteOut = sprintf('%s/%d/',pars.dirWriteOut,windowSize);
    
    for n = 1
        
        tic;
        subj = subjectNames{n};
        subjNum = subjects(n);
        
        % use bookkeeping file to choose right file, load data, and select output location
        %bookkeeping = readtable('/Volumes/FRNU_SVN/semanticSpan/data/bookkeeping.csv');
        %bookkeeping = table2cell(bookkeeping);
        
        pat_dir_out = [pars.dirWriteOut '/' subj];
        
        
        powDir = [pat_dir_out '/powDir'];
        pars.subjectWriteOut = powDir;
        
        pars.SR = 1000;
        pars.subj = subj;
        freqs = pars.freqs;
        
        % load events
        load(sprintf('%s/%s/behavioral/memoryRepeat/events.mat',pars.froot,subj));
        
        [events(find(strcmp({events.type},'STUDY_PAIR'))).type] = deal('STUDY_PAIR_START');
        [events(find(strcmp({events.type},'TEST_PROBE'))).type] = deal('PROBE_START');
        evsPull = events((union(find(strcmp({events.type},'STUDY_PAIR_START')),find(strcmp({events.type},'PROBE_START')))));
    
        num_events = length(evsPull);
        
        parfor e = 1:num_events
            ev = evsPull(e);
            toShift = 0;
            ev.eegfile = strrep(ev.eegfile,'noreref','processed');
            compute_power_event_fun_noMir(ev,e,pars,expType);
        end
        
        
        clear powDataAll avgPowData
        for e = 1:num_events
            ev = evsPull(e);
            fname = [powDir '/event' num2str(e) 'powData.mat'];
            load(fname,'powData')
            
            if pars.getERP
                fname_lfp = [powDir '/event' num2str(e) 'LFP.mat'];
                if exist(fname_lfp,'file')
                    load(fname_lfp,'save_lfp');
                else
                    warning('%s does not exist... filling w/ NaNs instead',fname_lfp)
                    save_lfp = NaN(size(save_lfp));
                end
            end
            
            if ~exist('powDataAll','var')
                powDataAll = single(nan([num_events, size(powData)]));
                if pars.getERP
                    lfpDataAll = single(nan([num_events, size(save_lfp)]));
                end
                
            end
            powDataAll(e,:,:,:) = powData;
            %powDataAll is num_events X electrode X freq X time
            
            if pars.getERP
                buff = size(lfpDataAll,3)-size(save_lfp,2);
                if buff>0
                    lfpDataAll(e,:,:) = [save_lfp NaN(size(save_lfp,1),buff)];
                else
                    lfpDataAll(e,:,:) = save_lfp(:,1:size(lfpDataAll,3));
                end
                %delete(fname_lfp)
            end
            % delete(fname)
            
        end
        
        
        save([powDir '/allPowData.mat'],'powDataAll','-v7.3')
        save([powDir '/pars'],'pars');
        save([powDir '/events.mat'],'evsPull');
        
        if pars.getERP
            save([powDir '/LFPData.mat'],'lfpDataAll','-v7.3')
            weightsReject = [1 1];
            [iChanKeep, iEvKeep, strOut] = cleanEEGevents(permute(lfpDataAll,[2 1 3]),weightsReject); % input: chan x events x time
            save([powDir '/rejectionInfo.mat'],'iChanKeep', 'iEvKeep', 'strOut', 'weightsReject');
        end
        
        fprintf('%s complete!\n',subj);
        
        toc;
        
    end
    
end