

pars = consolidation_setParams_ripples();
expType = 1;

subjects = [49 50 51 52 53 55 58];
subjectNames = {'NIH049', 'NIH050', 'NIH051', 'NIH052', 'NIH053', 'NIH055', 'NIH058'};
npats = [1 2 3 4 5 6 7];


% needs to be edited so that it uses the RT to pull out data centered
% around response times

for n = [1:7]
    tic;
    subj = subjectNames{n};
    subjNum = subjects(n);
    
    % use bookkeeping file to choose right file, load data, and select output location
    %bookkeeping = readtable('/Volumes/FRNU_SVN/semanticSpan/data/bookkeeping.csv');
    %bookkeeping = table2cell(bookkeeping);
    
    pat_dir_out = [pars.dirWriteOut '/' subj];
    
    
    rippleDir = [pat_dir_out '/rippleDir'];
    pars.subjectWriteOut = rippleDir;
    
    pars.SR = 1000;
    pars.subj = subj;
    freqs = pars.freqs;
    
    % load events
    load(sprintf('/Volumes/Shares/FRNU/data/eeg/%s/behavioral/palRam/events.mat',subj));
    evsPull = events((union(find(strcmp({events.type},'STUDY_PAIR_START')),find(strcmp({events.type},'PROBE_START')))));
    num_events = length(evsPull);
    
    eegRootDir = pars.eegRootDir;
    elecDir = sprintf('%s/%s/docs/leads_bp.txt',pars.eegRootDir,pars.subj);
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
    
    dataTmp = NaN(6000,length(electrodes),num_events);
    
    parfor e = 1:num_events
        ev = evsPull(e);
        toShift = 0;
        ev.eegfile = strrep(ev.eegfile,'noreref','processed');
        % compute_power_event_fun_noMir(ev,e,pars,expType);
        
        
        freqs = pars.freqs;
        msbuffer = pars.msbuffer;
        bp_flag = pars.bp_flag;
        adjust_points = round(1000/pars.SR);
        points_per_win = round(pars.points_per_win/adjust_points);
        points_per_slide = round(pars.points_per_slide/adjust_points);
        width = pars.width;
        
        
        ms_before = pars.ms_before(expType);
        ms_after = pars.ms_after(expType);
        
        
        mirror = ceil((pars.msmirror/1000)*pars.SR);
        total_points = round((ms_before+ms_after)/1000*pars.SR);
        if mirror > total_points, mirror = total_points-1; end
        
        try %--in case new range is outside of session
            [dataTmp_ev, has_error] = get_eeg_consolidation_ripples(eegRootDir,ev,ms_before+msbuffer,ms_after+msbuffer,bp_flag,0,[58 62],'stop',4,[],electrodes,pars);
        catch
            fprintf('\nfailed\n')
            continue
        end
        

        fprintf('\nevent%d done\n',e)
        
        dataTmp(:,:,e) = dataTmp_ev;
    end
    
    
    
    for chan=1:size(dataTmp,2)
        for ev=1:size(dataTmp,3)
            feeg = dataTmp(:,chan,ev);
            temp_eeg = hilbert(feeg')';
        henv = abs(temp_eeg);
        dataTmp(:,chan,ev) = henv;
        end
    end
    
    uniqueSessions = unique([evsPull.session]);
    for session=1:length(uniqueSessions)
        events = find([evsPull.session] == uniqueSessions(session));
        for chan=1:size(dataTmp,2)
            tempDistrib = reshape(dataTmp(:,chan,events),1,size(dataTmp,1)*length(events));
            tempMean = nanmean(tempDistrib);
            tempSte = nanstd(tempDistrib);
            for ev=1:length(events)
                dataTmp(:,chan,events(ev)) = (dataTmp(:,chan,events(ev)) - tempMean) / tempSte;
            end
        end
    end
    
    
    ripple_evs = cell(1,num_events);
    
    rippleDir = pars.subjectWriteOut;
    % now, lets find ripples and assign them into a cell
    for ev=1:size(dataTmp,3)
        feeg = squeeze(dataTmp(:,:,ev));
        %feeg = dataTmp';
        % zenv = zscore(henv',[ev.session],0);
        %save(fullfile(d.out_dir,'henv.mat'),'henv','zenv','-v7.3')
        
        % ripple boolean
        %zenv = zscore(henv,0,2);
        zenv = feeg;
        thresh = 3; minDur = 20;
        zenv(zenv<thresh)  = 0;
        zenv(zenv>=thresh) = 1;
        zenv = zenv';
        
        % initialize rips
        rips  = zeros(size(zenv),'single');
        rips_sta  = zeros(size(zenv),'single');
        %t_ind = -1*p.off_ms:(-1*p.off_ms+2000);
        for i = 1:size(zenv,1)
            series = zenv(i,:);
            rip(i).mid = NaN;
            rip(i).nr  = NaN;
            rip(i).dur = NaN;
            rip(i).isi = NaN;
            if sum(series)<2; continue; end
            
            % ripple stats
            [sta{i},en{i}] = consecutive_ones(series');
            dur = en{i}-sta{i};
            mid = round(sta{i}+dur/2);
            
            % number duration and isi of ripples
            rip(i).mid = mid;
            rip(i).nr  = sum(dur>minDur);
            rip(i).dur = dur(dur>minDur);
            rip(i).isi = diff(mid(dur>minDur));
            
            % ripple mid
            rips(i,mid) = 1;
            rips_sta(i,sta{i})= 1;
        end
        m.rips_sta = rips_sta;
        m.rips     = rips;
        m.rip      = rip;
        
        
        ripple_evs{1,ev} = m.rips;
        fprintf('\n%d done\n',ev)
    end
    
    if ~exist(rippleDir); mkdir(rippleDir); end
    save([rippleDir '/allRipples.mat'],'ripple_evs','events','-v7.3')
    
    fprintf('%s complete!\n',subj);
    
    toc;
    
end