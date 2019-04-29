function has_error = compute_power_event_fun_noMir(event,event_num,pars,expType)
    %Compue the power for a given event
    has_error = false; 
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
   
    freqs = pars.freqs;
    powDir = pars.subjectWriteOut;
    msbuffer = pars.msbuffer;
    bp_flag = pars.bp_flag;
    adjust_points = round(1000/pars.SR);
    points_per_win = round(pars.points_per_win/adjust_points);
    points_per_slide = round(pars.points_per_slide/adjust_points); 
    width = pars.width;
    
    if isfield(pars,'JUST_LFP') && pars.JUST_LFP
    	JUST_LFP = true();
    else
    	JUST_LFP = false();
    end	   
        
    if strcmpi(event.type,'RESPONSE')
    ms_before = pars.ms_before_ret(expType);
    ms_after = pars.ms_after_ret(expType);        
    else
    ms_before = pars.ms_before(expType);
    ms_after = pars.ms_after(expType);
    end
    
    
    mirror = ceil((pars.msmirror/1000)*pars.SR);
    total_points = round((ms_before+ms_after)/1000*pars.SR);
    if mirror > total_points, mirror = total_points-1; end
    
    fname = [powDir '/event' num2str(event_num) 'powData.mat'];%THIS IS THE ONE THAT SHOULD BE USED
    fname_lfp = [powDir '/event' num2str(event_num) 'LFP.mat'];
   
    if ~exist(fname,'file') || exist(fname,'file')

        
        try %--in case new range is outside of session
            [dataTmp, has_error] = get_eeg_consolidation(eegRootDir,event,ms_before+msbuffer,ms_after+msbuffer,bp_flag,0,[58 62],'stop',4,[],electrodes,pars);
        catch
            keyboard
        end
        
        dataTmp = dataTmp';
        dataMir = [fliplr(dataTmp(:,2:mirror+1)) dataTmp fliplr(dataTmp(:,end-mirror:end-1))];
        
        if JUST_LFP
        	save_lfp = dataTmp;
            save(fname_lfp,'save_lfp')
            return
        end
        
       
        for i = 1:size(dataMir,1)
            [~, wavPow(i,:,:)] = multiphasevec3(freqs,dataMir(i,:),pars.SR,width,1);
        end
        
        %wavPow is electrodes X freqs X time
        
        wavPow = wavPow(:,:,msbuffer+1+mirror:end-msbuffer-mirror); %- remove mirror
        wavPow = log10(wavPow);
        
        
        powData = windowed_average_with_freq(wavPow,points_per_win,points_per_slide);
        %powData = wavPow;
        
        dashes = find(fname=='/');
        if ~exist(fname(1:dashes(end))); mkdir(fname(1:dashes(end))); end
        save(fname,'powData','-v7.3')
        
        
        if pars.getERP
            
            save_lfp = dataTmp;
            save(fname_lfp,'save_lfp','-v7.3')
           
        end
    end

end