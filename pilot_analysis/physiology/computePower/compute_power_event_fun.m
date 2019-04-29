function compute_power_event_fun(event,event_num,pars,zz)
    %Compue the power for a given event
    
    eegRootDir = pars.eegRootDir;
    electrodes = pars.electrodes;
    freqs = pars.freqs;
    powDir = pars.powDir;
    msbuffer = pars.msbuffer;
    resamprate = pars.resamprate;
    bp_flag = pars.bp_flag;
    
    points_per_win = pars.points_per_win;
    points_per_slide = pars.points_per_slide;
    width = pars.width;
    
    if strcmp(event.type,'STUDY_PAIR'),
       ms_before = pars.ms_before(zz);
       ms_after = pars.ms_after(zz);
    else
       ms_before = pars.ms_before_ret(zz);
       ms_after = pars.ms_after_ret(zz);
    end
    
    %downsample_points = pars.downsample_points;
    mirror = ceil((pars.msmirror/1000)*pars.SR);
    total_points = round((ms_before+ms_after)/1000*pars.SR);
    if mirror > total_points, mirror = total_points-1; end
    
    fname = [powDir '/event' num2str(event_num) 'powData.mat'];%THIS IS THE ONE THAT SHOULD BE USED
    fname_lfp = [powDir '/event' num2str(event_num) 'LFP.mat'];
    
    if ~exist(fname,'file')

        %dataTmp=get_eeg_palram3(eegRootDir, event,ms_before+msbuffer,ms_after+msbuffer,bp_flag,0,[58 62],'stop',4,[],electrodes,pars)';
        dataTmp=get_eeg_palram3(eegRootDir, event,ms_before+msbuffer,ms_after+msbuffer,bp_flag,0,[58 62],'stop',4,[],electrodes,pars)';
        dataMir = [fliplr(dataTmp(:,2:mirror+1)) dataTmp fliplr(dataTmp(:,end-mirror:end-1))];
        
        for i = 1:size(dataMir,1),   
             [~, wavPow(i,:,:)] = multiphasevec2(freqs,dataMir(i,:),pars.SR,width);
        end
        %wavPow is electrodes X freqs X time
        
        
        wavPow = wavPow(:,:,msbuffer+1+mirror:end-msbuffer-mirror);
        wavPow = log10(wavPow);
        
        %powData is electrodes X freqs X time
        powData = mean(wavPow,3);
        save(fname,'powData')
        
        if pars.getERP
            save_lfp = downsample(dataTmp',pars.ERP_dsample)';
            save(fname_lfp,'save_lfp')
        end
    end
    
end