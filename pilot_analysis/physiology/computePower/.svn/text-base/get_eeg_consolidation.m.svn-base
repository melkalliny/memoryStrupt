function [EEG, has_error] = get_eeg_consolidation(eegRootDir, event,ms_before,ms_after,bpflag,BufferMS,filtfreq,filttype,filtorder,resamprate,electrodes,pars)
    %This function retrieves the EEG data associated with an event
    %pat_s: the patient identifier: NIH006, TJ027, etc
    %event: an event struct
    %ms_before: the number of ms to pull before the event
    %ms_after: the number of ms to pull after the event
    %bp_flag: 1 = bipolar montage, 0 = single channel data
       
    if ~exist('BufferMS','var')
        BufferMS = [];
    end
    if ~exist('filtfreq','var')
        filtfreq = [];
    end
    if ~exist('filttype','var')
        filttype = [];
    end
    if ~exist('filtorder','var')
        filtorder = [];
    end
        
    has_error = false;
    num_electrodes = size(electrodes,1);
    num_bp_electrodes = size(electrodes,1); 
    
    
    totalpoints = ms_before+ms_after;
    err_num = 0;
    if bpflag
        EEG = nan(ceil(totalpoints/1000*ceil(pars.SR)),num_electrodes);
        elec_count = 1;
        for bps = 1:num_electrodes
            try
                
                
                % see if monopolar or bipolar, and process appropriately
                if strcmpi(pars.BipolarOrMonopolar,'bipolar')
                    dashSplit = strfind(electrodes{bps},'-');
                    channelOne = electrodes{bps}(1:dashSplit-1);
                    channelTwo = electrodes{bps}(dashSplit+1:end);
                    EEG1 = gete_ms_mod(eegRootDir,channelOne,event,totalpoints,-ms_before,BufferMS);%,filtfreq,filttype,filtorder);
                    EEG2 = gete_ms_mod(eegRootDir,channelTwo,event,totalpoints,-ms_before,BufferMS);%,filtfreq,filttype,filtorder);
                    if isempty(filtorder) || filtorder == 0
                        EEG(:,bps) = EEG1;
                    else
                        EEG(:,bps) = buttfilt(EEG1-EEG2,filtfreq,pars.SR,filttype,filtorder);
                    end
                    
                else
                    dashSplit = strfind(electrodes{bps},'-');
                    channelOne = electrodes{bps}(1:dashSplit-1);
                    EEG1 = gete_ms_mod(eegRootDir,channelOne,event,totalpoints,-ms_before,BufferMS);%,filtfreq,filttype,filtorder);
                    EEG(:,bps) = buttfilt(EEG1,filtfreq,pars.SR,filttype,filtorder);
                end
                 
                % EEG(:,elec_count) = buttfilt(EEG1-EEG2,filtfreq,pars.SR,filttype,filtorder);
                
            catch err
                %display(['bad channel ' num2str(elec_count)])
                continue
            end
            elec_count = elec_count + 1;
        end
    end
end



















%                 if ~exist('filttype','var') || isempty(filttype); filttype = 'stop'; end
%                 if ~exist('filtorder','var') || isempty(filtorder); filtorder = 1; end
%                 sampleFreq = 1000;
%                 if( butterFiltSamplerate ~= sampleFreq )
%                     [fileEEG(fullEEGmask,1:readDuration),butterFilt] = buttfilt(fileEEG(fullEEGmask,1:readDuration),filtfreq,sampleFreq,filttype,filtorder);
%                     butterFiltSamplerate = sampleFreq;
%                 else
%                     fileEEG(fullEEGmask,1:readDuration) = buttfilt(fileEEG(fullEEGmask,1:readDuration),butterFilt);
%                 end
%                 
%                 for i = 1:length(shortEEGind)
%                     ind = shortEEGind(i);
%                     % no need to check that butterFilt is the correct sampling rate, since that was just done above
%                     fileEEG(ind,1:fileEEGlen(ind)) = buttfilt(fileEEG(ind,1:fileEEGlen(ind)),butterFilt);
%                 end
% old code here
%                 EEG1 = gete_ms(eegRootDir,electrodes(elec_count,1),event,totalpoints,-ms_before,BufferMS);%,filtfreq,filttype,filtorder);
%                 % subtract the common if we need to do that
%                 EEG2 = gete_ms(eegRootDir,electrodes(elec_count,2),event,totalpoints,-ms_before,BufferMS);%,filtfreq,filttype,filtorder);
